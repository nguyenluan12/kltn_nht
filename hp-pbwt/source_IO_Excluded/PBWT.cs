using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HP_PBWT
{
    public class PBWT
    {
        public class PDA
        {
            public int[] pArr;

            public int zeroCnt = 0;
            public int[] mLens;

            //public int nHap;
            public PDA(int nHapTotal)
            {
                //nHap = nHapTotal;
                pArr = new int[nHapTotal];

                mLens = new int[nHapTotal];

            }

        }

        public class Pal
        {
            public class UpLowRange
            {
                public int UpRange;
                public int LowerRange;
                public UpLowRange(int up, int low)
                {
                    UpRange = up;
                    LowerRange = low;
                }
            }

            static public PDA newPDA;
            static public PDA oldPDA;
            static PDA temPDA;
            static int[] psHolder;
            static int[] offsetsHolder;
            static int[] minDZ;
            static bool[] minDZ_Ready;


            static int blockSize = 0;

            static ParallelOptions pOp = new ParallelOptions();

            static int msD_Wait = 10;

            static int blockHasOneSignal = -10;

            public Pal()
            {
                Console.WriteLine(DateTime.Now + " Initializing HP-PBWT " + Program.nHap + " Samples (#haplotypes) " + Program.nThread + " Threads...");


                newPDA = new PDA(Program.nHap);
                oldPDA = new PDA(Program.nHap);

                psHolder = new int[Program.nHap];
                offsetsHolder = new int[Program.nThread];
                minDZ = new int[Program.nThread];
                minDZ_Ready = new bool[Program.nThread];


                pOp.MaxDegreeOfParallelism = Program.nThread;
                blockSize = Convert.ToInt32(Math.Ceiling(Convert.ToDouble(Program.nHap) / Convert.ToDouble(Program.nThread)));

                Console.WriteLine(DateTime.Now + " HP-PBWT Initialized.");
            }

            //need work
            public void coreP_Arr(int siteToSort)
            {
                BitArray oneSite = Program.panel[siteToSort];
                #region step 1 local sum
                Parallel.For(0, Program.nThread, (i) =>
                {

                    int s = i * blockSize;
                    if (s >= Program.nHap)
                    { return; }
                    int e = i * blockSize + blockSize;
                    if (e > Program.nHap)
                    {
                        e = Program.nHap;
                    }

                    if (oneSite[oldPDA.pArr[s]] == true)
                    {
                        psHolder[s] = 1;
                    }
                    else
                    {
                        psHolder[s] = 0;
                    }

                    for (int k = s + 1; k < e; k++)
                    {
                        //psHolder[k] = psHolder[k - 1] + site[oldPDA.pArr[k]];
                        if (oneSite[oldPDA.pArr[k]] == true)
                        {
                            psHolder[k] = psHolder[k - 1] + 1;
                        }
                        else
                        {
                            psHolder[k] = psHolder[k - 1];
                        }
                    }

                    offsetsHolder[i] = psHolder[e - 1];

                });
                #endregion

                #region step 2 seq offset handling
                for (int i = 1; i < Program.nThread; i++)
                {
                    offsetsHolder[i] = offsetsHolder[i - 1] + offsetsHolder[i];
                }

                #endregion

                #region step 3 apply offsets
                Parallel.For(1, Program.nThread, (i) =>
                {
                    int s = i * blockSize;
                    if (s >= Program.nHap)
                    { return; }
                    int e = i * blockSize + blockSize;
                    if (e > Program.nHap)
                    {
                        e = Program.nHap;
                    }

                    for (int k = s; k < e; k++)
                    {
                        psHolder[k] = psHolder[k] + offsetsHolder[i - 1];
                    }
                });
                #endregion

                #region step 4 pps -> index settle to new p arr 
                newPDA.zeroCnt = Program.nHap - psHolder[Program.nHap - 1];
                int oneOff = newPDA.zeroCnt - 1;
                Parallel.For(0, Program.nHap, (i) =>
                {
                    if (oneSite[oldPDA.pArr[i]] == true)
                    {
                        newPDA.pArr[psHolder[i] + oneOff] = oldPDA.pArr[i];
                    }
                    else
                    {
                        newPDA.pArr[i - psHolder[i]] = oldPDA.pArr[i];
                    }
                });
                #endregion

            }

            public void coreD_Arr(int siteToSort)
            {

                BitArray oneSite = Program.panel[siteToSort];

                //for (int i = 0; i < nThread; i++)
                Parallel.For(0, Program.nThread, (i) =>
                {
                    int s = i * blockSize;
                    if (s >= Program.nHap)
                    { return; }

                    int e = i * blockSize + blockSize;

                    if (e >= Program.nHap)
                    {
                        e = Program.nHap;
                    }

                    int prvLowM_Zero = -1;
                    int prvLowM_One = -1;
                    int hID = oldPDA.pArr[s];


                    #region first block
                    if (s == 0)
                    {
                        prvLowM_Zero = -1;
                        prvLowM_One = -1;
                        hID = oldPDA.pArr[s];

                        for (int k = s; k < e; k++)
                        {

                            hID = oldPDA.pArr[k];

                            //incoming is 0
                            if (oneSite[hID] == false)
                            {
                                newPDA.mLens[hID] = Math.Min(oldPDA.mLens[hID], prvLowM_One) + 1;
                                prvLowM_One = int.MaxValue;
                                prvLowM_Zero = Math.Min(prvLowM_Zero, oldPDA.mLens[hID]);
                            }
                            else//incoming is 1
                            {
                                newPDA.mLens[hID] = Math.Min(oldPDA.mLens[hID], prvLowM_Zero) + 1;
                                prvLowM_Zero = int.MaxValue;
                                prvLowM_One = Math.Min(prvLowM_One, oldPDA.mLens[hID]);
                            }
                        }


                    }
                    #endregion

                    #region other blocks
                    else
                    {

                        //look up to set prvLowM_One and prvLowM_Zero
                        prvLowM_Zero = -1;
                        prvLowM_One = -1;

                        bool minZeroSearch = true;
                        bool minOneSearch = true;

                        //special cases:
                        //A there is no 0 in upper blocks
                        //B there is no 1 in upper blocks
                        if (psHolder[s - 1] == s)//A there is no 0 in upper blocks
                        {
                            prvLowM_Zero = -1;
                            minZeroSearch = false;
                        }

                        if (psHolder[s - 1] == 0) //B there is no 1 in upper blocks
                        {
                            prvLowM_One = -1;
                            minOneSearch = false;
                        }
                        // oneSite[oldPDA.pArr[seakIndex]

                        //minZero and minOne should only happen once
                        if (oneSite[oldPDA.pArr[s - 1]] == false)
                        {
                            prvLowM_One = int.MaxValue;
                            minOneSearch = false;
                        }
                        else
                        {
                            prvLowM_Zero = int.MaxValue;
                            minZeroSearch = false;
                        }

                        
                        //regular search


                        if (minZeroSearch)
                        {

                            int seekIndex = s - 1;
                            //locate first upper zero
                            //seekHID = oldPDA.pArr[seekIndex];
                            while (seekIndex >= 0 && oneSite[oldPDA.pArr[seekIndex]] != false)
                            {
                                seekIndex--;
                            }
                            if (seekIndex >= 0)
                            {//ToDo: may need to change
                                prvLowM_Zero = oldPDA.mLens[oldPDA.pArr[seekIndex]];
                                while (seekIndex >= 0 && oneSite[oldPDA.pArr[seekIndex]] == false)
                                {
                                    prvLowM_Zero = Math.Min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_Zero);
                                    seekIndex--;
                                }

                                if (seekIndex == -1)
                                {
                                    prvLowM_Zero = -1;
                                }
                            }

                        }


                        if (minOneSearch)
                        {//go through a block of one find the min
                            int seekIndex = s - 1;
                            //locate first upper one
                            while (seekIndex >= 0 && oneSite[oldPDA.pArr[seekIndex]] == false)
                            {
                                seekIndex--;
                            }
                            if (seekIndex >= 0)
                            {
                                prvLowM_One = oldPDA.mLens[oldPDA.pArr[seekIndex]];
                                while (seekIndex >= 0 && oneSite[oldPDA.pArr[seekIndex]] != false)
                                {
                                    prvLowM_One = Math.Min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_One);
                                    seekIndex--;
                                }

                                if (seekIndex == -1)
                                {
                                    prvLowM_One = -1;
                                }
                            }
                        }



                        for (int k = s; k < e; k++)
                        {

                            hID = oldPDA.pArr[k];

                            //incoming is 0
                            if (oneSite[hID] == false)
                            {
                                newPDA.mLens[hID] = Math.Min(oldPDA.mLens[hID], prvLowM_One) + 1;
                                prvLowM_One = int.MaxValue;
                                prvLowM_Zero = Math.Min(prvLowM_Zero, oldPDA.mLens[hID]);


                            }
                            else//incoming is 1
                            {
                                newPDA.mLens[hID] = Math.Min(oldPDA.mLens[hID], prvLowM_Zero) + 1;
                                prvLowM_Zero = int.MaxValue;
                                prvLowM_One = Math.Min(prvLowM_One, oldPDA.mLens[hID]);
                            }
                        }
                    }
                    #endregion

                    //}
                });


            }


            public void ResetMinDZ()
            {
                for (int i = 0; i < Program.nThread; i++)
                {
                    minDZ[i] = int.MaxValue;
                    minDZ_Ready[i] = false;
                }
            }

            public void coreD_Arr_s3b(int siteToSort)
            {
                ResetMinDZ();
                BitArray oneSite = Program.panel[siteToSort];

                //for (int i = 0; i < nThread; i++)
                Parallel.For(0, Program.nThread, (i) =>
                {
                    int s = i * blockSize;
                    if (s >= Program.nHap)
                    { return; }

                    int e = i * blockSize + blockSize;

                    if (e >= Program.nHap)
                    {
                        e = Program.nHap;
                    }

                    int prvLowM_Zero = -1;
                    int prvLowM_One = -1;
                    int hID = oldPDA.pArr[s];

                    //minDZ[i] = int.MaxValue;


                    #region first block
                    if (s == 0)
                    {
                        if (psHolder[e - 1] != 0)
                        {
                            minDZ[i] = blockHasOneSignal;
                            minDZ_Ready[i] = true;
                        }

                        prvLowM_Zero = -1;
                        prvLowM_One = -1;
                        hID = oldPDA.pArr[s];

                        for (int k = s; k < e; k++)
                        {

                            hID = oldPDA.pArr[k];

                            //incoming is 0

                            if (oneSite[hID] == false)
                            {
                                newPDA.mLens[hID] = Math.Min(oldPDA.mLens[hID], prvLowM_One) + 1;
                                prvLowM_One = int.MaxValue;
                                prvLowM_Zero = Math.Min(prvLowM_Zero, oldPDA.mLens[hID]);

                                //minDZ[i] = Math.Min(newPDA.mLens[hID], minDZ[i]);
                            }
                            else//incoming is 1
                            {
                                newPDA.mLens[hID] = Math.Min(oldPDA.mLens[hID], prvLowM_Zero) + 1;
                                prvLowM_Zero = int.MaxValue;
                                prvLowM_One = Math.Min(prvLowM_One, oldPDA.mLens[hID]);

                                //minDZ[i] = -10;
                            }
                        }


                    }
                    #endregion

                    #region other blocks
                    else
                    {

                        if (psHolder[e - 1] != psHolder[s - 1])
                        {
                            minDZ[i] = blockHasOneSignal;
                            minDZ_Ready[i] = true;
                        }
                        else
                        {
                            //need to put this elsewhere
                            for (int j = s; j < e; j++)
                            {
                                minDZ[i] = Math.Min(minDZ[i], oldPDA.mLens[oldPDA.pArr[j]]);

                            }
                            minDZ_Ready[i] = true;
                        }

                        //look up to set prvLowM_One and prvLowM_Zero
                        prvLowM_Zero = -1;
                        prvLowM_One = -1;

                        bool minZeroSearch = true;
                        bool minOneSearch = true;

                        //special cases:
                        //A there is no 0 in upper blocks
                        //B there is no 1 in upper blocks
                        if (psHolder[s - 1] == s)//A there is no 0 in upper blocks
                        {
                            prvLowM_Zero = -1;
                            minZeroSearch = false;
                        }

                        if (psHolder[s - 1] == 0) //B there is no 1 in upper blocks
                        {
                            prvLowM_One = -1;
                            minOneSearch = false;
                        }

                        //minZero and minOne should only happen once

                        if (oneSite[oldPDA.pArr[s - 1]] == false)
                        {
                            //prvLowM_One = int.MaxValue;
                            minOneSearch = false;
                        }
                        else
                        {
                            //prvLowM_Zero = int.MaxValue;
                            minZeroSearch = false;
                        }

                        // search

                        if (minZeroSearch)
                        {

                            int seekIndex = s - 1;
                            // Benny, do work here!

                            //skip method to compute prvLowM_Zero

                            int blk_ID = i - 1;
                            //int bs = blk_ID * blockSize;

                            //int be = blk_ID * blockSize + blockSize;

                            while (blk_ID > 0)
                            {
                                if ((psHolder[blk_ID * blockSize - 1] != psHolder[(blk_ID + 1) * blockSize - 1]))
                                {
                                    break;
                                }

                                blk_ID--;
                            }
                            seekIndex = (blk_ID + 1) * blockSize - 1;
                            prvLowM_Zero = oldPDA.mLens[oldPDA.pArr[seekIndex]];

                            while (seekIndex >= 0 && oneSite[oldPDA.pArr[seekIndex]] == false)
                            {
                                prvLowM_Zero = Math.Min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_Zero);
                                seekIndex--;
                            }


                            if (seekIndex == -1)
                            {
                                prvLowM_Zero = -1;
                            }


                            for (int b = blk_ID + 1; b < i; b++)
                            {
                                while (minDZ_Ready[b] == false)
                                {
                                    continue;
                                }
                                prvLowM_Zero = Math.Min(minDZ[b], prvLowM_Zero);
                            }



                        }


                        if (minOneSearch)
                        {
                            int seekIndex = s - 1;
                            //locate first upper one

                            if (seekIndex >= 0)
                            {
                                prvLowM_One = oldPDA.mLens[oldPDA.pArr[seekIndex]];
                                while (seekIndex >= 0 && oneSite[oldPDA.pArr[seekIndex]] != false)
                                {
                                    prvLowM_One = Math.Min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_One);
                                    seekIndex--;
                                }

                                if (seekIndex == -1)
                                {
                                    prvLowM_One = -1;
                                }
                            }
                        }

                        if (oneSite[oldPDA.pArr[s - 1]] == false)
                        {
                            prvLowM_One = int.MaxValue;
                        }
                        else
                        {
                            prvLowM_Zero = int.MaxValue;
                        }


                        for (int k = s; k < e; k++)
                        {

                            hID = oldPDA.pArr[k];

                            //incoming is 0
                            if (oneSite[hID] == false)
                            {
                                newPDA.mLens[hID] = Math.Min(oldPDA.mLens[hID], prvLowM_One) + 1;
                                prvLowM_One = int.MaxValue;
                                prvLowM_Zero = Math.Min(prvLowM_Zero, oldPDA.mLens[hID]);


                            }
                            else//incoming is 1
                            {
                                newPDA.mLens[hID] = Math.Min(oldPDA.mLens[hID], prvLowM_Zero) + 1;
                                prvLowM_Zero = int.MaxValue;
                                prvLowM_One = Math.Min(prvLowM_One, oldPDA.mLens[hID]);
                            }
                        }
                    }
                    #endregion

                    //}
                });


            }




            public void ReportLongMatches_SL(int siteIndex)
            {

                Parallel.For(0, Program.nHap - 1, (h) =>
                {
                    int minVal = int.MaxValue;
                    for (int i = h + 1; i < Program.nHap; i++)
                    {


                        if (oldPDA.mLens[oldPDA.pArr[i]] < Program.LLM_Len)
                        { break; }

                        minVal = Math.Min(minVal, oldPDA.mLens[oldPDA.pArr[i]]);


                        if (Program.panel[siteIndex + 1][oldPDA.pArr[h]] != Program.panel[siteIndex + 1][oldPDA.pArr[i]])
                        {
                            //report
                            //Program.BW.Add(oldPDA.pArr[h] + "\t" + oldPDA.pArr[i] + "\t" + siteIndex + "\t" + minVal);
                            Program.Sums[oldPDA.pArr[h]] += minVal;
                            Program.Sums[oldPDA.pArr[h]] = Program.Sums[oldPDA.pArr[h]] % Program.nSite;
                        }
                    }
                }
                );


            }

            public void ReportLongMatches_SL_Tail(int siteIndex)
            {

                Parallel.For(0, Program.nHap - 1, (h) =>
                {
                    int minVal = int.MaxValue;
                    for (int i = h + 1; i < Program.nHap; i++)
                    {
                        if (oldPDA.mLens[oldPDA.pArr[i]] < Program.LLM_Len)
                        {
                            break;
                        }

                        minVal = Math.Min(minVal, oldPDA.mLens[oldPDA.pArr[i]]);
                        //report
                        Program.Sums[oldPDA.pArr[h]] += minVal;
                        Program.Sums[oldPDA.pArr[h]] = Program.Sums[oldPDA.pArr[h]] % Program.nSite;
                    }
                });

            }



            //create new P and D
            public void OneSort(int siteToSort)
            {
                coreP_Arr(siteToSort);

                //coreD_Arr(siteToSort);
                coreD_Arr_s3b(siteToSort);

                temPDA = oldPDA;
                oldPDA = newPDA;
                newPDA = temPDA;
            }

            public void InitialSort(BitArray oneSite)
            {

                #region step 1 local sum
                Parallel.For(0, Program.nThread, (i) =>
                {
                    //range by blocks
                    int s = i * blockSize;
                    if (s >= Program.nHap)
                    { return; }

                    int e = i * blockSize + blockSize;
                    if (e > Program.nHap)
                    {
                        e = Program.nHap;
                    }


                    if (oneSite[s] == true)
                    {
                        psHolder[s] = 1;
                    }
                    else
                    {
                        psHolder[s] = 0;
                    }



                    for (int k = s + 1; k < e; k++)
                    {


                        if (oneSite[k] == true)
                        {
                            psHolder[k] = psHolder[k - 1] + 1;
                        }
                        else
                        {
                            psHolder[k] = psHolder[k - 1];
                        }

                    }

                    offsetsHolder[i] = psHolder[e - 1];
                });
                #endregion

                #region step 2 seq offset handling
                for (int i = 1; i < Program.nThread; i++)
                {
                    offsetsHolder[i] = offsetsHolder[i - 1] + offsetsHolder[i];
                }

                #endregion

                #region step 3 apply offsets
                Parallel.For(1, Program.nThread, (i) =>
                {
                    int s = i * blockSize;
                    if (s >= Program.nHap)
                    { return; }
                    int e = i * blockSize + blockSize;
                    if (e > Program.nHap)
                    {
                        e = Program.nHap;
                    }

                    for (int k = s; k < e; k++)
                    {
                        psHolder[k] = psHolder[k] + offsetsHolder[i - 1];
                    }
                });
                #endregion

                #region step 4 pps -> index settle to new p arr 
                oldPDA.zeroCnt = Program.nHap - psHolder[Program.nHap - 1];
                int oneOff = oldPDA.zeroCnt - 1;

                Parallel.For(0, Program.nHap, (i) =>
                //for (int i = 0; i < Program.nHap; i++)
                {
                    if (oneSite[i] == true)
                    {
                        oldPDA.pArr[psHolder[i] + oneOff] = i;
                        oldPDA.mLens[psHolder[i] + oneOff] = 1;
                    }
                    else
                    {
                        oldPDA.pArr[i - psHolder[i]] = i;
                        oldPDA.mLens[i - psHolder[i]] = 1;
                    }
                });
                //D Arr adjust

                oldPDA.mLens[oldPDA.pArr.First()] = 0;
                oldPDA.mLens[oldPDA.pArr[oneOff + 1]] = 0;
                #endregion



            }

            public void Run()
            {
                InitialSort(Program.panel[0]);

                for (int doneIndex = 0; doneIndex < Program.nSite - 1; doneIndex++)
                {

                    if (doneIndex >= Program.LC_THD)
                    {
                        ReportLongMatches_SL(doneIndex);
                    }

                    OneSort(doneIndex + 1);
                }

                ReportLongMatches_SL_Tail(Program.nSite - 1);
            }

        }

        public class Seq
        {
            PDA newPDA;
            PDA oldPDA;
            PDA temPDA;
            BiArr biPDA;

            public class StatArr
            {
                public int cnt = 0;
                public int[] vals;
                public StatArr(int maxSize)
                {
                    vals = new int[maxSize];
                }

                public void Add(int val)
                {
                    vals[cnt] = val;
                    cnt++;
                }

                public void Set(int val, int index)
                {
                    vals[index] = val;
                }


            }

            public class BiArr
            {
                public StatArr zeros;
                public StatArr ones;

                public int nHap;
                public BiArr(int nHapTotal)
                {
                    nHap = nHapTotal;
                    zeros = new StatArr(nHap);
                    ones = new StatArr(nHap);
                }

                public void Clear()
                {
                    zeros.cnt = 0;
                    ones.cnt = 0;
                }

            }

            public Seq()
            {
                newPDA = new PDA(Program.nHap);
                oldPDA = new PDA(Program.nHap);
                biPDA = new BiArr(Program.nHap);
            }

            public void OneSort(PDA oldPDA, PDA newPDA, BitArray site)
            {

                biPDA.Clear();

                int prvLowM_One = -1;
                int prvLowM_Zero = -1;
                int hIndex = 0;




                for (int i = 0; i < Program.nHap; i++)
                {

                    if (i == oldPDA.zeroCnt)
                    {
                        prvLowM_One = -1;
                        if (biPDA.zeros.cnt > 0)
                        {
                            prvLowM_One = 0;
                        }
                        prvLowM_Zero = -1;
                        if (biPDA.ones.cnt > 0)
                        {
                            prvLowM_Zero = 0;
                        }
                    }


                    hIndex = oldPDA.pArr[i];
                    //incoming is 0
                    if (site[hIndex] == false)
                    {
                        biPDA.zeros.Add(hIndex);

                        newPDA.mLens[hIndex] = Math.Min(oldPDA.mLens[hIndex], prvLowM_One) + 1;
                        prvLowM_One = int.MaxValue;
                        prvLowM_Zero = Math.Min(prvLowM_Zero, oldPDA.mLens[hIndex]);


                    }
                    else//incoming is 1
                    {
                        biPDA.ones.Add(hIndex);

                        newPDA.mLens[hIndex] = Math.Min(oldPDA.mLens[hIndex], prvLowM_Zero) + 1;
                        prvLowM_Zero = int.MaxValue;
                        prvLowM_One = Math.Min(prvLowM_One, oldPDA.mLens[hIndex]);
                    }

                }

                for (int i = 0; i < biPDA.zeros.cnt; i++)
                {
                    newPDA.pArr[i] = biPDA.zeros.vals[i];
                }

                int oneIndex = 0;
                for (int i = biPDA.zeros.cnt; i < Program.nHap; i++)
                {
                    newPDA.pArr[i] = biPDA.ones.vals[oneIndex];
                    oneIndex++;
                }
                newPDA.zeroCnt = biPDA.zeros.cnt;


            }

            public void InitialSort()
            {
                for (int i = 0; i < Program.nHap; i++)
                {
                    if (Program.panel[0][i] == true)
                    {
                        biPDA.ones.Add(i);
                    }
                    else
                    {
                        biPDA.zeros.Add(i);
                    }
                    oldPDA.mLens[i] = 1;
                }

                for (int i = 0; i < biPDA.zeros.cnt; i++)
                {
                    oldPDA.pArr[i] = biPDA.zeros.vals[i];
                }
                int oneIndex = 0;
                for (int i = biPDA.zeros.cnt; i < Program.nHap; i++)
                {
                    oldPDA.pArr[i] = biPDA.ones.vals[oneIndex];
                    oneIndex++;
                }

                oldPDA.mLens[oldPDA.pArr[0]] = 0;
                oldPDA.mLens[biPDA.ones.vals[0]] = 0;
            }

            public void Run()
            {

                InitialSort();

                for (int i = 1; i < Program.nSite - 1; i++)
                {

                    OneSort(oldPDA, newPDA, Program.panel[i]);
                    if (i >= Program.LLM_Len)
                    {
                        ReportLongMatches(Program.LLM_Len, newPDA, i, Program.panel);
                    }

                    temPDA = oldPDA;
                    oldPDA = newPDA;
                    newPDA = temPDA;


                }

                ReportLongMatches(Program.LLM_Len, newPDA, Program.nSite, Program.panel);
            }


            /// <summary>
            /// 
            /// </summary>
            /// <param name="len"></param>
            /// <param name="pda">new pda</param>
            /// <param name="siteIndex">sorted to this index(included)</param>
            /// <param name="sites"></param>
            /// <returns></returns>
            public void ReportLongMatches(int len, PDA pda, int siteIndex, List<BitArray> sites)
            {


                BitArray site;
                if (siteIndex < sites.Count() - 1)
                {
                    site = sites[siteIndex + 1];
                }
                else
                {
                    site = null;
                }


                int zeroCnt = 0, oneCnt = 0;
                int topIndex = 0;
                int minVal = int.MaxValue;


                //need to consider special case: last site

                for (int h = 0; h < pda.pArr.Count(); h++)
                {
                    if (pda.mLens[pda.pArr[h]] < len)
                    {
                        if (zeroCnt != 0 && oneCnt != 0)
                        {
                            for (int i = topIndex; i < h - 1; i++)
                            {
                                minVal = int.MaxValue;
                                for (int k = i + 1; k < h; k++)
                                {
                                    minVal = Math.Min(minVal, pda.mLens[pda.pArr[k]]);

                                    if (site == null)
                                    {
                                        //Program.LM_Collection[pda.pArr[i] % Program.LMC_Size] += pda.pArr[k] + minVal;

                                        Program.Sums[pda.pArr[i]] += minVal;
                                        Program.Sums[pda.pArr[i]] = Program.Sums[pda.pArr[i]] % Program.nSite;
                                    }
                                    else
                                    {
                                        if (site[pda.pArr[i]] != site[pda.pArr[k]])
                                        {

                                            Program.Sums[pda.pArr[i]] += minVal;
                                            Program.Sums[pda.pArr[i]] = Program.Sums[pda.pArr[i]] % Program.nSite;
                                        }
                                    }

                                }
                            }

                        }

                        zeroCnt = 0;
                        oneCnt = 0;
                        topIndex = h;
                    }

                    if (site != null)
                    {
                        if (site[pda.pArr[h]] == true)
                        {
                            oneCnt++;
                        }
                        else
                        {
                            zeroCnt++;
                        }
                    }
                    else
                    {
                        oneCnt = 1;
                        zeroCnt = 1;
                    }



                }
                //tail case
                if (zeroCnt != 0 && oneCnt != 0)
                {
                    for (int i = topIndex; i < pda.pArr.Count(); i++)
                    {
                        minVal = int.MaxValue;
                        for (int k = i + 1; k < pda.pArr.Count(); k++)
                        {
                            minVal = Math.Min(minVal, pda.mLens[pda.pArr[k]]);

                            if (site == null)
                            {

                                Program.Sums[pda.pArr[i]] += minVal;
                                Program.Sums[pda.pArr[i]] = Program.Sums[pda.pArr[i]] % Program.nSite;
                            }
                            else
                            {
                                if (site[pda.pArr[i]] != site[pda.pArr[k]])
                                {

                                    Program.Sums[pda.pArr[i]] += minVal;
                                    Program.Sums[pda.pArr[i]] = Program.Sums[pda.pArr[i]] % Program.nSite;
                                }
                            }

                        }
                    }

                }





            }


        }
    }
}
