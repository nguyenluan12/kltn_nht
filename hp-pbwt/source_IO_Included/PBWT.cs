using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HP_PBWT
{
    public class PBWT
    {



        public class Pal
        {
            public class PDA
            {
                public int[] pArr;

                public int zeroCnt = 0;

                public int[] mLens;

                //public int nHap;
                public PDA(int nHapTotal)
                {

                    pArr = new int[nHapTotal];

                    mLens = new int[nHapTotal];

                }

            }

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

            static PDA newPDA;
            static PDA oldPDA;
            static PDA temPDA;
            static int[] psHolder;
            static int[] minDZ;
            static bool[] minDZ_Ready;

            static int[] offsetsHolder;

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



            /// <summary>
            /// us psHolder create new PDA
            /// </summary>
            /// <param name="map"></param>
            public void coreP_Arr(utl.RoundTableReaderV2.HapMap map)
            {
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

                    if (map.Get_HapVal(oldPDA.pArr[s]) == '0')
                    {
                        psHolder[s] = 0;
                    }
                    else
                    {
                        psHolder[s] = 1;
                    }

                    for (int k = s + 1; k < e; k++)
                    {
                        //psHolder[k] = psHolder[k - 1] + site[oldPDA.pArr[k]];
                        if (map.Get_HapVal(oldPDA.pArr[k]) == '0')
                        {
                            psHolder[k] = psHolder[k - 1];

                        }
                        else
                        {
                            psHolder[k] = psHolder[k - 1] + 1;
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
                    if (map.Get_HapVal(oldPDA.pArr[i]) == '0')
                    {
                        newPDA.pArr[i - psHolder[i]] = oldPDA.pArr[i];

                    }
                    else
                    {
                        newPDA.pArr[psHolder[i] + oneOff] = oldPDA.pArr[i];
                    }
                });
                #endregion

            }


            public void coreD_Arr(utl.RoundTableReaderV2.HapMap map)
            {


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
                    int seekHID;

                    //minDZ[i] = int.MaxValue;


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
                            if (map.Get_HapVal(hID) == '0')
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

                        //ToDo: minZero and minOne should only happen once
                        if (map.Get_HapVal(oldPDA.pArr[s - 1]) == '0')
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
                            while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) != '0')
                            {
                                seekIndex--;
                            }
                            if (seekIndex >= 0)
                            {//ToDo: may need to change
                                prvLowM_Zero = oldPDA.mLens[oldPDA.pArr[seekIndex]];
                                while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) == '0')
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
                            while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) == '0')
                            {
                                seekIndex--;
                            }
                            if (seekIndex >= 0)
                            {
                                prvLowM_One = oldPDA.mLens[oldPDA.pArr[seekIndex]];
                                while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) != '0')
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
                            if (map.Get_HapVal(hID) == '0')
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


            public void coreD_Arr_LowMAF(utl.RoundTableReaderV2.HapMap map, int doneIndex)
            {
                ResetMinDZ();

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

                    //Console.WriteLine("Blk"+i + " [" + s + "," + e + ") started.");
                    int prvLowM_Zero = -1;
                    int prvLowM_One = -1;
                    int hID;
                    int seekHID;



                    //minDZ[i] = int.MaxValue;

                    int temMinDZ = int.MaxValue;

                    #region first block 
                    //no up search
                    if (s == 0)
                    {
                        if (psHolder[e - 1] != 0)
                        {//block has one

                            minDZ[i] = blockHasOneSignal;
                            minDZ_Ready[i] = true;
                            //Console.WriteLine("Blk" + i + " [" + s + "," + e + ") MinDZ="+minDZ[i]);
                            for (int k = s; k < e; k++)
                            {

                                hID = oldPDA.pArr[k];

                                //incoming is 0
                                if (map.Get_HapVal(hID) == '0')
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
                        else// no one in the block
                        {
                            newPDA.mLens[oldPDA.pArr[0]] = 0;
                            for (int k = s + 1; k < e; k++)
                            {
                                hID = oldPDA.pArr[k];
                                newPDA.mLens[hID] = oldPDA.mLens[hID] + 1;
                            }

                        }
                        //Console.WriteLine("Blk"+i+"->"+minDZ[i]+" [" + s + "," + e + ") ended.");
                        return;
                    }
                    #endregion


                    bool minZeroSearch = true;
                    bool minOneSearch = true;

                    #region not first block, this block has only 0s
                    if (psHolder[s - 1] == psHolder[e - 1])
                    {
                        #region first row in this block
                        //upper search for a '0'
                        hID = oldPDA.pArr[s];
                        if (psHolder[s - 1] == s)//the first zero
                        {
                            newPDA.mLens[hID] = 0;
                        }
                        else
                        {
                            prvLowM_Zero = oldPDA.mLens[oldPDA.pArr[s]];
                            int seekIndex = s - 1;

                            while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) != '0')
                            {
                                prvLowM_Zero = Math.Min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_Zero);
                                seekIndex--;
                            }

                        }
                        newPDA.mLens[hID] = prvLowM_Zero + 1;

                        temMinDZ = newPDA.mLens[hID];
                        #endregion
                        //all other rows in this block
                        for (int k = s + 1; k < e; k++)
                        {
                            hID = oldPDA.pArr[k];
                            newPDA.mLens[hID] = oldPDA.mLens[hID] + 1;
                            temMinDZ = Math.Min(temMinDZ, newPDA.mLens[hID]);
                        }

                        minDZ[i] = temMinDZ;
                        minDZ_Ready[i] = true;

                        //Console.WriteLine("Blk" + i + "->" + minDZ[i] + " [" + s + "," + e + ") ended.");
                        return;
                    }
                    #endregion

                    #region not first block, this block has 1s

                    minDZ[i] = blockHasOneSignal;
                    minDZ_Ready[i] = true;
                    //special cases:
                    //A there is no 0 in upper blocks
                    //B there is no 1 in upper blocks
                    if (psHolder[s - 1] == s) //A there is no 0 in upper blocks
                    {
                        prvLowM_Zero = -1;
                        minZeroSearch = false;
                    }

                    if (psHolder[s - 1] == 0) //B there is no 1 in upper blocks
                    {
                        prvLowM_One = -1;
                        minOneSearch = false;
                    }

                    //upper search

                    if (minZeroSearch)
                    {
                        //locate first upper zero

                        prvLowM_Zero = oldPDA.mLens[oldPDA.pArr[s]];
                        int seekIndex = s - 1;
                        while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) != '0')
                        {
                            prvLowM_Zero = Math.Min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_Zero);
                            seekIndex--;
                        }

                    }

                    if (minOneSearch)
                    {
                        int seekIndex = s - 1;
                        int search_BlockID = i - 1;

                        prvLowM_One = oldPDA.mLens[oldPDA.pArr[seekIndex]];

                        //jump blocks
                        //while (minDZ[search_BlockID] != blockHasOneSignal)
                        //{
                        //    if (minDZ[search_BlockID]==int.MaxValue)
                        //    {//not ready
                        //        Console.WriteLine(doneIndex +" Blk" + i + " [" + s + "," + e + ") checking " +search_BlockID+" ->"+ minDZ[search_BlockID] +  " wait.");
                        //        Thread.Sleep(msD_Wait);
                        //        continue;
                        //    }

                        //    prvLowM_One = Math.Min(minDZ[search_BlockID], prvLowM_One);

                        //    search_BlockID--;
                        //}


                        while (true)
                        {
                            if (minDZ_Ready[search_BlockID] == false)
                            {//not ready
                                //Console.WriteLine(doneIndex + " Blk" + i + " [" + s + "," + e + ") checking " + search_BlockID + " ->" + minDZ[search_BlockID] + " wait.");
                                //Thread.Sleep(msD_Wait);
                                continue;
                            }
                            if (minDZ[search_BlockID] != blockHasOneSignal)
                            {
                                prvLowM_One = Math.Min(minDZ[search_BlockID] - 1, prvLowM_One);
                                search_BlockID--;
                            }
                            else
                            {
                                break;
                            }

                        }

                        seekIndex = search_BlockID * blockSize + blockSize - 1;
                        //reached the block with 1, locate first upper one.
                        while (map.Get_HapVal(oldPDA.pArr[seekIndex]) == '0')
                        {
                            prvLowM_One = Math.Min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_One);
                            seekIndex--;
                        }



                    }


                    hID = oldPDA.pArr[s];

                    //first row 
                    if (map.Get_HapVal(oldPDA.pArr[s - 1]) == map.Get_HapVal(oldPDA.pArr[s]))
                    {
                        //match continue
                        //incoming is 0
                        if (map.Get_HapVal(hID) == '0')
                        {
                            newPDA.mLens[hID] = oldPDA.mLens[hID] + 1;
                            prvLowM_One = int.MaxValue;
                            prvLowM_Zero = Math.Min(prvLowM_Zero, oldPDA.mLens[hID]);

                            //temMinDZ = Math.Min(newPDA.mLens[hID], temMinDZ);
                        }
                        else//incoming is 1
                        {
                            newPDA.mLens[hID] = oldPDA.mLens[hID] + 1;
                            prvLowM_Zero = int.MaxValue;
                            prvLowM_One = Math.Min(prvLowM_One, oldPDA.mLens[hID]);

                            //temMinDZ = blockHasOneSignal;
                        }

                    }
                    else
                    {
                        //No match
                        //incoming is 0
                        if (map.Get_HapVal(hID) == '0')
                        {
                            newPDA.mLens[hID] = Math.Min(prvLowM_One, oldPDA.mLens[hID]) + 1;
                            prvLowM_One = int.MaxValue;
                            prvLowM_Zero = Math.Min(prvLowM_Zero, oldPDA.mLens[hID]);

                            //temMinDZ = Math.Min(newPDA.mLens[hID], minDZ[i]);
                        }
                        else//incoming is 1
                        {
                            newPDA.mLens[hID] = Math.Min(prvLowM_Zero, oldPDA.mLens[hID]) + 1;
                            prvLowM_Zero = int.MaxValue;
                            prvLowM_One = Math.Min(prvLowM_One, oldPDA.mLens[hID]);

                            //temMinDZ = blockHasOneSignal;
                        }
                    }


                    for (int k = s + 1; k < e; k++)
                    {

                        hID = oldPDA.pArr[k];

                        //incoming is 0
                        if (map.Get_HapVal(hID) == '0')
                        {
                            newPDA.mLens[hID] = Math.Min(oldPDA.mLens[hID], prvLowM_One) + 1;
                            prvLowM_One = int.MaxValue;
                            prvLowM_Zero = Math.Min(prvLowM_Zero, oldPDA.mLens[hID]);

                            //temMinDZ = Math.Min(newPDA.mLens[hID], temMinDZ);
                        }
                        else//incoming is 1
                        {
                            newPDA.mLens[hID] = Math.Min(oldPDA.mLens[hID], prvLowM_Zero) + 1;
                            prvLowM_Zero = int.MaxValue;
                            prvLowM_One = Math.Min(prvLowM_One, oldPDA.mLens[hID]);

                            //temMinDZ = blockHasOneSignal;
                        }
                    }

                    //Console.WriteLine("Blk" + i + "->" + minDZ[i] + " [" + s + "," + e + ") ended.");
                    //minDZ[i] = temMinDZ;

                    #endregion


                    //}
                });


            }

            public void coreD_Arr_BU_a(utl.RoundTableReaderV2.HapMap map)
            {
                ResetMinDZ();

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

                    //Console.WriteLine("Blk"+i + " [" + s + "," + e + ") started.");
                    int prvLowM_Zero = -1;
                    int prvLowM_One = -1;
                    int hID;
                    int seekHID;



                    //minDZ[i] = int.MaxValue;

                    int temMinDZ = int.MaxValue;

                    #region first block 
                    //no up search
                    if (s == 0)
                    {
                        if (psHolder[e - 1] != 0)
                        {//block has one

                            minDZ[i] = blockHasOneSignal;
                            minDZ_Ready[i] = true;
                            //Console.WriteLine("Blk" + i + " [" + s + "," + e + ") MinDZ="+minDZ[i]);
                            for (int k = s; k < e; k++)
                            {

                                hID = oldPDA.pArr[k];

                                //incoming is 0
                                if (map.Get_HapVal(hID) == '0')
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
                        else// no one in the block
                        {
                            newPDA.mLens[oldPDA.pArr[0]] = 0;
                            for (int k = s + 1; k < e; k++)
                            {
                                hID = oldPDA.pArr[k];
                                newPDA.mLens[hID] = oldPDA.mLens[hID] + 1;
                            }

                        }
                        //Console.WriteLine("Blk"+i+"->"+minDZ[i]+" [" + s + "," + e + ") ended.");
                        return;
                    }
                    #endregion


                    bool minZeroSearch = true;
                    bool minOneSearch = true;

                    #region not first block, this block has only 0s
                    if (psHolder[s - 1] == psHolder[e - 1])
                    {
                        #region first row in this block
                        //upper search for a '0'
                        hID = oldPDA.pArr[s];
                        if (psHolder[s - 1] == s)//the first zero
                        {
                            newPDA.mLens[hID] = 0;
                        }
                        else
                        {
                            prvLowM_Zero = oldPDA.mLens[oldPDA.pArr[s]];
                            int seekIndex = s - 1;

                            while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) != '0')
                            {
                                prvLowM_Zero = Math.Min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_Zero);
                                seekIndex--;
                            }

                        }
                        newPDA.mLens[hID] = prvLowM_Zero + 1;

                        temMinDZ = newPDA.mLens[hID];
                        #endregion
                        //all other rows in this block
                        for (int k = s + 1; k < e; k++)
                        {
                            hID = oldPDA.pArr[k];
                            newPDA.mLens[hID] = oldPDA.mLens[hID] + 1;
                            temMinDZ = Math.Min(temMinDZ, newPDA.mLens[hID]);
                        }

                        minDZ[i] = temMinDZ;
                        minDZ_Ready[i] = true;

                        //Console.WriteLine("Blk" + i + "->" + minDZ[i] + " [" + s + "," + e + ") ended.");
                        return;
                    }
                    #endregion

                    #region not first block, this block has 1s

                    minDZ[i] = blockHasOneSignal;
                    minDZ_Ready[i] = true;

                    #region go down block
                    int firstprvLowM_One = int.MaxValue;
                    int firstprvLowM_Zero = int.MaxValue;

                    int firstZeroIndex = -1, firstOneIndex = -1;
                    bool firstOneFound = false;
                    bool firstZeroFound = false;

                    for (int k = s; k < e; k++)
                    {

                        hID = oldPDA.pArr[k];

                        //incoming is 0
                        if (map.Get_HapVal(hID) == '0')
                        {
                            if (!firstZeroFound)
                            {
                                firstZeroFound = true;
                                firstZeroIndex = k;
                                //continue;
                            }
                            if (!firstOneFound)
                            {
                                firstprvLowM_One = Math.Min(oldPDA.mLens[hID], firstprvLowM_One);
                            }


                            newPDA.mLens[hID] = Math.Min(oldPDA.mLens[hID], prvLowM_One) + 1;
                            prvLowM_One = int.MaxValue;
                            prvLowM_Zero = Math.Min(prvLowM_Zero, oldPDA.mLens[hID]);

                            //temMinDZ = Math.Min(newPDA.mLens[hID], temMinDZ);
                        }
                        else//incoming is 1
                        {
                            if (!firstOneFound)
                            {
                                firstOneFound = true;
                                firstOneIndex = k;
                                //continue;
                            }
                            if (!firstZeroFound)
                            {
                                firstprvLowM_Zero = Math.Min(oldPDA.mLens[hID], firstprvLowM_Zero);
                            }
                            newPDA.mLens[hID] = Math.Min(oldPDA.mLens[hID], prvLowM_Zero) + 1;
                            prvLowM_Zero = int.MaxValue;
                            prvLowM_One = Math.Min(prvLowM_One, oldPDA.mLens[hID]);

                            //temMinDZ = blockHasOneSignal;
                        }
                    }
                    #endregion


                    #region compute first 0 and first 1 divergence values
                    //special cases:
                    //A there is no 0 in upper blocks
                    //B there is no 1 in upper blocks
                    if (psHolder[s - 1] == s) //A there is no 0 in upper blocks
                    {
                        prvLowM_Zero = -1;
                        minZeroSearch = false;
                    }

                    if (psHolder[s - 1] == 0) //B there is no 1 in upper blocks
                    {
                        prvLowM_One = -1;
                        minOneSearch = false;
                    }

                    //upper search

                    if (minZeroSearch)
                    {
                        //locate first upper zero

                        prvLowM_Zero = oldPDA.mLens[oldPDA.pArr[s]];
                        int seekIndex = s - 1;
                        while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) != '0')
                        {
                            prvLowM_Zero = Math.Min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_Zero);
                            seekIndex--;
                        }

                    }




                    hID = oldPDA.pArr[s];

                    //compute firstprvLowM_one and firstprvLowM_Zero


                    if (minOneSearch)
                    {
                        int seekIndex = s - 1;
                        int search_BlockID = i - 1;

                        prvLowM_One = oldPDA.mLens[oldPDA.pArr[seekIndex]];


                        while (true)
                        {
                            if (minDZ_Ready[search_BlockID] == false)
                            {//not ready
                                //Console.WriteLine(doneIndex + " Blk" + i + " [" + s + "," + e + ") checking " + search_BlockID + " ->" + minDZ[search_BlockID] + " wait.");
                                Thread.Sleep(msD_Wait);
                                continue;
                            }


                            if (minDZ[search_BlockID] != blockHasOneSignal)
                            {
                                prvLowM_One = Math.Min(minDZ[search_BlockID] - 1, prvLowM_One);
                                search_BlockID--;
                            }
                            else
                            {
                                break;
                            }

                        }

                        seekIndex = search_BlockID * blockSize + blockSize - 1;
                        //reached the block with 1, locate first upper one.
                        while (map.Get_HapVal(oldPDA.pArr[seekIndex]) == '0')
                        {
                            prvLowM_One = Math.Min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_One);
                            seekIndex--;
                        }
                    }


                    if (firstOneIndex != -1)
                    {
                        prvLowM_One = Math.Min(prvLowM_One, firstprvLowM_One);
                        newPDA.mLens[oldPDA.pArr[firstOneIndex]] = Math.Min(oldPDA.mLens[oldPDA.pArr[firstOneIndex]], prvLowM_One) + 1;
                    }
                    else
                    {
                        Console.WriteLine("Exception D1.");
                    }

                    if (firstZeroIndex != -1)
                    {
                        prvLowM_Zero = Math.Min(prvLowM_Zero, firstprvLowM_Zero);
                        newPDA.mLens[oldPDA.pArr[firstZeroIndex]] = Math.Min(oldPDA.mLens[oldPDA.pArr[firstZeroIndex]], prvLowM_Zero) + 1;
                    }

                    #endregion
                    #endregion


                    //}
                });


            }

            //
            /// <summary>
            /// 
            /// </summary>
            /// <param name="siteIndex">sorted, 0 based</param>
            /// <param name="map">site + 1's data  </param>
            public void ReportLongMatches_SL(int siteIndex, utl.RoundTableReaderV2.HapMap map)
            {
                Parallel.For(0, Program.nHap - 1, (h) =>
                //for (int h = 0; h < Program.nHap-1; h++)
                {
                    int minVal = int.MaxValue;
                    for (int i = h + 1; i < Program.nHap; i++)
                    {

                        if (oldPDA.mLens[oldPDA.pArr[i]] < Program.LLM_Len)
                        { break; }

                        minVal = Math.Min(minVal, oldPDA.mLens[oldPDA.pArr[i]]);


                        if (map.Get_HapVal(oldPDA.pArr[h]) != map.Get_HapVal(oldPDA.pArr[i]))
                        {
                            //report
                            Program.BW.Add(oldPDA.pArr[h] + "\t" + oldPDA.pArr[i] + "\t" + siteIndex + "\t" + minVal);
                        }
                    }
                }
                );

            }

            public void ReportLongMatches_SL_Tail(int siteIndex)
            {
                Parallel.For(0, Program.nHap - 1, (h) =>
                //for (int h = 0; h < Program.nHap - 1; h++)
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
                        Program.BW.Add(oldPDA.pArr[h] + "\t" + oldPDA.pArr[i] + "\t" + siteIndex + "\t" + minVal);


                    }
                });

            }


            public void ReportSetMaxMatch(int siteIndex, utl.RoundTableReaderV2.HapMap map)
            {
                Parallel.For(0, Program.nHap, pOp, (oneChosenIndex) =>
                {
                    int rangeUp = oneChosenIndex, rangeDown = oneChosenIndex;
                    int currHID = oldPDA.pArr[oneChosenIndex];
                    int len;
                    char val = map.Get_HapVal(currHID);
                    int scanIndex;

                    bool runThough = false;


                    if (oneChosenIndex == 0)
                    {
                        len = oldPDA.mLens[oldPDA.pArr[oneChosenIndex + 1]];
                    }
                    else if (oneChosenIndex == Program.nHap - 1)
                    {
                        len = oldPDA.mLens[oldPDA.pArr[oneChosenIndex]];
                    }
                    else
                    {
                        len = Math.Max(oldPDA.mLens[oldPDA.pArr[oneChosenIndex]], oldPDA.mLens[oldPDA.pArr[oneChosenIndex + 1]]);
                    }

                    //check down
                    if (oneChosenIndex != Program.nHap - 1)
                    {
                        scanIndex = oneChosenIndex + 1;


                        while (scanIndex < Program.nHap && len <= oldPDA.mLens[oldPDA.pArr[scanIndex]])
                        {
                            if (val == map.Get_HapVal(oldPDA.pArr[scanIndex]))
                            {
                                runThough = true;
                                break;
                            }
                            scanIndex++;
                        }
                        rangeDown = scanIndex - 1;
                    }

                    if (runThough == true)
                    { return; }


                    if (oneChosenIndex != 0)
                    {//check up 
                        scanIndex = oneChosenIndex - 1;

                        //scanIndex--;
                        while (scanIndex >= 0 && len <= oldPDA.mLens[oldPDA.pArr[scanIndex + 1]])
                        {
                            if (val == map.Get_HapVal(oldPDA.pArr[scanIndex]))
                            {
                                runThough = true;
                                break;
                            }
                            scanIndex--;
                        }

                        rangeUp = scanIndex + 1;
                    }

                    if (runThough == true)
                    { return; }

                    //report 

                    for (int i = rangeUp; i <= rangeDown; i++)
                    {
                        if (i == oneChosenIndex)
                        { continue; }
                        Program.BW.Add(currHID + "\t" + oldPDA.pArr[i] + "\t" + siteIndex + " " + len);
   

                    }
                });


            }

            public void ReportSetMaxMatch_Tail(int siteIndex)
            {


                Parallel.For(0, Program.nHap, pOp, (oneChosenIndex) =>
                {
                    int rangeUp = oneChosenIndex, rangeDown = oneChosenIndex;
                    int currHID = oldPDA.pArr[oneChosenIndex];
                    int len;
                    int scanIndex;


                    if (oneChosenIndex == 0)
                    {
                        len = oldPDA.mLens[oldPDA.pArr[oneChosenIndex + 1]];
                    }
                    else if (oneChosenIndex == Program.nHap - 1)
                    {
                        len = oldPDA.mLens[oldPDA.pArr[oneChosenIndex]];
                    }
                    else
                    {
                        len = Math.Max(oldPDA.mLens[oldPDA.pArr[oneChosenIndex]], oldPDA.mLens[oldPDA.pArr[oneChosenIndex + 1]]);
                    }

                    //check down
                    if (oneChosenIndex != Program.nHap - 1)
                    {
                        scanIndex = oneChosenIndex + 1;

                        //scanIndex++;
                        while (scanIndex < Program.nHap && len <= oldPDA.mLens[oldPDA.pArr[scanIndex]])
                        {
                            scanIndex++;
                        }
                        rangeDown = scanIndex - 1;
                    }

                    if (oneChosenIndex != 0)
                    {//check up 
                        scanIndex = oneChosenIndex - 1;

                        //scanIndex--;
                        while (scanIndex >= 0 && len <= oldPDA.mLens[oldPDA.pArr[scanIndex + 1]])
                        {
                            scanIndex--;
                        }

                        rangeUp = scanIndex + 1;
                    }

                    //report 

                    for (int i = rangeUp; i <= rangeDown; i++)
                    {
                        if (i == oneChosenIndex)
                        { continue; }
                        Program.BW.Add(currHID + "\t" + oldPDA.pArr[i] + "\t" + siteIndex + " " + len);
                    }
                });





            }

            //create new P and D
            public void OneSort(utl.RoundTableReaderV2.HapMap map, int doneIndex)
            {
                coreP_Arr(map);



                coreD_Arr(map);

                temPDA = oldPDA;
                oldPDA = newPDA;
                newPDA = temPDA;

            }



            /// <summary>
            /// making the first P and D
            /// </summary>
            /// <param name="map"></param>
            public void InitialSort(utl.RoundTableReaderV2.HapMap map)
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


                    if (map.Get_HapVal(s) == '0')
                    {
                        psHolder[s] = 0;
                    }
                    else
                    {
                        psHolder[s] = 1;
                    }



                    for (int k = s + 1; k < e; k++)
                    {


                        if (map.Get_HapVal(k) == '0')
                        {
                            psHolder[k] = psHolder[k - 1];

                        }
                        else
                        {
                            psHolder[k] = psHolder[k - 1] + 1;
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
                    if (map.Get_HapVal(i) == '0')
                    {
                        oldPDA.pArr[i - psHolder[i]] = i;
                        oldPDA.mLens[i - psHolder[i]] = 1;
                    }
                    else
                    {
                        oldPDA.pArr[psHolder[i] + oneOff] = i;
                        oldPDA.mLens[psHolder[i] + oneOff] = 1;

                    }
                });
                //D Arr adjust

                oldPDA.mLens[oldPDA.pArr.First()] = 0;
                if (oldPDA.zeroCnt != Program.nHap)
                {
                    oldPDA.mLens[oldPDA.pArr[oneOff + 1]] = 0;
                }
                #endregion



            }



        }

    }
}
