using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO.Compression;

namespace HP_PBWT
{
    public class utl
    {

        public static string IntTo_BMK_Str(int n)
        {
            if (n >= 1000000000)
            {
                return (n / 1000000000).ToString() + "b";
            }
            else if (n >= 1000000)
            {
                return (n / 1000000).ToString() + "m";
            }
            else if (n >= 1000)
            {
                return (n / 1000).ToString() + "k";
            }
            else
            {
                return n.ToString();
            }


        }

        public class RoundTableReaderV2
        {

            public class HapMap
            {

                int LastHID_FirstBlk = -1;
                int charOff_SecBlk = -1;
                int secBlkID = -1;
                lineTag tag;

                public HapMap(lineTag tag, int nHap)
                {
                    this.tag = tag;
                    if (tag.dataCross == false)
                    {
                        LastHID_FirstBlk = nHap - 1;
                        return;
                    }

                    LastHID_FirstBlk = (Program.rtrBlockSize - tag.charIndex - 1) / 2;
                    charOff_SecBlk = (Program.rtrBlockSize - tag.charIndex) % 2;
                    secBlkID = (tag.blockID + 1) % Program.rtrN_Block;
                }

                public char Get_HapVal(int hID)
                {
                    if (tag.dataCross == false || hID <= LastHID_FirstBlk)
                    {
                        return Program.BufferArr[tag.blockID][tag.charIndex + hID * 2];
                    }

                    return Program.BufferArr[secBlkID][charOff_SecBlk + (hID - LastHID_FirstBlk - 1) * 2];
                }

            }

            public class lineTag
            {

                public int blockID;
                public int charIndex;
                public bool dataCross;
                public bool headerCross;

                public lineTag(int blockID, int charIndex, bool dataCross, bool headerCross = false)
                {
                    this.blockID = blockID;
                    this.charIndex = charIndex;
                    this.dataCross = dataCross;
                    this.headerCross = headerCross;

                }
            }

            string VCF_Path;

            ConcurrentQueue<lineTag> Lines = new ConcurrentQueue<lineTag>();


            public int ReadLine_BlkIndex1;
            public int ReadLine_BlkIndex2;
            public int Blk1_rPTR;
            public int Blk2_rPTR;

            int L1_blockSize = Program.rtrBlockSize;
            int L1_nBlock = 20;
            int L2_nBlock = 1000;
            int[] L1_nRead;
            int[] L1_Stat;//-1 ready to add, 0 ready to parse, -2 parser done
            int[] L1_Seq;
            int L1_rPTR = 0;
            int L1_Add_Index = 0;//currently working on
            int L1_Taker_Index = 0;//currently working on
            public List<char[]> BufferArr = new List<char[]>();

            int msAddWait = 200;
            int msTakeWait = 100;
            bool L1_ReadComplete = false;
            bool L2_ReadComplete = false;


            int dataLineLen;


            public RoundTableReaderV2(string path)
            {

                VCF_Path = path;
                Console.Write(DateTime.Now + " Initializing RTR, " + L1_blockSize + " by " + L1_nBlock + " ...");

                L1_nRead = new int[L1_nBlock];
                L1_Stat = new int[L1_nBlock];
                L1_Seq = new int[L1_nBlock];
                for (int i = 0; i < L1_nBlock; i++)
                {
                    L1_Stat[i] = -1;
                    L1_nRead[i] = -1;
                    L1_Seq[i] = -1;

                    BufferArr.Add(new char[L1_blockSize]);
                }
                Console.WriteLine("Done.");

                //Task t1 = Task.Run(() => BlockReader(path));

                //ProbTop();




                //Console.WriteLine(DateTime.Now + " Initializing HP-PBWT" + Program.nHap + " Samples (#haplotypes) " + Program.nThread + " Threads...");



                //Console.WriteLine(DateTime.Now + " HP-PBWT Initialized.");



                //Task t2 = Task.Run(() => LineParser_v2());

                //Task t3 = Task.Run(() => LineTaker_v2());

                //Task.WaitAll(t1, t2,t3);

            }

            public void BlockReader()
            {
                StreamReader sr = new StreamReader(VCF_Path);
                int nBlockRead = 0;

                int len = -1;
                while (len != 0)
                {
                    while (L1_Stat[L1_Add_Index] != -1)
                    {
                        Thread.Sleep(msAddWait);
                        continue;
                    }
                    len = sr.ReadBlock(BufferArr[L1_Add_Index]);
                    L1_nRead[L1_Add_Index] = len;
                    L1_Stat[L1_Add_Index] = 0;

                    L1_Seq[L1_Add_Index] = nBlockRead;
                    nBlockRead++;

                    L1_Add_Index++;
                    L1_Add_Index = L1_Add_Index % L1_nBlock;
                }

                sr.Close();
                L1_ReadComplete = true;
                Console.WriteLine();
                Console.WriteLine(DateTime.Now + " Block Reader Complete.");
            }

            public void ProbTop()
            {

                while (L1_Stat[0] != 0)
                {
                    Thread.Sleep(msTakeWait);
                    continue;
                }


                //skip top headers
                while (BufferArr[0][L1_rPTR] == '#' && BufferArr[0][L1_rPTR + 1] == '#')
                {
                    while (BufferArr[0][L1_rPTR] != '\n')
                    {
                        L1_rPTR++;
                    }
                    L1_rPTR++;
                }
                //count nIndv


                Program.nIndv = 0;
                while (BufferArr[0][L1_rPTR] != '\n')
                {
                    if (BufferArr[0][L1_rPTR] == '\t')
                    {
                        Program.nIndv++;
                    }
                    L1_rPTR++;
                }
                Program.nIndv = Program.nIndv - 8;
                Program.nHap = Program.nIndv * 2;
                dataLineLen = Program.nIndv * 4;
                Console.WriteLine("nIndv. " + Program.nIndv + " dataLineLen " + dataLineLen);
                L1_rPTR++;


            }


            public void LineParser()
            {
                bool headerCross = false;

                int tabCnt = 0;
                int nextBlock;
                int i;

                while (L1_ReadComplete == false || L1_Stat[L1_Taker_Index] == 0)
                {//loop through blocks

                    if (L1_Stat[L1_Taker_Index] != 0 || Lines.Count > L2_nBlock)
                    {
                        Thread.Sleep(msTakeWait);
                        continue;
                    }

                    nextBlock = (L1_Taker_Index + 1) % L1_nBlock;

                    for (i = L1_rPTR; i < L1_nRead[L1_Taker_Index]; i++)
                    {//loop in a block

                        if (BufferArr[L1_Taker_Index][i] == '\t')
                        {
                            tabCnt++;

                            if (tabCnt == 9)
                            {
                                tabCnt = 0;
                                i++;
                                //cross block or end block
                                if (i + dataLineLen >= L1_nRead[L1_Taker_Index])
                                {

                                    L1_rPTR = (i + dataLineLen) % L1_nRead[L1_Taker_Index];

                                    while (L1_Stat[nextBlock] != 0)
                                    {
                                        Thread.Sleep(msTakeWait);
                                        continue;
                                    }

                                    Lines.Enqueue(new lineTag(L1_Taker_Index, i, true, false));


                                    break;

                                }
                                else
                                {//within 

                                    if (headerCross == true)
                                    {

                                        Lines.Enqueue(new lineTag(L1_Taker_Index, i, false, true));

                                        headerCross = false;
                                    }
                                    else
                                    {
                                        Lines.Enqueue(new lineTag(L1_Taker_Index, i, false, false));
                                    }


                                    //jump
                                    i += dataLineLen;


                                }
                            }
                        }

                    }

                    //move to next block
                    while (L1_Stat[nextBlock] != 0)
                    {
                        if (L1_ReadComplete)
                        {
                            L2_ReadComplete = true;


                            Console.WriteLine();
                            Console.WriteLine(DateTime.Now + " L2 Read Complete.");

                            return;
                        }
                        Thread.Sleep(msTakeWait);
                        continue;
                    }

                    if (tabCnt != 0)
                    {
                        L1_rPTR = 0;
                        //L1_Stat[L1_Taker_Index] = -2;
                        headerCross = true;
                    }
                    else
                    {
                        L1_Stat[L1_Taker_Index] = -2;
                    }

                    L1_Taker_Index = nextBlock;


                }

                L2_ReadComplete = true;

                Console.WriteLine();
                Console.WriteLine(DateTime.Now + " L2 Read Complete.");


            }

            public void LineTaker_LLM()
            {

                lineTag tag;
                bool dq;

                int doneIndex = 0;
                //first site
                while (L2_ReadComplete == false || Lines.Count > 0)
                {
                    while (Lines.Count == 0)
                    {
                        Thread.Sleep(msTakeWait);
                        continue;
                    }

                    dq = Lines.TryDequeue(out tag);
                    if (dq == false)
                    { continue; }


                    HapMap map = new HapMap(tag, Program.nHap);

                    Program.pal.InitialSort(map);



                    break;
                }
                //Console.WriteLine(doneIndex + " D 19:" + Program.pal.oldPDA.mLens[19]);

                //toWorkIndex++;

                //mid sites

                while (L2_ReadComplete == false || Lines.Count > 0)
                {
                    while (Lines.Count == 0)
                    {
                        Thread.Sleep(msTakeWait);
                        continue;
                    }

                    dq = Lines.TryDequeue(out tag);
                    if (dq == false)
                    { continue; }


                    LineWork_LLM(doneIndex, tag);

                    if (tag.dataCross == true)
                    {
                        L1_Stat[tag.blockID] = -1;
                    }

                    if (tag.headerCross == true)
                    {
                        L1_Stat[(tag.blockID - 1 + L1_nBlock) % L1_nBlock] = -1;
                    }

                    doneIndex++;
                    //Console.WriteLine("site "+doneIndex);

                }


                Console.WriteLine(DateTime.Now + " L3 Line Taker Complete.");

                //PBWT last site process
                //last site already done PBWT need to output LM or SSM

                Program.pal.ReportLongMatches_SL_Tail(doneIndex);

                Program.BW.DoneAdding();
            }

            void LineWork_LLM(int doneIndex, lineTag tag)
            {
                //Console.WriteLine("Line: " + lineCnt);

                HapMap map = new HapMap(tag, Program.nHap);
                //Console.WriteLine("LW " + doneIndex);
                if (doneIndex >= Program.LC_THD)
                {
                    Program.pal.ReportLongMatches_SL(doneIndex, map);
                }

                Program.pal.OneSort(map, doneIndex);
            }





            public void LineTaker_SMM()
            {

                lineTag tag;
                bool dq;

                int siteCnt = 1;
                //first site
                while (L2_ReadComplete == false || Lines.Count > 0)
                {
                    while (Lines.Count == 0)
                    {
                        Thread.Sleep(msTakeWait);
                        continue;
                    }

                    dq = Lines.TryDequeue(out tag);
                    if (dq == false)
                    { continue; }


                    HapMap map = new HapMap(tag, Program.nHap);

                    Program.pal.InitialSort(map);

                    break;
                }

                siteCnt++;

                //mid sites
                while (L2_ReadComplete == false || Lines.Count > 0)
                {
                    while (Lines.Count == 0)
                    {
                        Thread.Sleep(msTakeWait);
                        continue;
                    }

                    dq = Lines.TryDequeue(out tag);
                    if (dq == false)
                    { continue; }


                    LineWork_SMM(siteCnt, tag);

                    if (tag.dataCross == true)
                    {
                        L1_Stat[tag.blockID] = -1;
                    }

                    if (tag.headerCross == true)
                    {
                        L1_Stat[(tag.blockID - 1 + L1_nBlock) % L1_nBlock] = -1;
                    }

                    siteCnt++;
                }


                Console.WriteLine();
                Console.WriteLine(DateTime.Now + " L3 Line Taker Complete.");

                //PBWT last site process
                //last site already done PBWT need to output LM or SSM
                Program.pal.ReportSetMaxMatch_Tail(siteCnt);

                Program.BW.DoneAdding();
            }

            void LineWork_SMM(int lineCnt, lineTag tag)
            {
                HapMap map = new HapMap(tag, Program.nHap);

                if (lineCnt >= Program.LLM_Len)
                {
                    Program.pal.ReportSetMaxMatch(lineCnt, map);
                }

                Program.pal.OneSort(map, 0);

            }

            public void LineTaker_DBG()
            {

                string line;

                while ((line = Program.dbgSR.ReadLine()) != null && line.StartsWith("##"))
                {
                    continue;
                }


                lineTag tag;
                bool dq, pk;

                int cnt = 0;

                while (L2_ReadComplete == false || Lines.Count > 0)
                {
                    while (Lines.Count == 0)
                    {
                        Thread.Sleep(msTakeWait);
                        continue;
                    }

                    dq = Lines.TryDequeue(out tag);
                    if (dq == false)
                    { continue; }

                    Console.WriteLine("Line " + cnt);

                    if (cnt == 9830)
                    {
                        int a = 0;
                    }

                    LineWork_DBG(tag);


                    cnt++;

                    if (tag.dataCross == true)
                    {
                        L1_Stat[tag.blockID] = -1;
                    }

                    if (tag.headerCross == true)
                    {
                        L1_Stat[(tag.blockID - 1 + L1_nBlock) % L1_nBlock] = -1;
                    }
                }


                Console.WriteLine();
                Console.WriteLine(DateTime.Now + " L3 Line Taker Complete.");

                Program.dbgSR.Close();
            }

            void LineWork_DBG(lineTag tag)
            {
                RoundTableReaderV2.HapMap map = new HapMap(tag, Program.nHap);

                string line;
                string[] parts;
                line = Program.dbgSR.ReadLine();
                parts = line.Split('\t');

                List<char> hapChars = new List<char>();
                for (int i = 9; i < parts.Count(); i++)
                {
                    hapChars.Add(parts[i][0]);
                    hapChars.Add(parts[i][2]);
                }
                Parallel.For(0, hapChars.Count(), (i) =>
                //for (int i = 0; i < hapChars.Count(); i++)
                {

                    if (hapChars[i] != map.Get_HapVal(i))
                    {
                        Console.WriteLine(i + ": " + hapChars[i] + " vs " + map.Get_HapVal(i));
                        char c = map.Get_HapVal(i);
                        int a = 0;
                        c = map.Get_HapVal(i);
                    }
                });
            }





        }

        public class BufferWriter
        {
            public ConcurrentQueue<string> buffer = new ConcurrentQueue<string>();
            StreamWriter sw;
            string outPath;
            bool AddComplete = false;
            int msWait = 500;

            public void Add(string s)
            {
                buffer.Enqueue(s);
            }

            public void DoneAdding()
            {
                AddComplete = true;
            }

            public BufferWriter(string outPath)
            {
                this.outPath = outPath;
            }

            public void Run()
            {
                sw = new StreamWriter(outPath);
                sw.NewLine = "\n";
                string sOut;
                while (AddComplete == false || buffer.Count() > 0)
                {
                    if (buffer.TryDequeue(out sOut) == false)
                    {
                        Thread.Sleep(msWait);
                    }
                    if (String.IsNullOrEmpty(sOut))
                    { continue; }
                    sw.WriteLine(sOut);
                }

                sw.Close();

                Console.WriteLine(DateTime.Now + " Buffer Writer Complete.");
                Environment.Exit(0);
            }


        }

        
    }
}
