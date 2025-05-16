// utl.h
#ifndef UTL_H
#define UTL_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>

namespace HP_PBWT::utl {

    std::string IntTo_BMK_Str(int n);

    class RoundTableReaderV2 {
    public:
        class lineTag {
        public:
            int blockID;
            int charIndex;
            bool dataCross;
            bool headerCross;
            lineTag(int blockID, int charIndex, bool dataCross, bool headerCross = false);
        };

        class HapMap {
        public:
            HapMap(lineTag tag, int nHap);
            char Get_HapVal(int hID);
        private:
            int LastHID_FirstBlk;
            int charOff_SecBlk;
            int secBlkID;
            lineTag tag;
        };

        RoundTableReaderV2(const std::string& path);
        void BlockReader();
        void ProbTop();
        void LineParser();
        void LineTaker_LLM();

    private:
        std::string VCF_Path;
        std::queue<lineTag> Lines;
        int ReadLine_BlkIndex1;
        int ReadLine_BlkIndex2;
        int Blk1_rPTR;
        int Blk2_rPTR;
        int L1_blockSize;
        int L1_nBlock;
        int L2_nBlock;
        std::vector<int> L1_nRead;
        std::vector<int> L1_Stat;
        std::vector<int> L1_Seq;
        int L1_rPTR;
        int L1_Add_Index;
        int L1_Taker_Index;
        std::vector<std::vector<char>> BufferArr;
        int msAddWait;
        int msTakeWait;
        bool L1_ReadComplete;
        bool L2_ReadComplete;
        int dataLineLen;
    };

    class BufferWriter {
        public:
            BufferWriter(const std::string& outPath);
            void Add(const std::string& s);
            void DoneAdding();
            void Run();
        
        private:
            std::queue<std::string> buffer;
            std::ofstream sw;
            std::string outPath;
            bool addComplete;
            std::mutex mtx;
            std::condition_variable cv;
        };

} // namespace HP_PBWT::utl

#endif