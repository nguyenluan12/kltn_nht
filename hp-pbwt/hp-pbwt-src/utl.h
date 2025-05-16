#ifndef UTL_H
#define UTL_H

#include <iostream>
#include <vector>
#include <string>
#include <bitset>
#include <fstream>
#include <sstream>
#include <queue>
#include <thread>
#include <chrono>
#include <algorithm>
#include <omp.h>

namespace HP_PBWT {

class utl {
public:
    static std::string IntTo_BMK_Str(int n);

    class RoundTableReaderV2 {
    public:
        class lineTag {
        public:
            int blockID;
            int charIndex;
            bool dataCross;
            bool headerCross;

            lineTag();
            lineTag(int blockID, int charIndex, bool dataCross, bool headerCross = false);
        };

        class HapMap {
        private:
            int LastHID_FirstBlk;
            int charOff_SecBlk;
            int secBlkID;
            lineTag tag;

        public:
            HapMap(const lineTag& tag, int nHap);
            char Get_HapVal(int hID);
        };

    private:
        std::string VCF_Path;
        std::vector<int> L1_Stat;
        int msTakeWait;
        std::vector<std::vector<char>> BufferArr;
        int L1_rPTR;
        int dataLineLen;
        bool L1_ReadComplete;
        int L1_Taker_Index;
        std::queue<lineTag> Lines;
        int L2_nBlock;
        int L1_nBlock;
        bool L2_ReadComplete;
        std::vector<int> L1_nRead;

    public:
        RoundTableReaderV2(const std::string& path);
        void BlockReader();
        void ProbTop();
        void LineParser();
        void LineTaker_LLM();
        void LineTaker_SMM();
        void LineTaker_DBG();
        void LineWork_LLM(int doneIndex, lineTag tag);
        void LineWork_SMM(int lineCnt, lineTag tag);
        void LineWork_DBG(lineTag tag);
    };

    class BufferWriter {
    private:
        std::queue<std::string> buffer;
        std::ofstream sw;
        std::string outPath;
        bool AddComplete;
        int msWait;

    public:
        BufferWriter(const std::string& outPath);
        void Add(const std::string& s);
        void DoneAdding();
        void Run();
    };

    static std::vector<std::bitset<2>> LoadVCF(const std::string& filePath);
};

} // namespace HP_PBWT

#endif // UTL_H 