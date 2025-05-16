#ifndef PROGRAM_H
#define PROGRAM_H

#include <vector>
#include <string>
#include <queue>

namespace HP_PBWT {

    class PBWT {
    public:
        class Pal;
    };

    class Program {
    public:
        static std::vector<int> LM_Collection;
        static int nHap;
        static int nIndv;
        static std::vector<char[]> BufferArr;
        static bool semiRun;
        static PBWT::Pal pal;
        static int THD_nHapSemiRun;
        static int rtrBlockSize;
        static int rtrN_Block;
        static int nThread;
        static int LLM_Len;
        static int LC_THD;
    };

}

#endif // PROGRAM_H 