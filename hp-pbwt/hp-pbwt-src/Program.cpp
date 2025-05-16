#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <thread>
#include <chrono>
#include <algorithm>
#include <omp.h> // For parallel processing
#include <queue>
#include <sstream>
#include <future>
#include <list>
#include <cstring>
#include "utl.cpp" // Include the converted utl.cpp file

using namespace HP_PBWT;

class Program {
public:
    static std::ifstream dbgSR;
    static std::list<int> LM_Collection;
    static int nHap;
    static int nIndv;
    static std::vector<std::vector<char>> BufferArr;
    static bool semiRun;
    static PBWT::Pal pal;
    static utl::BufferWriter BW;
    static int THD_nHapSemiRun;
    static int rtrBlockSize;
    static int rtrN_Block;
    static int nThread;
    static int LLM_Len;
    static int LC_THD;

    static void Run(const std::string& VCF_Path, const std::string& outPath) {
        utl::RoundTableReaderV2 rtr(VCF_Path);
        BufferArr = rtr.BufferArr;

        BW = utl::BufferWriter(outPath);

        auto t1 = std::async(std::launch::async, [&rtr]() { rtr.BlockReader(); });

        rtr.ProbTop();

        pal = PBWT::Pal();

        auto t2 = std::async(std::launch::async, [&rtr]() { rtr.LineParser(); });

        std::future<void> t3;
        if (LLM_Len < 1) {
            t3 = std::async(std::launch::async, [&rtr]() { rtr.LineTaker_SMM(); });
        } else {
            t3 = std::async(std::launch::async, [&rtr]() { rtr.LineTaker_LLM(); });
        }

        auto t4 = std::async(std::launch::async, []() { BW.Run(); });

        t1.get();
        t2.get();
        t3.get();
        t4.get();
    }

    static int main(int argc, char* argv[]) {
        if (argc != 5 && argc != 8) {
            std::cout << "HP-PBWT" << std::endl;
            std::cout << "<Input VCF path,string> <Output path,string> <# of Thread,int> <L-Long Match Length,int. <1 for SetMaxMatch.>" << std::endl;
            std::cout << "<Input VCF path,string> <Output path,string> <# of Thread,int> <L-Long Match Length,int. <1 for SetMaxMatch.> --SetBuffer <Buffer Size> <# of Buffer>" << std::endl;
            return 1;
        }

        std::string inVCF = argv[1];
        std::string outFile = argv[2];
        nThread = std::stoi(argv[3]);
        LLM_Len = std::stoi(argv[4]);

        std::cout << "HP-PBWT" << std::endl;
        std::cout << "VCF: " << inVCF << std::endl;
        std::cout << "Output: " << outFile << std::endl;
        std::cout << "nThread: " << nThread << std::endl;
        std::cout << "Length: " << LLM_Len << std::endl;

        if (argc == 8) {
            rtrBlockSize = std::stoi(argv[6]);
            rtrN_Block = std::stoi(argv[7]);
        }

        LC_THD = LLM_Len - 1;

        auto start = std::chrono::high_resolution_clock::now();

        Run(inVCF, outFile);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        std::cout << "MS: " << elapsed.count() << std::endl;

        return 0;
    }
};

// Initialize static members
std::ifstream Program::dbgSR;
std::list<int> Program::LM_Collection;
int Program::nHap = 0;
int Program::nIndv = 0;
std::vector<std::vector<char>> Program::BufferArr;
bool Program::semiRun = false;
PBWT::Pal Program::pal;
utl::BufferWriter Program::BW("");
int Program::THD_nHapSemiRun = 10000;
int Program::rtrBlockSize = 500000000;
int Program::rtrN_Block = 20;
int Program::nThread = 0;
int Program::LLM_Len = 0;
int Program::LC_THD = 0;