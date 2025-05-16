#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <chrono>
#include <sstream>
#include "PBWT.h"
#include "utl.h"

using namespace HP_PBWT;

// Global variables (similar to C# Program class)
int nHap = 0;
int nIndv = 0;
std::vector<char[]> BufferArr;
utl::BufferWriter BW;
int LLM_Len = 0;
int LC_THD = 0;
int rtrBlockSize = 500000000;
int rtrN_Block = 20;
int nThread = 0;

PBWT::Pal pal; // Declare an instance of PBWT::Pal

void Run(const std::string& VCF_Path, const std::string& outPath) {
    utl::RoundTableReaderV2 rtr(VCF_Path);
    BufferArr = rtr.BufferArr;
    BW = utl::BufferWriter(outPath);

    std::thread t1(&utl::RoundTableReaderV2::BlockReader, &rtr);
    std::thread t2(&utl::RoundTableReaderV2::ProbTop, &rtr);
    t1.join();
    t2.join();

    std::cout << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) << " Initializing HP-PBWT " << nHap << " Samples (#haplotypes) " << nThread << " Threads..." << std::endl;
    std::cout << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) << " HP-PBWT Initialized." << std::endl;

    std::thread t3(&utl::RoundTableReaderV2::LineParser, &rtr);
    std::thread t4(&utl::RoundTableReaderV2::LineTaker_LLM, &rtr);
    std::thread t5(&utl::BufferWriter::Run, &BW);

    t3.join();
    t4.join();
    t5.join();
}

int main(int argc, char* argv[]) {
    if (argc != 4 && argc != 7) {
        std::cerr << "Usage: " << argv[0] << " <Input VCF path> <Output path> <# of Threads> [<Long Match Length> <rtrBlockSize> <rtrN_Block>]" << std::endl;
        return 1; // Indicate an error
    }

    std::string inVCF = argv[1];
    std::string outFile = argv[2];
    nThread = std::stoi(argv[3]);

    if (argc >= 4) {
        LLM_Len = std::stoi(argv[4]);
        LC_THD = LLM_Len - 1;
    }

    if (argc == 7) {
        rtrBlockSize = std::stoi(argv[5]);
        rtrN_Block = std::stoi(argv[6]);
    }

    std::cout << "HP-PBWT" << std::endl;
    std::cout << "VCF: " << inVCF << std::endl;
    std::cout << "Output: " << outFile << std::endl;
    std::cout << "nThread: " << nThread << std::endl;
    std::cout << "Length: " << LLM_Len << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    Run(inVCF, outFile);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << std::endl << "MS: " << duration.count() << std::endl;

    return 0; // Indicate success
}