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
#include <omp.h> // For parallel processing
#include "Program.h" // Ensure Program is included

namespace HP_PBWT {

    class utl {
    public:
        static std::string IntTo_BMK_Str(int n) {
            if (n >= 1000000000) {
                return std::to_string(n / 1000000000) + "b";
            } else if (n >= 1000000) {
                return std::to_string(n / 1000000) + "m";
            } else if (n >= 1000) {
                return std::to_string(n / 1000) + "k";
            } else {
                return std::to_string(n);
            }
        }

        class RoundTableReaderV2 {
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
            class HapMap {
            private:
                int LastHID_FirstBlk = -1;
                int charOff_SecBlk = -1;
                int secBlkID = -1;
                lineTag tag;

            public:
                HapMap(const lineTag& tag, int nHap) : tag(tag) {
                    if (!tag.dataCross) {
                        LastHID_FirstBlk = nHap - 1;
                        return;
                    }

                    LastHID_FirstBlk = (Program::rtrBlockSize - tag.charIndex - 1) / 2;
                    charOff_SecBlk = (Program::rtrBlockSize - tag.charIndex) % 2;
                    secBlkID = (tag.blockID + 1) % Program::rtrN_Block;
                }

                char Get_HapVal(int hID) {
                    if (!tag.dataCross || hID <= LastHID_FirstBlk) {
                        return Program::BufferArr[tag.blockID][tag.charIndex + hID * 2];
                    }
                    return Program::BufferArr[secBlkID][charOff_SecBlk + (hID - LastHID_FirstBlk - 1) * 2];
                }
            };

            class lineTag {
            public:
                int blockID;
                int charIndex;
                bool dataCross;
                bool headerCross;

                lineTag() : blockID(0), charIndex(0), dataCross(false), headerCross(false) {}
                lineTag(int blockID, int charIndex, bool dataCross, bool headerCross = false)
                    : blockID(blockID), charIndex(charIndex), dataCross(dataCross), headerCross(headerCross) {}
            };

            RoundTableReaderV2(const std::string& path) : VCF_Path(path) {
                // Initialize member variables
                msTakeWait = 500;
                L1_rPTR = 0;
                L1_ReadComplete = false;
                L1_Taker_Index = 0;
                L2_ReadComplete = false;
                L2_nBlock = 20;
                L1_nBlock = 20;
                
                // Initialize vectors
                L1_Stat.resize(L1_nBlock);
                L1_nRead.resize(L1_nBlock);
                BufferArr.resize(L1_nBlock);
                for (auto& buffer : BufferArr) {
                    buffer.resize(Program::rtrBlockSize);
                }
            }

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
            BufferWriter(const std::string& outPath) : outPath(outPath), AddComplete(false), msWait(500) {}

            void Add(const std::string& s) {
                buffer.push(s);
            }

            void DoneAdding() {
                AddComplete = true;
            }

            void Run() {
                sw.open(outPath);
                sw << std::endl;
                std::string sOut;
                while (!AddComplete || !buffer.empty()) {
                    if (buffer.empty()) {
                        std::this_thread::sleep_for(std::chrono::milliseconds(msWait));
                    } else {
                        sOut = buffer.front();
                        buffer.pop();
                        if (!sOut.empty()) {
                            sw << sOut << std::endl;
                        }
                    }
                }

                sw.close();
                std::cout << "Buffer Writer Complete." << std::endl;
                std::exit(0);
            }
        };

        static std::vector<std::bitset<2>> LoadVCF(const std::string& filePath);
    };

}

// Ensure all necessary variables are declared
std::vector<int> L1_Stat;
int msTakeWait;
std::vector<std::vector<char>> BufferArr;
int L1_rPTR;
int dataLineLen;
bool L1_ReadComplete;
int L1_Taker_Index;
std::queue<RoundTableReaderV2::lineTag> Lines;
int L2_nBlock;
int L1_nBlock;
bool L2_ReadComplete;
std::vector<int> L1_nRead;

// Replace starts_with and ends_with with appropriate alternatives
// For example, use rfind or find to check for prefixes or suffixes

// Correct the constructor for RoundTableReaderV2
RoundTableReaderV2(const std::string& path) : VCF_Path(path) {
    // Initialization code here
}

// Correct the constructor for HapMap
HapMap(lineTag tag, int nHap) : tag(tag) {
    if (!tag.dataCross) {
        LastHID_FirstBlk = nHap - 1;
        return;
    }

    LastHID_FirstBlk = (Program::rtrBlockSize - tag.charIndex - 1) / 2;
    charOff_SecBlk = (Program::rtrBlockSize - tag.charIndex) % 2;
    secBlkID = (tag.blockID + 1) % Program::rtrN_Block;
}

// Correct the constructor for lineTag
lineTag() : blockID(0), charIndex(0), dataCross(false), headerCross(false) {}

// Correct the use of Program::pal and Program::BW
PBWT::Pal Program::pal;
utl::BufferWriter Program::BW;

// Replace unsupported methods
bool starts_with(const std::string& str, const std::string& prefix) {
    return str.rfind(prefix, 0) == 0;
}

bool ends_with(const std::string& str, const std::string& suffix) {
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}