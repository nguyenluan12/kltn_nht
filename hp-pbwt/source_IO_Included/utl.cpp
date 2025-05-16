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

        void ProbTop() {
            while (L1_Stat[0] != 0) {
                std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
                continue;
            }

            // Skip top headers
            while (BufferArr[0][L1_rPTR] == '#' && BufferArr[0][L1_rPTR + 1] == '#') {
                while (BufferArr[0][L1_rPTR] != '\n') {
                    L1_rPTR++;
                }
                L1_rPTR++;
            }

            // Count nIndv
            Program::nIndv = 0;
            while (BufferArr[0][L1_rPTR] != '\n') {
                if (BufferArr[0][L1_rPTR] == '\t') {
                    Program::nIndv++;
                }
                L1_rPTR++;
            }
            Program::nIndv -= 8;
            Program::nHap = Program::nIndv * 2;
            dataLineLen = Program::nIndv * 4;
            std::cout << "nIndv. " << Program::nIndv << " dataLineLen " << dataLineLen << std::endl;
            L1_rPTR++;
        }

        void LineParser() {
            bool headerCross = false;
            int tabCnt = 0;
            int nextBlock;
            int i;

            while (!L1_ReadComplete || L1_Stat[L1_Taker_Index] == 0) {
                if (L1_Stat[L1_Taker_Index] != 0 || Lines.size() > L2_nBlock) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
                    continue;
                }

                nextBlock = (L1_Taker_Index + 1) % L1_nBlock;

                for (i = L1_rPTR; i < L1_nRead[L1_Taker_Index]; i++) {
                    if (BufferArr[L1_Taker_Index][i] == '\t') {
                        tabCnt++;

                        if (tabCnt == 9) {
                            tabCnt = 0;
                            i++;
                            if (i + dataLineLen >= L1_nRead[L1_Taker_Index]) {
                                L1_rPTR = (i + dataLineLen) % L1_nRead[L1_Taker_Index];

                                while (L1_Stat[nextBlock] != 0) {
                                    std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
                                    continue;
                                }

                                Lines.push(lineTag(L1_Taker_Index, i, true, false));
                                break;
                            } else {
                                if (headerCross) {
                                    Lines.push(lineTag(L1_Taker_Index, i, false, true));
                                    headerCross = false;
                                } else {
                                    Lines.push(lineTag(L1_Taker_Index, i, false, false));
                                }
                                i += dataLineLen;
                            }
                        }
                    }
                }

                while (L1_Stat[nextBlock] != 0) {
                    if (L1_ReadComplete) {
                        L2_ReadComplete = true;
                        std::cout << "L2 Read Complete." << std::endl;
                        return;
                    }
                    std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
                    continue;
                }

                if (tabCnt != 0) {
                    L1_rPTR = 0;
                    headerCross = true;
                } else {
                    L1_Stat[L1_Taker_Index] = -2;
                }

                L1_Taker_Index = nextBlock;
            }

            L2_ReadComplete = true;
            std::cout << "L2 Read Complete." << std::endl;
        }

        void LineTaker_LLM() {
            lineTag tag;
            bool dq;
            int doneIndex = 0;

            while (!L2_ReadComplete || !Lines.empty()) {
                while (Lines.empty()) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
                    continue;
                }

                dq = !Lines.empty();
                if (!dq) continue;

                tag = Lines.front();
                Lines.pop();

                HapMap map(tag, Program::nHap);
                Program::pal.InitialSort(map);
                break;
            }

            while (!L2_ReadComplete || !Lines.empty()) {
                while (Lines.empty()) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
                    continue;
                }

                dq = !Lines.empty();
                if (!dq) continue;

                tag = Lines.front();
                Lines.pop();

                LineWork_LLM(doneIndex, tag);

                if (tag.dataCross) {
                    L1_Stat[tag.blockID] = -1;
                }

                if (tag.headerCross) {
                    L1_Stat[(tag.blockID - 1 + L1_nBlock) % L1_nBlock] = -1;
                }

                doneIndex++;
            }

            std::cout << "L3 Line Taker Complete." << std::endl;
            Program::pal.ReportLongMatches_SL_Tail(doneIndex);
            Program::BW.DoneAdding();
        }

        void LineTaker_SMM() {
            lineTag tag;
            bool dq;
            int siteCnt = 1;

            while (!L2_ReadComplete || !Lines.empty()) {
                while (Lines.empty()) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
                    continue;
                }

                dq = !Lines.empty();
                if (!dq) continue;

                tag = Lines.front();
                Lines.pop();

                HapMap map(tag, Program::nHap);
                Program::pal.InitialSort(map);
                break;
            }

            siteCnt++;

            while (!L2_ReadComplete || !Lines.empty()) {
                while (Lines.empty()) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
                    continue;
                }

                dq = !Lines.empty();
                if (!dq) continue;

                tag = Lines.front();
                Lines.pop();

                LineWork_SMM(siteCnt, tag);

                if (tag.dataCross) {
                    L1_Stat[tag.blockID] = -1;
                }

                if (tag.headerCross) {
                    L1_Stat[(tag.blockID - 1 + L1_nBlock) % L1_nBlock] = -1;
                }

                siteCnt++;
            }

            std::cout << "L3 Line Taker Complete." << std::endl;
            Program::pal.ReportSetMaxMatch_Tail(siteCnt);
            Program::BW.DoneAdding();
        }

        void LineTaker_DBG() {
            std::string line;

            while (std::getline(Program::dbgSR, line) && line.starts_with("##")) {
                continue;
            }

            lineTag tag;
            bool dq;
            int cnt = 0;

            while (!L2_ReadComplete || !Lines.empty()) {
                while (Lines.empty()) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
                    continue;
                }

                dq = !Lines.empty();
                if (!dq) continue;

                tag = Lines.front();
                Lines.pop();

                std::cout << "Line " << cnt << std::endl;

                if (cnt == 9830) {
                    int a = 0; // Debugging breakpoint
                }

                LineWork_DBG(tag);

                cnt++;

                if (tag.dataCross) {
                    L1_Stat[tag.blockID] = -1;
                }

                if (tag.headerCross) {
                    L1_Stat[(tag.blockID - 1 + L1_nBlock) % L1_nBlock] = -1;
                }
            }

            std::cout << "L3 Line Taker Complete." << std::endl;
            Program::dbgSR.close();
        }

        static std::vector<std::bitset<2>> LoadVCF(const std::string& filePath) {
            std::cout << "Loading panel from VCF: " << filePath << std::endl;

            std::ifstream reader;
            if (filePath.ends_with(".vcf.gz")) {
                // Handle gzipped VCF files
                // Note: Requires additional library support for gzip handling in C++
                // This is a placeholder for actual gzip handling code
                std::cerr << "Gzipped VCF handling not implemented." << std::endl;
                return {};
            } else {
                reader.open(filePath);
            }

            std::vector<std::bitset<2>> panel;
            int numSamples = 0;
            bool headerParsed = false;
            std::string line;

            while (std::getline(reader, line)) {
                if (line.starts_with("##")) continue;

                if (line.starts_with("#CHROM")) {
                    std::istringstream iss(line);
                    std::string token;
                    int tokenCount = 0;
                    while (std::getline(iss, token, '\t')) {
                        tokenCount++;
                    }
                    numSamples = (tokenCount - 9) * 2; // each sample has 2 haplotypes
                    headerParsed = true;
                    continue;
                }

                if (!headerParsed) continue;

                std::istringstream iss(line);
                std::vector<std::string> fields;
                std::string field;
                while (std::getline(iss, field, '\t')) {
                    fields.push_back(field);
                }

                std::bitset<2> site;
                for (size_t i = 9; i < fields.size(); ++i) {
                    if (fields[i].length() >= 3) {
                        char a = fields[i][0];
                        char b = fields[i][2];
                        site[2 * (i - 9)] = (a == '1');
                        site[2 * (i - 9) + 1] = (b == '1');
                    }
                }

                panel.push_back(site);
            }

            reader.close();
            std::cout << "VCF loaded with " << panel.size() << " sites and " << numSamples << " haplotypes." << std::endl;
            return panel;
        }

    private:
        void LineWork_LLM(int doneIndex, lineTag tag) {
            HapMap map(tag, Program::nHap);
            if (doneIndex >= Program::LC_THD) {
                Program::pal.ReportLongMatches_SL(doneIndex, map);
            }
            Program::pal.OneSort(map, doneIndex);
        }

        void LineWork_SMM(int lineCnt, lineTag tag) {
            HapMap map(tag, Program::nHap);

            if (lineCnt >= Program::LLM_Len) {
                Program::pal.ReportSetMaxMatch(lineCnt, map);
            }

            Program::pal.OneSort(map, 0);
        }

        void LineWork_DBG(lineTag tag) {
            HapMap map(tag, Program::nHap);

            std::string line;
            std::getline(Program::dbgSR, line);
            std::istringstream iss(line);
            std::vector<std::string> parts;
            std::string part;
            while (std::getline(iss, part, '\t')) {
                parts.push_back(part);
            }

            std::vector<char> hapChars;
            for (size_t i = 9; i < parts.size(); ++i) {
                hapChars.push_back(parts[i][0]);
                hapChars.push_back(parts[i][2]);
            }

            #pragma omp parallel for
            for (size_t i = 0; i < hapChars.size(); ++i) {
                if (hapChars[i] != map.Get_HapVal(i)) {
                    std::cout << i << ": " << hapChars[i] << " vs " << map.Get_HapVal(i) << std::endl;
                    char c = map.Get_HapVal(i);
                    int a = 0; // Debugging breakpoint
                    c = map.Get_HapVal(i);
                }
            }
        }

        std::string VCF_Path;
        std::queue<lineTag> Lines;

        int ReadLine_BlkIndex1;
        int ReadLine_BlkIndex2;
        int Blk1_rPTR;
        int Blk2_rPTR;

        int L1_blockSize = Program::rtrBlockSize;
        int L1_nBlock = 20;
        int L2_nBlock = 1000;
        std::vector<int> L1_nRead;
        std::vector<int> L1_Stat;
        std::vector<int> L1_Seq;
        int L1_rPTR = 0;
        int L1_Add_Index = 0;
        int L1_Taker_Index = 0;
        std::vector<std::vector<char>> BufferArr;

        int msAddWait = 200;
        int msTakeWait = 100;
        bool L1_ReadComplete = false;
        bool L2_ReadComplete = false;

        int dataLineLen;
    };

    class BufferWriter {
    public:
        std::queue<std::string> buffer;
        std::ofstream sw;
        std::string outPath;
        bool AddComplete = false;
        int msWait = 500;

        BufferWriter(const std::string& outPath) : outPath(outPath) {}

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

} 