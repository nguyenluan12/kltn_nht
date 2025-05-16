#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <omp.h> // For parallel processing

namespace HP_PBWT {

    class PBWT {
    public:

        class Pal {
        public:
            class PDA {
            public:
                std::vector<int> pArr;
                int zeroCnt = 0;
                std::vector<int> mLens;

                PDA(int nHapTotal) {
                    pArr.resize(nHapTotal);
                    mLens.resize(nHapTotal);
                }
            };

            class UpLowRange {
            public:
                int UpRange;
                int LowerRange;
                UpLowRange(int up, int low) : UpRange(up), LowerRange(low) {}
            };

            static PDA* newPDA;
            static PDA* oldPDA;
            static PDA* temPDA;
            static std::vector<int> psHolder;
            static std::vector<int> minDZ;
            static std::vector<bool> minDZ_Ready;
            static std::vector<int> offsetsHolder;

            static int blockSize;
            static int msD_Wait;
            static int blockHasOneSignal;

            Pal() {
                std::cout << "Initializing HP-PBWT " << Program::nHap << " Samples (#haplotypes) " << Program::nThread << " Threads..." << std::endl;

                newPDA = new PDA(Program::nHap);
                oldPDA = new PDA(Program::nHap);

                psHolder.resize(Program::nHap);
                offsetsHolder.resize(Program::nThread);
                minDZ.resize(Program::nThread);
                minDZ_Ready.resize(Program::nThread);

                blockSize = static_cast<int>(std::ceil(static_cast<double>(Program::nHap) / static_cast<double>(Program::nThread)));
                std::cout << "HP-PBWT Initialized." << std::endl;
            }

            // Add more methods here...
            void coreP_Arr(utl::RoundTableReaderV2::HapMap& map) {
        // Step 1: Local sum
        #pragma omp parallel for
        for (int i = 0; i < Program::nThread; ++i) {
            int s = i * blockSize;
            if (s >= Program::nHap) continue;
            int e = i * blockSize + blockSize;
            if (e > Program::nHap) e = Program::nHap;

            if (map.Get_HapVal(oldPDA->pArr[s]) == '0') {
                psHolder[s] = 0;
            } else {
                psHolder[s] = 1;
            }

            for (int k = s + 1; k < e; ++k) {
                if (map.Get_HapVal(oldPDA->pArr[k]) == '0') {
                    psHolder[k] = psHolder[k - 1];
                } else {
                    psHolder[k] = psHolder[k - 1] + 1;
                }
            }

            offsetsHolder[i] = psHolder[e - 1];
        }

        // Step 2: Sequential offset handling
        for (int i = 1; i < Program::nThread; ++i) {
            offsetsHolder[i] = offsetsHolder[i - 1] + offsetsHolder[i];
        }

        // Step 3: Apply offsets
        #pragma omp parallel for
        for (int i = 1; i < Program::nThread; ++i) {
            int s = i * blockSize;
            if (s >= Program::nHap) continue;
            int e = i * blockSize + blockSize;
            if (e > Program::nHap) e = Program::nHap;

            for (int k = s; k < e; ++k) {
                psHolder[k] += offsetsHolder[i - 1];
            }
        }

        // Step 4: pps -> index settle to new p arr
        newPDA->zeroCnt = Program::nHap - psHolder[Program::nHap - 1];
        int oneOff = newPDA->zeroCnt - 1;
        #pragma omp parallel for
        for (int i = 0; i < Program::nHap; ++i) {
            if (map.Get_HapVal(oldPDA->pArr[i]) == '0') {
                newPDA->pArr[i - psHolder[i]] = oldPDA->pArr[i];
            } else {
                newPDA->pArr[psHolder[i] + oneOff] = oldPDA->pArr[i];
            }
        }
    }
    void coreD_Arr(utl::RoundTableReaderV2::HapMap& map) {
    #pragma omp parallel for
    for (int i = 0; i < Program::nThread; ++i) {
        int s = i * blockSize;
        if (s >= Program::nHap) continue;

        int e = i * blockSize + blockSize;
        if (e >= Program::nHap) e = Program::nHap;

        int prvLowM_Zero = -1;
        int prvLowM_One = -1;
        int hID = oldPDA->pArr[s];

        // First block
        if (s == 0) {
            for (int k = s; k < e; ++k) {
                hID = oldPDA->pArr[k];
                if (map.Get_HapVal(hID) == '0') {
                    newPDA->mLens[hID] = std::min(oldPDA->mLens[hID], prvLowM_One) + 1;
                    prvLowM_One = INT_MAX;
                    prvLowM_Zero = std::min(prvLowM_Zero, oldPDA->mLens[hID]);
                } else {
                    newPDA->mLens[hID] = std::min(oldPDA->mLens[hID], prvLowM_Zero) + 1;
                    prvLowM_Zero = INT_MAX;
                    prvLowM_One = std::min(prvLowM_One, oldPDA->mLens[hID]);
                }
            }
        } else {
            // Other blocks
            prvLowM_Zero = -1;
            prvLowM_One = -1;
            bool minZeroSearch = true;
            bool minOneSearch = true;

            if (psHolder[s - 1] == s) {
                prvLowM_Zero = -1;
                minZeroSearch = false;
            }

            if (psHolder[s - 1] == 0) {
                prvLowM_One = -1;
                minOneSearch = false;
            }

            if (map.Get_HapVal(oldPDA->pArr[s - 1]) == '0') {
                prvLowM_One = INT_MAX;
                minOneSearch = false;
            } else {
                prvLowM_Zero = INT_MAX;
                minZeroSearch = false;
            }

            if (minZeroSearch) {
                int seekIndex = s - 1;
                while (seekIndex >= 0 && map.Get_HapVal(oldPDA->pArr[seekIndex]) != '0') {
                    seekIndex--;
                }
                if (seekIndex >= 0) {
                    prvLowM_Zero = oldPDA->mLens[oldPDA->pArr[seekIndex]];
                    while (seekIndex >= 0 && map.Get_HapVal(oldPDA->pArr[seekIndex]) == '0') {
                        prvLowM_Zero = std::min(oldPDA->mLens[oldPDA->pArr[seekIndex]], prvLowM_Zero);
                        seekIndex--;
                    }
                    if (seekIndex == -1) {
                        prvLowM_Zero = -1;
                    }
                }
            }

            if (minOneSearch) {
                int seekIndex = s - 1;
                while (seekIndex >= 0 && map.Get_HapVal(oldPDA->pArr[seekIndex]) == '0') {
                    seekIndex--;
                }
                if (seekIndex >= 0) {
                    prvLowM_One = oldPDA->mLens[oldPDA->pArr[seekIndex]];
                    while (seekIndex >= 0 && map.Get_HapVal(oldPDA->pArr[seekIndex]) != '0') {
                        prvLowM_One = std::min(oldPDA->mLens[oldPDA->pArr[seekIndex]], prvLowM_One);
                        seekIndex--;
                    }
                    if (seekIndex == -1) {
                        prvLowM_One = -1;
                    }
                }
            }

            for (int k = s; k < e; ++k) {
                hID = oldPDA->pArr[k];
                if (map.Get_HapVal(hID) == '0') {
                    newPDA->mLens[hID] = std::min(oldPDA->mLens[hID], prvLowM_One) + 1;
                    prvLowM_One = INT_MAX;
                    prvLowM_Zero = std::min(prvLowM_Zero, oldPDA->mLens[hID]);
                } else {
                    newPDA->mLens[hID] = std::min(oldPDA->mLens[hID], prvLowM_Zero) + 1;
                    prvLowM_Zero = INT_MAX;
                    prvLowM_One = std::min(prvLowM_One, oldPDA->mLens[hID]);
                }
            }
        }
    }
}
void coreD_Arr_LowMAF(utl::RoundTableReaderV2::HapMap& map, int doneIndex) {
    ResetMinDZ();

    #pragma omp parallel for
    for (int i = 0; i < Program::nThread; ++i) {
        int s = i * blockSize;
        if (s >= Program::nHap) continue;

        int e = i * blockSize + blockSize;
        if (e >= Program::nHap) e = Program::nHap;

        int prvLowM_Zero = -1;
        int prvLowM_One = -1;
        int hID;

        int temMinDZ = INT_MAX;

        // First block
        if (s == 0) {
            if (psHolder[e - 1] != 0) {
                minDZ[i] = blockHasOneSignal;
                minDZ_Ready[i] = true;
                for (int k = s; k < e; ++k) {
                    hID = oldPDA->pArr[k];
                    if (map.Get_HapVal(hID) == '0') {
                        newPDA->mLens[hID] = std::min(oldPDA->mLens[hID], prvLowM_One) + 1;
                        prvLowM_One = INT_MAX;
                        prvLowM_Zero = std::min(prvLowM_Zero, oldPDA->mLens[hID]);
                    } else {
                        newPDA->mLens[hID] = std::min(oldPDA->mLens[hID], prvLowM_Zero) + 1;
                        prvLowM_Zero = INT_MAX;
                        prvLowM_One = std::min(prvLowM_One, oldPDA->mLens[hID]);
                    }
                }
            } else {
                newPDA->mLens[oldPDA->pArr[0]] = 0;
                for (int k = s + 1; k < e; ++k) {
                    hID = oldPDA->pArr[k];
                    newPDA->mLens[hID] = oldPDA->mLens[hID] + 1;
                }
            }
            return;
        }

        bool minZeroSearch = true;
        bool minOneSearch = true;

        // Not first block, this block has only 0s
        if (psHolder[s - 1] == psHolder[e - 1]) {
            hID = oldPDA->pArr[s];
            if (psHolder[s - 1] == s) {
                newPDA->mLens[hID] = 0;
            } else {
                prvLowM_Zero = oldPDA->mLens[oldPDA->pArr[s]];
                int seekIndex = s - 1;
                while (seekIndex >= 0 && map.Get_HapVal(oldPDA->pArr[seekIndex]) != '0') {
                    prvLowM_Zero = std::min(oldPDA->mLens[oldPDA->pArr[seekIndex]], prvLowM_Zero);
                    seekIndex--;
                }
            }
            newPDA->mLens[hID] = prvLowM_Zero + 1;
            temMinDZ = newPDA->mLens[hID];
            for (int k = s + 1; k < e; ++k) {
                hID = oldPDA->pArr[k];
                newPDA->mLens[hID] = oldPDA->mLens[hID] + 1;
                temMinDZ = std::min(temMinDZ, newPDA->mLens[hID]);
            }
            minDZ[i] = temMinDZ;
            minDZ_Ready[i] = true;
            return;
        }

        // Not first block, this block has 1s
        minDZ[i] = blockHasOneSignal;
        minDZ_Ready[i] = true;

        if (psHolder[s - 1] == s) {
            prvLowM_Zero = -1;
            minZeroSearch = false;
        }

        if (psHolder[s - 1] == 0) {
            prvLowM_One = -1;
            minOneSearch = false;
        }

        if (minZeroSearch) {
            prvLowM_Zero = oldPDA->mLens[oldPDA->pArr[s]];
            int seekIndex = s - 1;
            while (seekIndex >= 0 && map.Get_HapVal(oldPDA->pArr[seekIndex]) != '0') {
                prvLowM_Zero = std::min(oldPDA->mLens[oldPDA->pArr[seekIndex]], prvLowM_Zero);
                seekIndex--;
            }
        }

        if (minOneSearch) {
            int seekIndex = s - 1;
            int search_BlockID = i - 1;
            prvLowM_One = oldPDA->mLens[oldPDA->pArr[seekIndex]];
            while (true) {
                if (!minDZ_Ready[search_BlockID]) continue;
                if (minDZ[search_BlockID] != blockHasOneSignal) {
                    prvLowM_One = std::min(minDZ[search_BlockID] - 1, prvLowM_One);
                    search_BlockID--;
                } else {
                    break;
                }
            }
            seekIndex = search_BlockID * blockSize + blockSize - 1;
            while (map.Get_HapVal(oldPDA->pArr[seekIndex]) == '0') {
                prvLowM_One = std::min(oldPDA->mLens[oldPDA->pArr[seekIndex]], prvLowM_One);
                seekIndex--;
            }
        }

        hID = oldPDA->pArr[s];
        if (map.Get_HapVal(oldPDA->pArr[s - 1]) == map.Get_HapVal(oldPDA->pArr[s])) {
            if (map.Get_HapVal(hID) == '0') {
                newPDA->mLens[hID] = oldPDA->mLens[hID] + 1;
                prvLowM_One = INT_MAX;
                prvLowM_Zero = std::min(prvLowM_Zero, oldPDA->mLens[hID]);
            } else {
                newPDA->mLens[hID] = oldPDA->mLens[hID] + 1;
                prvLowM_Zero = INT_MAX;
                prvLowM_One = std::min(prvLowM_One, oldPDA->mLens[hID]);
            }
        } else {
            if (map.Get_HapVal(hID) == '0') {
                newPDA->mLens[hID] = std::min(prvLowM_One, oldPDA->mLens[hID]) + 1;
                prvLowM_One = INT_MAX;
                prvLowM_Zero = std::min(prvLowM_Zero, oldPDA->mLens[hID]);
            } else {
                newPDA->mLens[hID] = std::min(prvLowM_Zero, oldPDA->mLens[hID]) + 1;
                prvLowM_Zero = INT_MAX;
                prvLowM_One = std::min(prvLowM_One, oldPDA->mLens[hID]);
            }
        }

        for (int k = s + 1; k < e; ++k) {
            hID = oldPDA->pArr[k];
            if (map.Get_HapVal(hID) == '0') {
                newPDA->mLens[hID] = std::min(oldPDA->mLens[hID], prvLowM_One) + 1;
                prvLowM_One = INT_MAX;
                prvLowM_Zero = std::min(prvLowM_Zero, oldPDA->mLens[hID]);
            } else {
                newPDA->mLens[hID] = std::min(oldPDA->mLens[hID], prvLowM_Zero) + 1;
                prvLowM_Zero = INT_MAX;
                prvLowM_One = std::min(prvLowM_One, oldPDA->mLens[hID]);
            }
        }
    }
}
void coreD_Arr_BU_a(utl::RoundTableReaderV2::HapMap& map) {
    ResetMinDZ();

    #pragma omp parallel for
    for (int i = 0; i < Program::nThread; ++i) {
        int s = i * blockSize;
        if (s >= Program::nHap) continue;

        int e = i * blockSize + blockSize;
        if (e >= Program::nHap) e = Program::nHap;

        int prvLowM_Zero = -1;
        int prvLowM_One = -1;
        int hID;

        int temMinDZ = INT_MAX;

        // First block
        if (s == 0) {
            if (psHolder[e - 1] != 0) {
                minDZ[i] = blockHasOneSignal;
                minDZ_Ready[i] = true;
                for (int k = s; k < e; ++k) {
                    hID = oldPDA->pArr[k];
                    if (map.Get_HapVal(hID) == '0') {
                        newPDA->mLens[hID] = std::min(oldPDA->mLens[hID], prvLowM_One) + 1;
                        prvLowM_One = INT_MAX;
                        prvLowM_Zero = std::min(prvLowM_Zero, oldPDA->mLens[hID]);
                    } else {
                        newPDA->mLens[hID] = std::min(oldPDA->mLens[hID], prvLowM_Zero) + 1;
                        prvLowM_Zero = INT_MAX;
                        prvLowM_One = std::min(prvLowM_One, oldPDA->mLens[hID]);
                    }
                }
            } else {
                newPDA->mLens[oldPDA->pArr[0]] = 0;
                for (int k = s + 1; k < e; ++k) {
                    hID = oldPDA->pArr[k];
                    newPDA->mLens[hID] = oldPDA->mLens[hID] + 1;
                }
            }
            return;
        }

        bool minZeroSearch = true;
        bool minOneSearch = true;

        // Not first block, this block has only 0s
        if (psHolder[s - 1] == psHolder[e - 1]) {
            hID = oldPDA->pArr[s];
            if (psHolder[s - 1] == s) {
                newPDA->mLens[hID] = 0;
            } else {
                prvLowM_Zero = oldPDA->mLens[oldPDA->pArr[s]];
                int seekIndex = s - 1;
                while (seekIndex >= 0 && map.Get_HapVal(oldPDA->pArr[seekIndex]) != '0') {
                    prvLowM_Zero = std::min(oldPDA->mLens[oldPDA->pArr[seekIndex]], prvLowM_Zero);
                    seekIndex--;
                }
            }
            newPDA->mLens[hID] = prvLowM_Zero + 1;
            temMinDZ = newPDA->mLens[hID];
            for (int k = s + 1; k < e; ++k) {
                hID = oldPDA->pArr[k];
                newPDA->mLens[hID] = oldPDA->mLens[hID] + 1;
                temMinDZ = std::min(temMinDZ, newPDA->mLens[hID]);
            }
            minDZ[i] = temMinDZ;
            minDZ_Ready[i] = true;
            return;
        }

        // Not first block, this block has 1s
        minDZ[i] = blockHasOneSignal;
        minDZ_Ready[i] = true;

        int firstprvLowM_One = INT_MAX;
        int firstprvLowM_Zero = INT_MAX;

        int firstZeroIndex = -1, firstOneIndex = -1;
        bool firstOneFound = false;
        bool firstZeroFound = false;

        for (int k = s; k < e; ++k) {
            hID = oldPDA->pArr[k];
            if (map.Get_HapVal(hID) == '0') {
                if (!firstZeroFound) {
                    firstZeroFound = true;
                    firstZeroIndex = k;
                }
                if (!firstOneFound) {
                    firstprvLowM_One = std::min(oldPDA->mLens[hID], firstprvLowM_One);
                }
                newPDA->mLens[hID] = std::min(oldPDA->mLens[hID], prvLowM_One) + 1;
                prvLowM_One = INT_MAX;
                prvLowM_Zero = std::min(prvLowM_Zero, oldPDA->mLens[hID]);
            } else {
                if (!firstOneFound) {
                    firstOneFound = true;
                    firstOneIndex = k;
                }
                if (!firstZeroFound) {
                    firstprvLowM_Zero = std::min(oldPDA->mLens[hID], firstprvLowM_Zero);
                }
                newPDA->mLens[hID] = std::min(oldPDA->mLens[hID], prvLowM_Zero) + 1;
                prvLowM_Zero = INT_MAX;
                prvLowM_One = std::min(prvLowM_One, oldPDA->mLens[hID]);
            }
        }

        if (psHolder[s - 1] == s) {
            prvLowM_Zero = -1;
            minZeroSearch = false;
        }

        if (psHolder[s - 1] == 0) {
            prvLowM_One = -1;
            minOneSearch = false;
        }

        if (minZeroSearch) {
            prvLowM_Zero = oldPDA->mLens[oldPDA->pArr[s]];
            int seekIndex = s - 1;
            while (seekIndex >= 0 && map.Get_HapVal(oldPDA->pArr[seekIndex]) != '0') {
                prvLowM_Zero = std::min(oldPDA->mLens[oldPDA->pArr[seekIndex]], prvLowM_Zero);
                seekIndex--;
            }
        }

        if (minOneSearch) {
            int seekIndex = s - 1;
            int search_BlockID = i - 1;
            prvLowM_One = oldPDA->mLens[oldPDA->pArr[seekIndex]];
            while (true) {
                if (!minDZ_Ready[search_BlockID]) continue;
                if (minDZ[search_BlockID] != blockHasOneSignal) {
                    prvLowM_One = std::min(minDZ[search_BlockID] - 1, prvLowM_One);
                    search_BlockID--;
                } else {
                    break;
                }
            }
            seekIndex = search_BlockID * blockSize + blockSize - 1;
            while (map.Get_HapVal(oldPDA->pArr[seekIndex]) == '0') {
                prvLowM_One = std::min(oldPDA->mLens[oldPDA->pArr[seekIndex]], prvLowM_One);
                seekIndex--;
            }
        }

        if (firstOneIndex != -1) {
            prvLowM_One = std::min(prvLowM_One, firstprvLowM_One);
            newPDA->mLens[oldPDA->pArr[firstOneIndex]] = std::min(oldPDA->mLens[oldPDA->pArr[firstOneIndex]], prvLowM_One) + 1;
        }

        if (firstZeroIndex != -1) {
            prvLowM_Zero = std::min(prvLowM_Zero, firstprvLowM_Zero);
            newPDA->mLens[oldPDA->pArr[firstZeroIndex]] = std::min(oldPDA->mLens[oldPDA->pArr[firstZeroIndex]], prvLowM_Zero) + 1;
        }
    }
}
void ReportLongMatches_SL(int siteIndex, utl::RoundTableReaderV2::HapMap& map) {
    #pragma omp parallel for
    for (int h = 0; h < Program::nHap - 1; ++h) {
        int minVal = INT_MAX;
        for (int i = h + 1; i < Program::nHap; ++i) {
            if (oldPDA->mLens[oldPDA->pArr[i]] < Program::LLM_Len) {
                break;
            }
            minVal = std::min(minVal, oldPDA->mLens[oldPDA->pArr[i]]);
            if (map.Get_HapVal(oldPDA->pArr[h]) != map.Get_HapVal(oldPDA->pArr[i])) {
                // Report
                #pragma omp critical
                {
                    Program::BW.push_back(std::to_string(oldPDA->pArr[h]) + "\t" + std::to_string(oldPDA->pArr[i]) + "\t" + std::to_string(siteIndex) + "\t" + std::to_string(minVal));
                }
            }
        }
    }
}
void ReportLongMatches_SL_Tail(int siteIndex) {
    #pragma omp parallel for
    for (int h = 0; h < Program::nHap - 1; ++h) {
        int minVal = INT_MAX;
        for (int i = h + 1; i < Program::nHap; ++i) {
            if (oldPDA->mLens[oldPDA->pArr[i]] < Program::LLM_Len) {
                break;
            }
            minVal = std::min(minVal, oldPDA->mLens[oldPDA->pArr[i]]);
            // Report
            #pragma omp critical
            {
                Program::BW.push_back(std::to_string(oldPDA->pArr[h]) + "\t" + std::to_string(oldPDA->pArr[i]) + "\t" + std::to_string(siteIndex) + "\t" + std::to_string(minVal));
            }
        }
    }
}
void ReportSetMaxMatch(int siteIndex, utl::RoundTableReaderV2::HapMap& map) {
    #pragma omp parallel for
    for (int oneChosenIndex = 0; oneChosenIndex < Program::nHap; ++oneChosenIndex) {
        int rangeUp = oneChosenIndex, rangeDown = oneChosenIndex;
        int currHID = oldPDA->pArr[oneChosenIndex];
        int len;
        char val = map.Get_HapVal(currHID);
        int scanIndex;

        bool runThough = false;

        if (oneChosenIndex == 0) {
            len = oldPDA->mLens[oldPDA->pArr[oneChosenIndex + 1]];
        } else if (oneChosenIndex == Program::nHap - 1) {
            len = oldPDA->mLens[oldPDA->pArr[oneChosenIndex]];
        } else {
            len = std::max(oldPDA->mLens[oldPDA->pArr[oneChosenIndex]], oldPDA->mLens[oldPDA->pArr[oneChosenIndex + 1]]);
        }

        // Check down
        if (oneChosenIndex != Program::nHap - 1) {
            scanIndex = oneChosenIndex + 1;
            while (scanIndex < Program::nHap && len <= oldPDA->mLens[oldPDA->pArr[scanIndex]]) {
                if (val == map.Get_HapVal(oldPDA->pArr[scanIndex])) {
                    runThough = true;
                    break;
                }
                scanIndex++;
            }
            rangeDown = scanIndex - 1;
        }

        if (runThough) continue;

        // Check up
        if (oneChosenIndex != 0) {
            scanIndex = oneChosenIndex - 1;
            while (scanIndex >= 0 && len <= oldPDA->mLens[oldPDA->pArr[scanIndex + 1]]) {
                if (val == map.Get_HapVal(oldPDA->pArr[scanIndex])) {
                    runThough = true;
                    break;
                }
                scanIndex--;
            }
            rangeUp = scanIndex + 1;
        }

        if (runThough) continue;

        // Report
        for (int i = rangeUp; i <= rangeDown; ++i) {
            if (i == oneChosenIndex) continue;
            #pragma omp critical
            {
                Program::BW.push_back(std::to_string(currHID) + "\t" + std::to_string(oldPDA->pArr[i]) + "\t" + std::to_string(siteIndex) + " " + std::to_string(len));
            }
        }
    }
}
void ReportSetMaxMatch_Tail(int siteIndex) {
    #pragma omp parallel for
    for (int oneChosenIndex = 0; oneChosenIndex < Program::nHap; ++oneChosenIndex) {
        int rangeUp = oneChosenIndex, rangeDown = oneChosenIndex;
        int currHID = oldPDA->pArr[oneChosenIndex];
        int len;
        int scanIndex;

        if (oneChosenIndex == 0) {
            len = oldPDA->mLens[oldPDA->pArr[oneChosenIndex + 1]];
        } else if (oneChosenIndex == Program::nHap - 1) {
            len = oldPDA->mLens[oldPDA->pArr[oneChosenIndex]];
        } else {
            len = std::max(oldPDA->mLens[oldPDA->pArr[oneChosenIndex]], oldPDA->mLens[oldPDA->pArr[oneChosenIndex + 1]]);
        }

        // Check down
        if (oneChosenIndex != Program::nHap - 1) {
            scanIndex = oneChosenIndex + 1;
            while (scanIndex < Program::nHap && len <= oldPDA->mLens[oldPDA->pArr[scanIndex]]) {
                scanIndex++;
            }
            rangeDown = scanIndex - 1;
        }

        // Check up
        if (oneChosenIndex != 0) {
            scanIndex = oneChosenIndex - 1;
            while (scanIndex >= 0 && len <= oldPDA->mLens[oldPDA->pArr[scanIndex + 1]]) {
                scanIndex--;
            }
            rangeUp = scanIndex + 1;
        }

        // Report
        for (int i = rangeUp; i <= rangeDown; ++i) {
            if (i == oneChosenIndex) continue;
            #pragma omp critical
            {
                Program::BW.push_back(std::to_string(currHID) + "\t" + std::to_string(oldPDA->pArr[i]) + "\t" + std::to_string(siteIndex) + " " + std::to_string(len));
            }
        }
    }
}
void OneSort(utl::RoundTableReaderV2::HapMap& map, int doneIndex) {
    coreP_Arr(map);
    coreD_Arr(map);

    PDA* temp = oldPDA;
    oldPDA = newPDA;
    newPDA = temp;
}

void InitialSort(utl::RoundTableReaderV2::HapMap& map) {
    // Step 1: Local sum
    #pragma omp parallel for
    for (int i = 0; i < Program::nThread; ++i) {
        int s = i * blockSize;
        if (s >= Program::nHap) continue;

        int e = i * blockSize + blockSize;
        if (e > Program::nHap) e = Program::nHap;

        if (map.Get_HapVal(s) == '0') {
            psHolder[s] = 0;
        } else {
            psHolder[s] = 1;
        }

        for (int k = s + 1; k < e; ++k) {
            if (map.Get_HapVal(k) == '0') {
                psHolder[k] = psHolder[k - 1];
            } else {
                psHolder[k] = psHolder[k - 1] + 1;
            }
        }

        offsetsHolder[i] = psHolder[e - 1];
    }

    // Step 2: Sequential offset handling
    for (int i = 1; i < Program::nThread; ++i) {
        offsetsHolder[i] = offsetsHolder[i - 1] + offsetsHolder[i];
    }

    // Step 3: Apply offsets
    #pragma omp parallel for
    for (int i = 1; i < Program::nThread; ++i) {
        int s = i * blockSize;
        if (s >= Program::nHap) continue;
        int e = i * blockSize + blockSize;
        if (e > Program::nHap) e = Program::nHap;

        for (int k = s; k < e; ++k) {
            psHolder[k] += offsetsHolder[i - 1];
        }
    }

    // Step 4: pps -> index settle to new p arr
    oldPDA->zeroCnt = Program::nHap - psHolder[Program::nHap - 1];
    int oneOff = oldPDA->zeroCnt - 1;

    #pragma omp parallel for
    for (int i = 0; i < Program::nHap; ++i) {
        if (map.Get_HapVal(i) == '0') {
            oldPDA->pArr[i - psHolder[i]] = i;
            oldPDA->mLens[i - psHolder[i]] = 1;
        } else {
            oldPDA->pArr[psHolder[i] + oneOff] = i;
            oldPDA->mLens[psHolder[i] + oneOff] = 1;
        }
    }

    // D Arr adjust
    oldPDA->mLens[oldPDA->pArr[0]] = 0;
    if (oldPDA->zeroCnt != Program::nHap) {
        oldPDA->mLens[oldPDA->pArr[oneOff + 1]] = 0;
    }
}
        };
    };
}