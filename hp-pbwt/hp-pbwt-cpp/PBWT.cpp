#include "PBWT.h"
#include <iostream>
#include <algorithm>
#include <execution> // If using parallel algorithms

namespace HP_PBWT {

    // Initialize static members of Pal::PDA
    PBWT::Pal::PDA PBWT::Pal::newPDA(0); // Dummy initialization, will be overwritten
    PBWT::Pal::PDA PBWT::Pal::oldPDA(0); // Dummy initialization, will be overwritten
    PBWT::Pal::PDA PBWT::Pal::temPDA(0); // Dummy initialization, will be overwritten
    std::vector<int> PBWT::Pal::psHolder;
    std::vector<int> PBWT::Pal::minDZ;
    std::vector<bool> PBWT::Pal::minDZ_Ready;
    std::vector<int> PBWT::Pal::offsetsHolder;
    int PBWT::Pal::blockSize = 0;
    //std::execution::parallel_policy PBWT::Pal::pOp; // C++17 Parallel Algorithms
    int PBWT::Pal::msD_Wait = 10;
    int PBWT::Pal::blockHasOneSignal = -10;

    // Pal::PDA Constructor
    PBWT::Pal::PDA::PDA(int nHapTotal) : pArr(nHapTotal), mLens(nHapTotal) {}

    PBWT::Pal::UpLowRange::UpLowRange(int up, int low) : UpRange(up), LowerRange(low) {}

    // Pal Constructor
    PBWT::Pal::Pal() {
        std::cout << std::chrono::system_clock::now().time_since_epoch().count() << " Initializing HP-PBWT " << ::nHap << " Samples (#haplotypes) " << ::nThread << " Threads..." << std::endl;

        newPDA = PDA(::nHap);
        oldPDA = PDA(::nHap);

        psHolder.resize(::nHap);
        offsetsHolder.resize(::nThread);
        minDZ.resize(::nThread);
        minDZ_Ready.resize(::nThread, false); // Initialize all to false

        //pOp.max_degree_of_parallelism() = Program.nThread; // C++17 way to set parallelism
        blockSize = static_cast<int>(std::ceil(static_cast<double>(::nHap) / ::nThread));

        std::cout << std::chrono::system_clock::now().time_since_epoch().count() << " HP-PBWT Initialized." << std::endl;
    }

    

    void PBWT::Pal::coreD_Arr(utl::RoundTableReaderV2::HapMap& map) {
    std::for_each(std::execution::par, std::execution::seq,
        [](int i) {
            int s = i * blockSize;
            if (s >= ::nHap) return;
            int e = i * blockSize + blockSize;
            if (e >= ::nHap) e = ::nHap;

            int prvLowM_Zero = -1;
            int prvLowM_One = -1;
            int hID;

            if (s == 0) {
                for (int k = s; k < e; k++) {
                    hID = oldPDA.pArr[k];
                    if (map.Get_HapVal(hID) == '0') {
                        newPDA.mLens[hID] = std::min(oldPDA.mLens[hID], prvLowM_One) + 1;
                        prvLowM_One = std::numeric_limits<int>::max();
                        prvLowM_Zero = std::min(prvLowM_Zero, oldPDA.mLens[hID]);
                    } else {
                        newPDA.mLens[hID] = std::min(oldPDA.mLens[hID], prvLowM_Zero) + 1;
                        prvLowM_Zero = std::numeric_limits<int>::max();
                        prvLowM_One = std::min(prvLowM_One, oldPDA.mLens[hID]);
                    }
                }
            } else {
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

                if (map.Get_HapVal(oldPDA.pArr[s - 1]) == '0') {
                    prvLowM_One = std::numeric_limits<int>::max();
                    minOneSearch = false;
                } else {
                    prvLowM_Zero = std::numeric_limits<int>::max();
                    minZeroSearch = false;
                }

                if (minZeroSearch) {
                    int seekIndex = s - 1;
                    while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) != '0') {
                        seekIndex--;
                    }
                    if (seekIndex >= 0) {
                        prvLowM_Zero = oldPDA.mLens[oldPDA.pArr[seekIndex]];
                        while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) == '0') {
                            prvLowM_Zero = std::min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_Zero);
                            seekIndex--;
                        }
                        if (seekIndex == -1) {
                            prvLowM_Zero = -1;
                        }
                    }
                }

                if (minOneSearch) {
                    int seekIndex = s - 1;
                    while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) != '0') {
                        seekIndex--;
                    }
                    if (seekIndex >= 0) {
                        prvLowM_One = oldPDA.mLens[oldPDA.pArr[seekIndex]];
                        while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) != '0') {
                            prvLowM_One = std::min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_One);
                            seekIndex--;
                        }
                        if (seekIndex == -1) {
                            prvLowM_One = -1;
                        }
                    }
                }

                for (int k = s; k < e; k++) {
                    hID = oldPDA.pArr[k];
                    if (map.Get_HapVal(hID) == '0') {
                        newPDA.mLens[hID] = std::min(oldPDA.mLens[hID], prvLowM_One) + 1;
                        prvLowM_One = std::numeric_limits<int>::max();
                        prvLowM_Zero = std::min(prvLowM_Zero, oldPDA.mLens[hID]);
                    } else {
                        newPDA.mLens[hID] = std::min(oldPDA.mLens[hID], prvLowM_Zero) + 1;
                        prvLowM_Zero = std::numeric_limits<int>::max();
                        prvLowM_One = std::min(prvLowM_One, oldPDA.mLens[hID]);
                    }
                }
            }
        }, 0, ::nThread);
        std::swap(oldPDA, newPDA);
    }

    void PBWT::Pal::ResetMinDZ() {
        for (int i = 0; i < ::nThread; i++) {
            minDZ[i] = std::numeric_limits<int>::max();
            minDZ_Ready[i] = false;
        }
    }

    void PBWT::Pal::coreD_Arr_LowMAF(utl::RoundTableReaderV2::HapMap& map, int doneIndex) {
        ResetMinDZ();

        std::for_each(std::execution::par, std::execution::seq,
            [](int i) {
                int s = i * blockSize;
                if (s >= ::nHap) return;
                int e = i * blockSize + blockSize;
                if (e >= ::nHap) e = ::nHap;

                int prvLowM_Zero = -1;
                int prvLowM_One = -1;
                int hID;

                int temMinDZ = std::numeric_limits<int>::max();

                if (s == 0) {
                    if (psHolder[e - 1] != 0) {
                        minDZ[i] = blockHasOneSignal;
                        minDZ_Ready[i] = true;
                        for (int k = s; k < e; k++) {
                            hID = oldPDA.pArr[k];
                            if (map.Get_HapVal(hID) == '0') {
                                newPDA.mLens[hID] = std::min(oldPDA.mLens[hID], prvLowM_One) + 1;
                                prvLowM_One = std::numeric_limits<int>::max();
                                prvLowM_Zero = std::min(prvLowM_Zero, oldPDA.mLens[hID]);
                            } else {
                                newPDA.mLens[hID] = std::min(oldPDA.mLens[hID], prvLowM_Zero) + 1;
                                prvLowM_Zero = std::numeric_limits<int>::max();
                                prvLowM_One = std::min(prvLowM_One, oldPDA.mLens[hID]);
                            }
                        }
                    } else {
                        newPDA.mLens[oldPDA.pArr[0]] = 0;
                        for (int k = s + 1; k < e; k++) {
                            hID = oldPDA.pArr[k];
                            newPDA.mLens[hID] = oldPDA.mLens[hID] + 1;
                        }
                    }
                    return;
                }

                bool minZeroSearch = true;
                bool minOneSearch = true;

                if (psHolder[s - 1] == psHolder[e - 1]) {
                    hID = oldPDA.pArr[s];
                    if (psHolder[s - 1] == s) {
                        newPDA.mLens[hID] = 0;
                    } else {
                        prvLowM_Zero = oldPDA.mLens[oldPDA.pArr[s]];
                        int seekIndex = s - 1;
                        while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) != '0') {
                            prvLowM_Zero = std::min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_Zero);
                            seekIndex--;
                        }
                    }
                    newPDA.mLens[hID] = prvLowM_Zero + 1;
                    temMinDZ = newPDA.mLens[hID];
                    for (int k = s + 1; k < e; k++) {
                        hID = oldPDA.pArr[k];
                        newPDA.mLens[hID] = newPDA.mLens[hID] + 1;
                        temMinDZ = std::min(temMinDZ, newPDA.mLens[hID]);
                    }

                    minDZ[i] = temMinDZ;
                    minDZ_Ready[i] = true;
                    return;
                }

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
                    prvLowM_Zero = oldPDA.mLens[oldPDA.pArr[s]];
                    int seekIndex = s - 1;
                    while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) != '0') {
                        prvLowM_Zero = std::min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_Zero);
                        seekIndex--;
                    }
                }

                if (minOneSearch) {
                    int seekIndex = s - 1;
                    int search_BlockID = i - 1;

                    prvLowM_One = oldPDA.mLens[oldPDA.pArr[seekIndex]];

                    while (true) {
                        if (minDZ_Ready[search_BlockID] == false) {
                            continue;
                        }
                        if (minDZ[search_BlockID] != blockHasOneSignal) {
                            prvLowM_One = std::min(minDZ[search_BlockID] - 1, prvLowM_One);
                            search_BlockID--;
                        } else {
                            break;
                        }
                        if (search_BlockID < 0) break; // Added boundary check
                    }
                    if (search_BlockID >= 0) { // Added boundary check
                        seekIndex = search_BlockID * blockSize + blockSize - 1;
                        while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) == '0') {
                            prvLowM_One = std::min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_One);
                            seekIndex--;
                            if (seekIndex < 0) break; // Added boundary check
                        }
                    }
                }

                hID = oldPDA.pArr[s];
                if (s > 0 && map.Get_HapVal(oldPDA.pArr[s - 1]) == map.Get_HapVal(oldPDA.pArr[s])) {
                    if (map.Get_HapVal(hID) == '0') {
                        newPDA.mLens[hID] = oldPDA.mLens[hID] + 1;
                        prvLowM_One = std::numeric_limits<int>::max();
                        prvLowM_Zero = std::min(prvLowM_Zero, oldPDA.mLens[hID]);
                    } else {
                        newPDA.mLens[hID] = oldPDA.mLens[hID] + 1;
                        prvLowM_Zero = std::numeric_limits<int>::max();
                        prvLowM_One = std::min(prvLowM_One, oldPDA.mLens[hID]);
                    }
                } else {
                    if (map.Get_HapVal(hID) == '0') {
                        newPDA.mLens[hID] = std::min(prvLowM_One, oldPDA.mLens[hID]) + 1;
                        prvLowM_One = std::numeric_limits<int>::max();
                        prvLowM_Zero = std::min(prvLowM_Zero, oldPDA.mLens[hID]);
                    } else {
                        newPDA.mLens[hID] = std::min(prvLowM_Zero, oldPDA.mLens[hID]) + 1;
                        prvLowM_Zero = std::numeric_limits<int>::max();
                        prvLowM_One = std::min(prvLowM_One, oldPDA.mLens[hID]);
                    }
                }

                for (int k = s + 1; k < e; k++) {
                    hID = oldPDA.pArr[k];
                    if (map.Get_HapVal(hID) == '0') {
                        newPDA.mLens[hID] = std::min(oldPDA.mLens[hID], prvLowM_One) + 1;
                        prvLowM_One = std::numeric_limits<int>::max();
                        prvLowM_Zero = std::min(prvLowM_Zero, oldPDA.mLens[hID]);
                    } else {
                        newPDA.mLens[hID] = std::min(oldPDA.mLens[hID], prvLowM_Zero) + 1;
                        prvLowM_Zero = std::numeric_limits<int>::max();
                        prvLowM_One = std::min(prvLowM_One, oldPDA.mLens[hID]);
                    }
                }
            }, 0, ::nThread);
    }

    void PBWT::Pal::coreD_Arr_BU_a(utl::RoundTableReaderV2::HapMap& map) {
        ResetMinDZ();

        std::for_each(std::execution::par, std::execution::seq,
            [](int i) {
                int s = i * blockSize;
                if (s >= ::nHap) return;
                int e = i * blockSize + blockSize;
                if (e >= ::nHap) e = ::nHap;

                int prvLowM_Zero = -1;
                int prvLowM_One = -1;
                int hID;

                int temMinDZ = std::numeric_limits<int>::max();

                if (s == 0) {
                    if (psHolder[e - 1] != 0) {
                        minDZ[i] = blockHasOneSignal;
                        minDZ_Ready[i] = true;
                        for (int k = s; k < e; k++) {
                            hID = oldPDA.pArr[k];
                            if (map.Get_HapVal(hID) == '0') {
                                newPDA.mLens[hID] = std::min(oldPDA.mLens[hID], prvLowM_One) + 1;
                                prvLowM_One = std::numeric_limits<int>::max();
                                prvLowM_Zero = std::min(prvLowM_Zero, oldPDA.mLens[hID]);
                            } else {
                                newPDA.mLens[hID] = std::min(oldPDA.mLens[hID], prvLowM_Zero) + 1;
                                prvLowM_Zero = std::numeric_limits<int>::max();
                                prvLowM_One = std::min(prvLowM_One, oldPDA.mLens[hID]);
                            }
                        }
                    } else {
                        newPDA.mLens[oldPDA.pArr[0]] = 0;
                        for (int k = s + 1; k < e; k++) {
                            hID = oldPDA.pArr[k];
                            newPDA.mLens[hID] = newPDA.mLens[hID] + 1;
                        }
                    }
                    return;
                }

                bool minZeroSearch = true;
                bool minOneSearch = true;

                if (psHolder[s - 1] == psHolder[e - 1]) {
                    hID = oldPDA.pArr[s];
                    if (psHolder[s - 1] == s) {
                        newPDA.mLens[hID] = 0;
                    } else {
                        prvLowM_Zero = oldPDA.mLens[oldPDA.pArr[s]];
                        int seekIndex = s - 1;
                        while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) != '0') {
                            prvLowM_Zero = std::min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_Zero);
                            seekIndex--;
                        }
                    }
                    newPDA.mLens[hID] = prvLowM_Zero + 1;
                    temMinDZ = newPDA.mLens[hID];
                    for (int k = s + 1; k < e; k++) {
                        hID = oldPDA.pArr[k];
                        newPDA.mLens[hID] = newPDA.mLens[hID] + 1;
                        temMinDZ = std::min(temMinDZ, newPDA.mLens[hID]);
                    }

                    minDZ[i] = temMinDZ;
                    minDZ_Ready[i] = true;
                    return;
                }

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
                    prvLowM_Zero = oldPDA.mLens[oldPDA.pArr[s]];
                    int seekIndex = s - 1;
                    while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) != '0') {
                        prvLowM_Zero = std::min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_Zero);
                        seekIndex--;
                    }
                }

                if (minOneSearch) {
                    int seekIndex = s - 1;
                    int search_BlockID = i - 1;

                    prvLowM_One = oldPDA.mLens[oldPDA.pArr[seekIndex]];

                    while (true) {
                        if (minDZ_Ready[search_BlockID] == false) {
                            continue;
                        }
                        if (minDZ[search_BlockID] != blockHasOneSignal) {
                            prvLowM_One = std::min(minDZ[search_BlockID] - 1, prvLowM_One);
                            search_BlockID--;
                        } else {
                            break;
                        }
                        if (search_BlockID < 0) break;
                    }
                    if (search_BlockID >= 0) {
                        seekIndex = search_BlockID * blockSize + blockSize - 1;
                        while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) == '0') {
                            prvLowM_One = std::min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_One);
                            seekIndex--;
                            if (seekIndex < 0) break;
                        }
                    }
                }

                hID = oldPDA.pArr[s];
                if (s > 0 && map.Get_HapVal(oldPDA.pArr[s - 1]) == map.Get_HapVal(oldPDA.pArr[s])) {
                    if (map.Get_HapVal(hID) == '0') {
                        newPDA.mLens[hID] = oldPDA.mLens[hID] + 1;
                        prvLowM_One = std::numeric_limits<int>::max();
                        prvLowM_Zero = std::min(prvLowM_Zero, oldPDA.mLens[hID]);
                    } else {
                        newPDA.mLens[hID] = oldPDA.mLens[hID] + 1;
                        prvLowM_Zero = std::numeric_limits<int>::max();
                        prvLowM_One = std::numeric_limits<int>::max();
                    }
                } else {
                    if (map.Get_HapVal(hID) == '0') {
                        newPDA.mLens[hID] = std::min(prvLowM_One, oldPDA.mLens[hID]) + 1;
                        prvLowM_One = std::numeric_limits<int>::max();
                        prvLowM_Zero = std::min(prvLowM_Zero, oldPDA.mLens[hID]);
                    } else {
                        newPDA.mLens[hID] = std::min(prvLowM_Zero, oldPDA.mLens[hID]) + 1;
                        prvLowM_Zero = std::numeric_limits<int>::max();
                        prvLowM_One = std::min(prvLowM_One, oldPDA.mLens[hID]);
                    }
                }

                for (int k = s + 1; k < e; k++) {
                    hID = oldPDA.pArr[k];
                    if (map.Get_HapVal(hID) == '0') {
                        newPDA.mLens[hID] = std::min(oldPDA.mLens[hID], prvLowM_One) + 1;
                        prvLowM_One = std::numeric_limits<int>::max();
                        prvLowM_Zero = std::min(prvLowM_Zero, oldPDA.mLens[hID]);
                    } else {
                        newPDA.mLens[hID] = std::min(oldPDA.mLens[hID], prvLowM_Zero) + 1;
                        prvLowM_Zero = std::numeric_limits<int>::max();
                        prvLowM_One = std::min(prvLowM_One, oldPDA.mLens[hID]);
                    }
                }
            }, 0, ::nThread);
    }

    void PBWT::Pal::InitialSort(utl::RoundTableReaderV2::HapMap& map) {
        std::vector<int> local_offsets(::nThread);
    
        std::for_each(std::execution::par, std::execution::seq,
            [&](int i) {
                int s = i * blockSize;
                if (s >= ::nHap) return;
                int e = i * blockSize + blockSize;
                if (e > ::nHap) e = ::nHap;
    
                int zeroCnt = 0;
                for (int k = s; k < e; k++) {
                    if (map.Get_HapVal(k) == '0') { // Changed from oldPDA.pArr[k] to k
                        zeroCnt++;
                    }
                }
                local_offsets[i] = zeroCnt;
            }, 0, ::nThread);
    
        // Calculate sequential offsets
        for (size_t i = 1; i < local_offsets.size(); ++i) {
            local_offsets[i] += local_offsets[i - 1];
        }
    
        std::for_each(std::execution::par, std::execution::seq,
            [&](int i) {
                int s = i * blockSize;
                if (s >= ::nHap) return;
                int e = i * blockSize + blockSize;
                if (e > ::nHap) e = ::nHap;
    
                int zeroWrite = s;
                int oneWrite = s + local_offsets[i] - local_offsets[0];
                if (i > 0) {
                    oneWrite -= local_offsets[0];
                }
    
                for (int k = s; k < e; k++) {
                    if (map.Get_HapVal(k) == '0') { // Changed from oldPDA.pArr[k] to k
                        oldPDA.pArr[zeroWrite++] = k;  // Store k, not oldPDA.pArr[k]
                        oldPDA.mLens[k] = 1;
                    } else {
                        oldPDA.pArr[oneWrite++] = k;  // Store k, not oldPDA.pArr[k]
                        oldPDA.mLens[k] = 1;
                    }
                }
            }, 0, ::nHap);
    
        oldPDA.mLens[oldPDA.pArr[0]] = 0;
        if (oldPDA.zeroCnt != ::nHap) {
            oldPDA.mLens[oldPDA.pArr[oldPDA.zeroCnt]] = 0;
        }
    }

    void PBWT::Pal::Swap_new_old() {
        temPDA = oldPDA;
        oldPDA = newPDA;
        newPDA = temPDA;
    }
    void PBWT::Pal::LineWork_LLM(int doneIndex, utl::RoundTableReaderV2::lineTag& tag) {
        utl::RoundTableReaderV2::HapMap map(tag, ::nHap);

        ResetMinDZ();

        std::for_each(std::execution::par, std::execution::seq,
            [&](int i) {
                int s = i * blockSize;
                if (s >= ::nHap) return;
                int e = i * blockSize + blockSize;
                if (e >= ::nHap) e = ::nHap;

                int prvLowM_Zero = -1;
                int prvLowM_One = -1;
                int hID;
                int temMinDZ = std::numeric_limits<int>::max();

                if (s == 0) {
                    // (Xử lý khối đầu tiên - tương tự như coreD_Arr_LowMAF)
                     if (psHolder[e - 1] != 0) {
                         minDZ[i] = blockHasOneSignal;
                         minDZ_Ready[i] = true;
                         for (int k = s; k < e; k++) {
                             hID = oldPDA.pArr[k];
                             if (map.Get_HapVal(hID) == '0') {
                                 newPDA.mLens[hID] = std::min(oldPDA.mLens[hID], prvLowM_One) + 1;
                                 prvLowM_One = std::numeric_limits<int>::max();
                                 prvLowM_Zero = std::min(prvLowM_Zero, oldPDA.mLens[hID]);
                             } else {
                                 newPDA.mLens[hID] = std::min(oldPDA.mLens[hID], prvLowM_Zero) + 1;
                                 prvLowM_Zero = std::numeric_limits<int>::max();
                                 prvLowM_One = std::min(prvLowM_One, oldPDA.mLens[hID]);
                             }
                         }
                     } else {
                         newPDA.mLens[oldPDA.pArr[0]] = 0;
                         for (int k = s + 1; k < e; k++) {
                             hID = oldPDA.pArr[k];
                             newPDA.mLens[hID] = newPDA.mLens[hID] + 1;
                         }
                     }
                    return;
                }

                bool minZeroSearch = true;
                bool minOneSearch = true;

                if (psHolder[s - 1] == psHolder[e - 1]) {
                    // (Xử lý khối có cùng giá trị)
                    hID = oldPDA.pArr[s];
                    if (psHolder[s - 1] == s) {
                        newPDA.mLens[hID] = 0;
                    } else {
                        prvLowM_Zero = oldPDA.mLens[oldPDA.pArr[s]];
                        int seekIndex = s - 1;
                        while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) != '0') {
                            prvLowM_Zero = std::min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_Zero);
                            seekIndex--;
                        }
                    }
                    newPDA.mLens[hID] = prvLowM_Zero + 1;
                    temMinDZ = newPDA.mLens[hID];
                    for (int k = s + 1; k < e; k++) {
                        hID = oldPDA.pArr[k];
                        newPDA.mLens[hID] = newPDA.mLens[hID] + 1;
                        temMinDZ = std::min(temMinDZ, newPDA.mLens[hID]);
                    }

                    minDZ[i] = temMinDZ;
                    minDZ_Ready[i] = true;
                    return;
                }

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
                    prvLowM_Zero = oldPDA.mLens[oldPDA.pArr[s]];
                    int seekIndex = s - 1;
                    while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) != '0') {
                        prvLowM_Zero = std::min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_Zero);
                        seekIndex--;
                    }
                }

                if (minOneSearch) {
                    int seekIndex = s - 1;
                    int search_BlockID = i - 1;

                    prvLowM_One = oldPDA.mLens[oldPDA.pArr[seekIndex]];

                    while (true) {
                        if (minDZ_Ready[search_BlockID] == false) {
                            std::this_thread::sleep_for(std::chrono::milliseconds(msD_Wait)); // Simulate Thread.Sleep
                            continue;
                        }
                        if (minDZ[search_BlockID] != blockHasOneSignal) {
                            prvLowM_One = std::min(minDZ[search_BlockID] - 1, prvLowM_One);
                            search_BlockID--;
                        } else {
                            break;
                        }
                        if (search_BlockID < 0) break;
                    }
                    if (search_BlockID >= 0) {
                        seekIndex = search_BlockID * blockSize + blockSize - 1;
                        while (seekIndex >= 0 && map.Get_HapVal(oldPDA.pArr[seekIndex]) == '0') {
                            prvLowM_One = std::min(oldPDA.mLens[oldPDA.pArr[seekIndex]], prvLowM_One);
                            seekIndex--;
                        }
                    }
                }

                hID = oldPDA.pArr[s];
                if (s > 0 && map.Get_HapVal(oldPDA.pArr[s - 1]) == map.Get_HapVal(oldPDA.pArr[s])) {
                    if (map.Get_HapVal(hID) == '0') {
                        newPDA.mLens[hID] = oldPDA.mLens[hID] + 1;
                        prvLowM_One = std::numeric_limits<int>::max();
                        prvLowM_Zero = std::min(prvLowM_Zero, oldPDA.mLens[hID]);
                    } else {
                        newPDA.mLens[hID] = oldPDA.mLens[hID] + 1;
                        prvLowM_Zero = std::numeric_limits<int>::max();
                        prvLowM_One = std::numeric_limits<int>::max();
                    }
                } else {
                    if (map.Get_HapVal(hID) == '0') {
                        newPDA.mLens[hID] = std::min(prvLowM_One, oldPDA.mLens[hID]) + 1;
                        prvLowM_One = std::numeric_limits<int>::max();
                        prvLowM_Zero = std::min(prvLowM_Zero, oldPDA.mLens[hID]);
                    } else {
                        newPDA.mLens[hID] = std::min(prvLowM_Zero, oldPDA.mLens[hID]) + 1;
                        prvLowM_Zero = std::numeric_limits<int>::max();
                        prvLowM_One = std::min(prvLowM_One, oldPDA.mLens[hID]);
                    }
                }

                for (int k = s + 1; k < e; k++) {
                    hID = oldPDA.pArr[k];
                    if (map.Get_HapVal(hID) == '0') {
                        newPDA.mLens[hID] = std::min(oldPDA.mLens[hID], prvLowM_One) + 1;
                        prvLowM_One = std::numeric_limits<int>::max();
                        prvLowM_Zero = std::min(prvLowM_Zero, oldPDA.mLens[hID]);
                    } else {
                        newPDA.mLens[hID] = std::min(oldPDA.mLens[hID], prvLowM_Zero) + 1;
                        prvLowM_Zero = std::numeric_limits<int>::max();
                        prvLowM_One = std::min(prvLowM_One, oldPDA.mLens[hID]);
                    }
                }

                // (Phần còn lại của LineWork_LLM - xử lý minDZ, cập nhật globalMinDZ, etc.)
                // ...
                if (local_minDZ[i] != blockHasOneSignal) {
                    std::lock_guard<std::mutex> lock(mtx);
                    globalMinDZ = std::min(globalMinDZ.load(), local_minDZ[i]);
                }
            }, 0, ::nThread);

        // (Phần cuối của LineWork_LLM - xử lý kết quả, ghi ra output, etc.)
        // ...
    }
    void PBWT::Pal::coreP_Arr(utl::RoundTableReaderV2::HapMap map) {
        std::vector<int> local_offsets(::nThread);  // Sử dụng vector để lưu trữ offset cục bộ
    
        std::for_each(std::execution::par, std::execution::seq,
            [&](int i) {
                int s = i * blockSize;
                if (s >= ::nHap) return;
                int e = i * blockSize + blockSize;
                if (e > ::nHap) e = ::nHap;
    
                if (map.Get_HapVal(oldPDA.pArr[s]) == '0') {
                    psHolder[s] = 0;
                } else {
                    psHolder[s] = 1;
                }
    
                for (int k = s + 1; k < e; k++) {
                    if (map.Get_HapVal(oldPDA.pArr[k]) == '0') {
                        psHolder[k] = psHolder[k - 1];
                    } else {
                        psHolder[k] = psHolder[k - 1] + 1;
                    }
                }
                local_offsets[i] = psHolder[e - 1];
            }, 0, ::nThread);
    
        // Tính toán offset tuần tự
        for (size_t i = 1; i < local_offsets.size(); ++i) {
            local_offsets[i] += local_offsets[i - 1];
        }
    
        // Áp dụng offset
        std::for_each(std::execution::par, std::execution::seq,
            [&](int i) {
                int s = i * blockSize;
                if (s >= ::nHap) return;
                int e = i * blockSize + blockSize;
                if (e > ::nHap) e = ::nHap;
    
                for (int k = s; k < e; k++) {
                    if (i > 0) {
                        psHolder[k] += local_offsets[i - 1];
                    }
                }
            }, 1, ::nThread);
    
        // Sắp xếp lại pArr
        newPDA.zeroCnt = ::nHap - psHolder[::nHap - 1];
        int oneOff = newPDA.zeroCnt - 1;
    
        std::for_each(std::execution::par, std::execution::seq,
            [&](int i) {
                if (map.Get_HapVal(oldPDA.pArr[i]) == '0') {
                    newPDA.pArr[i - psHolder[i]] = oldPDA.pArr[i];
                } else {
                    newPDA.pArr[psHolder[i] + oneOff] = oldPDA.pArr[i];
                }
            }, 0, ::nHap);
    
        //  oldPDA.pArr = newPDA.pArr;  // C++ cần sao chép, không gán tham chiếu trực tiếp
        std::copy(newPDA.pArr, newPDA.pArr + ::nHap, oldPDA.pArr);
    }

    void PBWT::Pal::OneSort(utl::RoundTableReaderV2::HapMap& map, int doneIndex) {
        std::vector<int> local_offsets(::nThread);
    
        std::for_each(std::execution::par, std::execution::seq,
            [&](int i) {
                int s = i * blockSize;
                if (s >= ::nHap) return;
                int e = i * blockSize + blockSize;
                if (e > ::nHap) e = ::nHap;
    
                int zeroCnt = 0;
                for (int k = s; k < e; k++) {
                    if (map.Get_HapVal(oldPDA.pArr[k]) == '0') {
                        zeroCnt++;
                    }
                }
                local_offsets[i] = zeroCnt;
            }, 0, ::nThread);
    
        // Calculate sequential offsets
        for (size_t i = 1; i < local_offsets.size(); ++i) {
            local_offsets[i] += local_offsets[i - 1];
        }
    
        std::for_each(std::execution::par, std::execution::seq,
            [&](int i) {
                int s = i * blockSize;
                if (s >= ::nHap) return;
                int e = i * blockSize + blockSize;
                if (e > ::nHap) e = ::nHap;
    
                int zeroWrite = s;
                int oneWrite = s + local_offsets[i] - local_offsets[0]; // Corrected calculation
                if (i > 0) {
                    oneWrite -= local_offsets[0];
                }
    
                for (int k = s; k < e; k++) {
                    if (map.Get_HapVal(oldPDA.pArr[k]) == '0') {
                        newPDA.pArr[zeroWrite++] = oldPDA.pArr[k];
                        newPDA.mLens[oldPDA.pArr[k]] = oldPDA.mLens[oldPDA.pArr[k]] + 1;
                    } else {
                        newPDA.pArr[oneWrite++] = oldPDA.pArr[k];
                        newPDA.mLens[oldPDA.pArr[k]] = oldPDA.mLens[oldPDA.pArr[k]] + 1;
                    }
                }
            }, 0, ::nThread);
    
        std::swap(oldPDA, newPDA);
    }
    void PBWT::Pal::ReportLongMatches_SL(int doneIndex, utl::RoundTableReaderV2::HapMap& map) {
        for (int i = 0; i < ::nHap; i++) {
            if (oldPDA.mLens[oldPDA.pArr[i]] >= ::LLM_Len) {
                std::string sOut = std::to_string(doneIndex) + "\t" + std::to_string(oldPDA.pArr[i]) + "\t" + std::to_string(oldPDA.mLens[oldPDA.pArr[i]]);
                ::BW.Add(sOut);
            }
        }
    }
    
    void PBWT::Pal::ReportLongMatches_SL_Tail(int doneIndex) {
        for (int i = 0; i < ::nHap; i++) {
            if (oldPDA.mLens[oldPDA.pArr[i]] > 1) { // Changed from >= LLM_Len to > 1
                std::string sOut = std::to_string(doneIndex) + "\t" + std::to_string(oldPDA.pArr[i]) + "\t" + std::to_string(oldPDA.mLens[oldPDA.pArr[i]]);
                ::BW.Add(sOut);
            }
        }
    }
    void RoundTableReaderV2::LineTaker_LLM() {
        lineTag tag;
        bool dq;
        int doneIndex = 0;
    
        // First site
        while (!L2_ReadComplete || Lines.empty()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
        }
    
        tag = Lines.front();
        Lines.pop();
        HapMap map(tag, Program::nHap);
        Program::pal.InitialSort(map);
    
        doneIndex++;
    
        // Mid sites
        while (!L2_ReadComplete || !Lines.empty()) {
            while (Lines.empty()) {
                if (L2_ReadComplete) return;
                std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
            }
    
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
    
        std::cout << std::endl << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) << " L3 Line Taker Complete." << std::endl;
    
        // PBWT last site process
        Program::pal.ReportLongMatches_SL_Tail(doneIndex);
        ::BW.DoneAdding();
    }
    
    void RoundTableReaderV2::LineWork_LLM(int doneIndex, lineTag tag) {
        HapMap map(tag, Program::nHap);
    
        if (doneIndex >= ::LLM_Len) {
            Program::pal.ReportLongMatches_SL(doneIndex, map);
        }
    
        Program::pal.OneSort(map, doneIndex);
    }
    // You will need to implement the remaining methods (e.g., LineWork_LLM, etc.)
    // following the same pattern.  Remember to handle parallelism and data access carefully.
}