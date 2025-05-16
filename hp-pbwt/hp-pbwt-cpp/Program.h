#ifndef HP_PBWT_PBWT_H
#define HP_PBWT_PBWT_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <execution> // For parallel algorithms (if available)
#include "utl.h"
#include "Program.h"

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

                PDA(int nHapTotal);
            };

            class UpLowRange {
            public:
                int UpRange;
                int LowerRange;

                UpLowRange(int up, int low);
            };

        private:
            static PDA newPDA;
            static PDA oldPDA;
            static PDA temPDA;
            static std::vector<int> psHolder;
            static std::vector<int> minDZ;
            static std::vector<bool> minDZ_Ready;
            static std::vector<int> offsetsHolder;
            static int blockSize;
            //static std::execution::parallel_policy pOp; // C++17 Parallel Algorithms
            static int msD_Wait;
            static int blockHasOneSignal;

        public:
            Pal();

            void coreP_Arr(utl::RoundTableReaderV2::HapMap& map);
            void coreD_Arr(utl::RoundTableReaderV2::HapMap& map);
            void ResetMinDZ();
            void coreD_Arr_LowMAF(utl::RoundTableReaderV2::HapMap& map, int doneIndex);
            void coreD_Arr_BU_a(utl::RoundTableReaderV2::HapMap& map);
        };
    };

} // namespace HP_PBWT

#endif // HP_PBWT_PBWT_H