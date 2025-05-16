#ifndef PBWT_H
#define PBWT_H

#include <vector>
#include <queue>
#include "utl.cpp"

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

        Pal();
        void coreP_Arr(utl::RoundTableReaderV2::HapMap& map);
        void coreD_Arr(utl::RoundTableReaderV2::HapMap& map);
        void coreD_Arr_LowMAF(utl::RoundTableReaderV2::HapMap& map, int doneIndex);
        void coreD_Arr_BU_a(utl::RoundTableReaderV2::HapMap& map);
        void ReportLongMatches_SL(int siteIndex, utl::RoundTableReaderV2::HapMap& map);
        void ReportLongMatches_SL_Tail(int siteIndex);
        void ReportSetMaxMatch(int siteIndex, utl::RoundTableReaderV2::HapMap& map);
        void ReportSetMaxMatch_Tail(int siteIndex);
        void OneSort(utl::RoundTableReaderV2::HapMap& map, int doneIndex);
        void InitialSort(utl::RoundTableReaderV2::HapMap& map);
    };
};

} // namespace HP_PBWT

#endif // PBWT_H 