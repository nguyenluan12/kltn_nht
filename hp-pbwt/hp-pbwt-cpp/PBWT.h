#ifndef HP_PBWT_PBWT_H
#define HP_PBWT_PBWT_H

#include <vector>
#include <iostream>
#include <thread>
#include <future>
#include <mutex>
#include <queue>
#include <chrono>

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

            static int msD_Wait;

            static int blockHasOneSignal;

        public:
            Pal();
            // Các phương thức của Pal sẽ được thêm sau
        };

        // Các phương thức của PBWT (nếu có)
    };

} // namespace HP_PBWT

#endif // HP_PBWT_PBWT_H