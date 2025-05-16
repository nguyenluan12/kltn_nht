#include "utl.h"
#include <sstream>
#include <iomanip>
#include "Program.h" // Include Program.h to access global variables

namespace HP_PBWT {

    std::string utl::IntTo_BMK_Str(int n) {
        if (n >= 1000000000) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(0) << (double)n / 1000000000.0 << "b";
            return oss.str();
        } else if (n >= 1000000) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(0) << (double)n / 1000000.0 << "m";
            return oss.str();
        } else if (n >= 1000) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(0) << (double)n / 1000.0 << "k";
            return oss.str();
        } else {
            return std::to_string(n);
        }
    }

    utl::RoundTableReaderV2::HapMap::HapMap(utl::RoundTableReaderV2::lineTag tag, int nHap) : tag(tag) {
        if (tag.dataCross == false) {
            LastHID_FirstBlk = nHap - 1;
            return;
        }

        LastHID_FirstBlk = (::rtrBlockSize - tag.charIndex - 1) / 2;
        charOff_SecBlk = (::rtrBlockSize - tag.charIndex) % 2;
        secBlkID = (tag.blockID + 1) % ::rtrN_Block;
    }

    char utl::RoundTableReaderV2::HapMap::Get_HapVal(int hID) {
        if (tag.dataCross == false || hID <= LastHID_FirstBlk) {
            return ::BufferArr[tag.blockID][tag.charIndex + hID * 2];
        }

        return ::BufferArr[secBlkID][charOff_SecBlk + (hID - LastHID_FirstBlk - 1) * 2];
    }

    utl::RoundTableReaderV2::lineTag::lineTag(int blockID, int charIndex, bool dataCross, bool headerCross)
        : blockID(blockID), charIndex(charIndex), dataCross(dataCross), headerCross(headerCross) {}

    utl::RoundTableReaderV2::RoundTableReaderV2(std::string path) : VCF_Path(path) {
        std::cout << std::chrono::system_clock::now().time_since_epoch().count() << " Initializing RTR, " << ::rtrBlockSize << " by " << ::rtrN_Block << " ..." << std::endl;

        L1_nRead.resize(L1_nBlock);
        L1_Stat.resize(L1_nBlock);
        L1_Seq.resize(L1_nBlock);
        BufferArr.resize(L1_nBlock);

        for (int i = 0; i < L1_nBlock; i++) {
            L1_Stat[i] = -1;
            L1_nRead[i] = -1;
            L1_Seq[i] = -1;
            BufferArr[i].resize(L1_blockSize);
        }
        std::cout << "Done." << std::endl;
    }

    void utl::RoundTableReaderV2::BlockReader() {
        std::ifstream sr(VCF_Path);
        int nBlockRead = 0;

        int len = -1;
        while (len != 0) {
            while (L1_Stat[L1_Add_Index] != -1) {
                std::this_thread::sleep_for(std::chrono::milliseconds(msAddWait));
                continue;
            }

            len = 0;
            for (int i = 0; i < L1_blockSize && sr.get(BufferArr[L1_Add_Index][i]); ++i) {
                len++;
            }

            L1_nRead[L1_Add_Index] = len;
            L1_Stat[L1_Add_Index] = 0;

            L1_Seq[L1_Add_Index] = nBlockRead;
            nBlockRead++;

            L1_Add_Index++;
            L1_Add_Index = L1_Add_Index % L1_nBlock;
        }

        sr.close();
        L1_ReadComplete = true;
        std::cout << std::endl;
        std::cout << std::chrono::system_clock::now().time_since_epoch().count() << " Block Reader Complete." << std::endl;
    }

    void utl::RoundTableReaderV2::ProbTop() {
        while (L1_Stat[0] != 0) {
            std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
            continue;
        }

        while (BufferArr[0][L1_rPTR] == '#' && BufferArr[0][L1_rPTR + 1] == '#') {
            while (BufferArr[0][L1_rPTR] != '\n') {
                L1_rPTR++;
            }
            L1_rPTR++;
        }

        ::nIndv = 0;
        while (BufferArr[0][L1_rPTR] != '\n') {
            if (BufferArr[0][L1_rPTR] == '\t') {
                ::nIndv++;
            }
            L1_rPTR++;
        }
        ::nIndv = ::nIndv - 8;
        ::nHap = ::nIndv * 2;
        dataLineLen = ::nIndv * 4;
        std::cout << "nIndv. " << ::nIndv << " dataLineLen " << dataLineLen << std::endl;
        L1_rPTR++;
    }

    void utl::RoundTableReaderV2::LineParser() {
        bool headerCross = false;

        int tabCnt = 0;
        int nextBlock;
        int i;

        while (!L1_ReadComplete || L1_Stat[L1_Taker_Index] == 0) {

            if (L1_Stat[L1_Taker_Index] != 0 || Lines.size() > 1000) {
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

                            std::lock_guard<std::mutex> lock(lines_mutex);
                            Lines.push(lineTag(L1_Taker_Index, i, true, false));

                            break;

                        } else {

                            if (headerCross) {
                                std::lock_guard<std::mutex> lock(lines_mutex);
                                Lines.push(lineTag(L1_Taker_Index, i, false, true));

                                headerCross = false;
                            } else {
                                std::lock_guard<std::mutex> lock(lines_mutex);
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

                    std::cout << std::endl;
                    std::cout << std::chrono::system_clock::now().time_since_epoch().count() << " L2 Read Complete." << std::endl;

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

        std::cout << std::endl;
        std::cout << std::chrono::system_clock::now().time_since_epoch().count() << " L2 Read Complete." << std::endl;
    }

    void utl::RoundTableReaderV2::LineTaker_LLM() {
        lineTag tag;
        bool dq;

        int doneIndex = 0;
        while (!L2_ReadComplete || !Lines.empty()) { // C++: !Lines.empty() instead of Lines.Count > 0

            {
                std::lock_guard<std::mutex> lock(lines_mutex);
                dq = Lines.empty();
                if (!dq) {
                    tag = Lines.front();
                    Lines.pop();
                }
            }

            if (dq) {
                std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
                continue;
            }

            LineWork_LLM(doneIndex, tag);
            doneIndex++;
        }

        std::cout << std::endl;
        std::cout << std::chrono::system_clock::now().time_since_epoch().count() << " L3 Read Complete." << std::endl;
    }

    void utl::RoundTableReaderV2::LineTaker_SMM() {
        lineTag tag;
        bool dq;

        int lineCnt = 0;
        while (!L2_ReadComplete || !Lines.empty()) {
            {
                std::lock_guard<std::mutex> lock(lines_mutex);
                dq = Lines.empty();
                if (!dq) {
                    tag = Lines.front();
                    Lines.pop();
                }
            }

            if (dq) {
                std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
                continue;
            }

            LineWork_SMM(lineCnt, tag);
            lineCnt++;
        }

        std::cout << std::endl;
        std::cout << std::chrono::system_clock::now().time_since_epoch().count() << " L3 Read Complete." << std::endl;
    }

    void utl::RoundTableReaderV2::LineTaker_DBG() {
        lineTag tag;
        bool dq;

        while (!L2_ReadComplete || !Lines.empty()) {
            {
                std::lock_guard<std::mutex> lock(lines_mutex);
                dq = Lines.empty();
                if (!dq) {
                    tag = Lines.front();
                    Lines.pop();
                }
            }

            if (dq) {
                std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
                continue;
            }

            LineWork_DBG(tag);
        }

        std::cout << std::endl;
        std::cout << std::chrono::system_clock::now().time_since_epoch().count() << " L3 Read Complete." << std::endl;
    }

    void utl::RoundTableReaderV2::LineWork_LLM(int doneIndex, lineTag tag) {
        // Implementation for LineWork_LLM
    }

    void utl::RoundTableReaderV2::LineWork_SMM(int lineCnt, lineTag tag) {
        // Implementation for LineWork_SMM
    }

    void utl::RoundTableReaderV2::LineWork_DBG(lineTag tag) {
        // Implementation for LineWork_DBG
    }

    utl::BufferWriter::BufferWriter(std::string outPath) : outPath(outPath) {
        sw.open(outPath);
        if (!sw.is_open()) {
            std::cerr << "Error opening output file: " << outPath << std::endl;
            // Handle error (e.g., throw exception, exit program)
        }
    }

    void utl::BufferWriter::Add(std::string s) {
        {
            std::lock_guard<std::mutex> lock(buffer_mutex);
            buffer.push(s);
        }
        buffer_cv.notify_one(); // Notify the writer thread that data is available
    }

    void utl::BufferWriter::DoneAdding() {
        {
            std::lock_guard<std::mutex> lock(buffer_mutex);
            AddComplete = true;
        }
        buffer_cv.notify_one(); // Notify the writer thread to exit
    }

    void utl::BufferWriter::Run() {
        std::string sOut;
        while (!AddComplete || !buffer.empty()) {
            {
                std::unique_lock<std::mutex> lock(buffer_mutex);
                buffer_cv.wait(lock, [this] { return AddComplete || !buffer.empty(); }); // Wait with condition
                if (!buffer.empty()) {
                    sOut = buffer.front();
                    buffer.pop();
                } else {
                    continue;
                }
            }
            if (!sOut.empty()) {
                sw << sOut << "\n";
            }
        }
        sw.close();
        std::cout << std::chrono::system_clock::now().time_since_epoch().count() << " Buffer Writer Complete." << std::endl;
        // Exit is generally not recommended from a thread. Consider better thread management.
        // std::exit(0);
    }

    std::string IntTo_BMK_Str(int n) {
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

    RoundTableReaderV2::lineTag::lineTag(int blockID, int charIndex, bool dataCross, bool headerCross)
        : blockID(blockID), charIndex(charIndex), dataCross(dataCross), headerCross(headerCross) {}

    RoundTableReaderV2::HapMap::HapMap(lineTag tag, int nHap) : tag(tag) {
        if (!tag.dataCross) {
            LastHID_FirstBlk = nHap - 1;
            return;
        }
        LastHID_FirstBlk = (Program::rtrBlockSize - tag.charIndex - 1) / 2;
        charOff_SecBlk = (Program::rtrBlockSize - tag.charIndex) % 2;
        secBlkID = (tag.blockID + 1) % Program::rtrN_Block;
    }

    char RoundTableReaderV2::HapMap::Get_HapVal(int hID) {
        if (!tag.dataCross || hID <= LastHID_FirstBlk) {
            return Program::BufferArr[tag.blockID][tag.charIndex + hID * 2];
        }
        return Program::BufferArr[secBlkID][charOff_SecBlk + (hID - LastHID_FirstBlk - 1) * 2];
    }

    RoundTableReaderV2::RoundTableReaderV2(const std::string& path) : VCF_Path(path),
        L1_blockSize(Program::rtrBlockSize), L1_nBlock(20), L2_nBlock(1000),
        L1_nRead(L1_nBlock, 0),
        L1_Stat(L1_nBlock, -1),
        L1_Seq(L1_nBlock, -1),
        BufferArr(L1_nBlock),
        L1_rPTR(0),
        L1_Add_Index(0),
        L1_Taker_Index(0),
        L1_ReadComplete(false),
        L2_ReadComplete(false),
        dataLineLen(0) {
        for (int i = 0; i < L1_nBlock; ++i) {
            BufferArr[i].resize(L1_blockSize);
        }
        std::cout << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) << " Initializing RTR, " << L1_blockSize << " by " << L1_nBlock << " ..." << std::endl;
    }

    void RoundTableReaderV2::BlockReader() {
        std::ifstream sr(VCF_Path, std::ios::binary); // Open in binary mode for efficiency
        if (!sr.is_open()) {
            throw std::runtime_error("Could not open file: " + VCF_Path);
        }
        int nBlockRead = 0;
        size_t len = 1; // Initialize to non-zero to enter the loop

        while (len > 0) {
            while (L1_Stat[L1_Add_Index] != -1) {
                std::this_thread::sleep_for(std::chrono::milliseconds(msAddWait));
            }
            sr.read(BufferArr[L1_Add_Index].data(), L1_blockSize);
            len = sr.gcount();
            L1_nRead[L1_Add_Index] = len;
            L1_Stat[L1_Add_Index] = 0;
            L1_Seq[L1_Add_Index] = nBlockRead++;
            L1_Add_Index = (L1_Add_Index + 1) % L1_nBlock;
        }

        sr.close();
        L1_ReadComplete = true;
        std::cout << std::endl << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) << " Block Reader Complete." << std::endl;
    }

    void RoundTableReaderV2::ProbTop() {
        while (L1_Stat[0] != 0) {
            std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
        }

        // Skip header lines
        while (L1_rPTR < L1_nRead[0] && BufferArr[0][L1_rPTR] == '#' && BufferArr[0][L1_rPTR + 1] == '#') {
            while (L1_rPTR < L1_nRead[0] && BufferArr[0][L1_rPTR] != '\n') {
                L1_rPTR++;
            }
            if (L1_rPTR < L1_nRead[0]) L1_rPTR++;
        }

        // Count individuals
        Program::nIndv = 0;
        while (L1_rPTR < L1_nRead[0] && BufferArr[0][L1_rPTR] != '\n') {
            if (BufferArr[0][L1_rPTR] == '\t') {
                Program::nIndv++;
            }
            L1_rPTR++;
        }
        Program::nIndv = std::max(0, Program::nIndv - 8);
        Program::nHap = Program::nIndv * 2;
        dataLineLen = Program::nIndv * 4;
        std::cout << "nIndv. " << Program::nIndv << " dataLineLen " << dataLineLen << std::endl;
        if (L1_rPTR < L1_nRead[0]) L1_rPTR++;
    }

    void RoundTableReaderV2::LineParser() {
        bool headerCross = false;
        int tabCnt = 0;
        int nextBlock;
        size_t i;

        while (!L1_ReadComplete || L1_Stat[L1_Taker_Index] == 0) {
            if (L1_Stat[L1_Taker_Index] != 0 || Lines.size() > L2_nBlock) {
                std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
                continue;
            }

            nextBlock = (L1_Taker_Index + 1) % L1_nBlock;

            for (i = L1_rPTR; i < L1_nRead[L1_Taker_Index]; ++i) {
                if (BufferArr[L1_Taker_Index][i] == '\t') {
                    tabCnt++;
                    if (tabCnt == 9) {
                        tabCnt = 0;
                        i++;

                        // Cross block or end block
                        if (i + dataLineLen >= L1_nRead[L1_Taker_Index]) {
                            L1_rPTR = (i + dataLineLen) % L1_nRead[L1_Taker_Index];
                            while (L1_Stat[nextBlock] != 0) {
                                std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
                            }
                            Lines.push(lineTag(L1_Taker_Index, i, true, false));
                            break;
                        } else {
                            // Within block
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

            // Move to next block
            while (L1_Stat[nextBlock] != 0) {
                if (L1_ReadComplete) {
                    L2_ReadComplete = true;
                    std::cout << std::endl << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) << " L2 Read Complete." << std::endl;
                    return;
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
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
        std::cout << std::endl << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) << " L2 Read Complete." << std::endl;
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

        // Mid sites
        while (!L2_ReadComplete || !Lines.empty()) {
            while (Lines.empty()) {
                if (L2_ReadComplete) return;  // Ensure exit if reading is done
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
    }
    BufferWriter::BufferWriter(const std::string& outPath) : outPath(outPath), addComplete(false) {}

void BufferWriter::Add(const std::string& s) {
    {
        std::lock_guard<std::mutex> lock(mtx);
        buffer.push(s);
    }
    cv.notify_one(); // Notify the writer thread that data is available
}

void BufferWriter::DoneAdding() {
    {
        std::lock_guard<std::mutex> lock(mtx);
        addComplete = true;
    }
    cv.notify_one(); // Ensure the writer thread wakes up to check for completion
}

void BufferWriter::Run() {
    sw.open(outPath);
    if (!sw.is_open()) {
        std::cerr << "Error opening output file: " << outPath << std::endl;
        return; // Or throw an exception
    }
    sw.exceptions(std::ofstream::failbit | std::ofstream::badbit); // Enable exceptions for file errors

    std::string sOut;
    while (!addComplete || !buffer.empty()) {
        {
            std::unique_lock<std::mutex> lock(mtx);
            cv.wait(lock, [this]{ return addComplete || !buffer.empty(); }); // Wait with condition
            if (!buffer.empty()) {
                sOut = buffer.front();
                buffer.pop();
            } else {
                continue; // Spurious wakeup, check condition again
            }
        }
        try {
            sw << sOut << "\n";
        } catch (const std::ofstream::failure& e) {
            std::cerr << "Exception writing to file: " << e.what() << std::endl;
            break;
        }
    }
    sw.close();
    std::cout << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) << " Buffer Writer Complete." << std::endl;
}

void RoundTableReaderV2::LineTaker_SMM() {
    lineTag tag;
    bool dq;
    int siteCnt = 1;

    // First site
    while (!L2_ReadComplete || Lines.empty()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
    }

    tag = Lines.front();
    Lines.pop();
    HapMap map(tag, Program::nHap);
    Program::pal.InitialSort(map);

    siteCnt++;

    // Mid sites
    while (!L2_ReadComplete || !Lines.empty()) {
        while (Lines.empty()) {
            if (L2_ReadComplete) return;
            std::this_thread::sleep_for(std::chrono::milliseconds(msTakeWait));
        }

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

    std::cout << std::endl << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) << " L3 Line Taker Complete." << std::endl;

    // PBWT last site process
    Program::pal.ReportSetMaxMatch_Tail(siteCnt);
    ::BW.DoneAdding();
}

void RoundTableReaderV2::LineWork_SMM(int lineCnt, lineTag tag) {
    HapMap map(tag, Program::nHap);

    if (lineCnt >= ::LLM_Len) {
        Program::pal.ReportSetMaxMatch(lineCnt, map);
    }

    Program::pal.OneSort(map, 0); // Note: doneIndex is 0 here, might need adjustment
}
} // namespace HP_PBWT