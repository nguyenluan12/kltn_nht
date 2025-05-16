#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <fstream>
#include <iterator>
#include <bits/stdc++.h>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace std;
using namespace chrono;
using namespace boost::iostreams;

void build_prefix_and_divergence_arrays(const vector<vector<int>>& X, vector<int>& ppa, vector<int>& div) {
    // Number of haplotypes
    int M = X.size();

    // Number of variable sites
    int N = X[0].size();

    // Initialize positional prefix array and divergence array
    ppa.resize(M);
    div.resize(M, 0);

    // Fill initial ppa with indices
    for (int i = 0; i < M; ++i) {
        ppa[i] = i;
    }

    // Iterate over variants
    for (int k = 0; k < N - 1; ++k) {
        // Temporary vectors to store intermediates
        vector<int> a, b, d, e;
        int p = k + 1, q = k + 1;

        // Iterate over haplotypes in reverse prefix sorted order
        for (int i = 0; i < M; ++i) {
            int index = ppa[i];
            int match_start = div[i];
            const vector<int>& haplotype = X[index];
            int allele = haplotype[k];

            // Update intermediates
            if (match_start > p) {
                p = match_start;
            }
            if (match_start > q) {
                q = match_start;
            }

            if (allele == 0) {
                a.push_back(index);
                d.push_back(p);
                p = 0;
            } else {
                b.push_back(index);
                e.push_back(q);
                q = 0;
            }
        }

        // Construct the new arrays for k+1 by concatenating intermediates
        ppa.clear();
        ppa.reserve(a.size() + b.size());
        ppa.insert(ppa.end(), a.begin(), a.end());
        ppa.insert(ppa.end(), b.begin(), b.end());

        div.clear();
        div.reserve(d.size() + e.size());
        div.insert(div.end(), d.begin(), d.end());
        div.insert(div.end(), e.begin(), e.end());
    }
}

vector<vector<int>> report_long_matches(const vector<vector<int>>& X, int L) {
    int M = X.size();       // Number of haplotypes
    int N = X[0].size();    // Number of variable sites

    // Initialize positional prefix array and divergence array
    vector<int> ppa(M);
    iota(ppa.begin(), ppa.end(), 0);  // Fill ppa with 0, 1, 2, ..., M-1
    vector<int> div(M, 0);

    vector<vector<int>> result; // To store results

    // Iterate over variants
    for (int k = 0; k < N; ++k) {
        vector<int> a, b, d, e;
        int p = k + 1;
        int q = k + 1;
        vector<int> ma, mb;

        // Iterate over haplotypes in reverse prefix sorted order
        for (size_t i = 0; i < ppa.size(); ++i) {
            int index = ppa[i];
            int match_start = div[i];

            // Report matches
            if (match_start > k - L) {
                if (!ma.empty() && !mb.empty()) {
                    // Store position k, ma, and mb in the result vector
                    vector<int> result_entry;
                    result_entry.push_back(k);
                    result_entry.insert(result_entry.end(), ma.begin(), ma.end());
                    result_entry.push_back(-1); // Separator for mb
                    result_entry.insert(result_entry.end(), mb.begin(), mb.end());
                    result.push_back(result_entry);
                }
                ma.clear();
                mb.clear();
            }

            // Current haplotype
            const vector<int>& haplotype = X[index];
            int allele = haplotype[k];

            // Update intermediates
            if (match_start > p) p = match_start;
            if (match_start > q) q = match_start;

            // Update intermediates
            if (allele == 0) {
                a.push_back(index);
                d.push_back(p);
                p = 0;
                ma.push_back(index);
            } else {
                b.push_back(index);
                e.push_back(q);
                q = 0;
                mb.push_back(index);
            }
        }

        // Report any remaining matches including final haplotype
        if (!ma.empty() && !mb.empty()) {
            vector<int> result_entry;
            result_entry.push_back(k);
            result_entry.insert(result_entry.end(), ma.begin(), ma.end());
            result_entry.push_back(-1); // Separator for mb
            result_entry.insert(result_entry.end(), mb.begin(), mb.end());
            result.push_back(result_entry);
        }

        // Construct the new arrays for k+1
        if (k < N - 1) {
            vector<int> new_ppa;
            vector<int> new_div;

            new_ppa.reserve(a.size() + b.size());
            new_div.reserve(d.size() + e.size());

            new_ppa.insert(new_ppa.end(), a.begin(), a.end());
            new_ppa.insert(new_ppa.end(), b.begin(), b.end());

            new_div.insert(new_div.end(), d.begin(), d.end());
            new_div.insert(new_div.end(), e.begin(), e.end());

            ppa = move(new_ppa);
            div = move(new_div);
        }
    }

    return result;
}

vector<vector<int> > read_hap(const string &filename) {
    vector<vector<int> > res;
    ifstream ifs(filename);       // open the file
    string tempstr;
    while (getline(ifs, tempstr)) {
        stringstream lineStream(tempstr);
        vector<int> numbers(istream_iterator<int>(lineStream),{});
        if (res.empty()) {
            for (auto number: numbers) {
                res.push_back(vector<int>{number});
            }
        }
        else {
            for (int i = 0; i < numbers.size(); ++i) {
                res[i].push_back(numbers[i]);
            }
        }
    }
    ifs.close();
    return res;
}

vector<vector<int> > read_hap_gz(const string &filename){
    vector<vector<int> > res;
    ifstream file(filename, ios_base::in | ios_base::binary);
    filtering_streambuf<input> inbuf;
    inbuf.push(gzip_decompressor());
    inbuf.push(file);
    //Convert streambuf to istream
    istream instream(&inbuf);
    //Iterate lines
    string line;
    while(getline(instream, line)) {
        stringstream lineStream(line);
        vector<int> numbers(istream_iterator<int>(lineStream), {});
        if (res.empty()) {
            for (auto number: numbers) {
                res.push_back(vector<int>{number});
            }
        } else {
            for (int i = 0; i < numbers.size(); ++i) {
                res[i].push_back(numbers[i]);
            }
        }
    }
    file.close();
    return res;
}
vector<vector<int>> read_vcf_gz_to_hapmatrix(const string &filename) {
    ifstream file(filename, ios_base::in | ios_base::binary);
    filtering_streambuf<input> inbuf;
    inbuf.push(gzip_decompressor());
    inbuf.push(file);
    istream instream(&inbuf);

    string line;
    vector<vector<int>> haplotypes;
    bool header_found = false;
    int num_samples = 0;

    while (getline(instream, line)) {
        if (line.rfind("##", 0) == 0) continue;  // Skip meta-information lines
        if (line.rfind("#CHROM", 0) == 0) {
            header_found = true;
            istringstream ss(line);
            vector<string> tokens{istream_iterator<string>{ss}, {}};
            num_samples = tokens.size() - 9;
            haplotypes.resize(num_samples * 2); // Two haplotypes per sample
            continue;
        }

        if (!header_found) continue;

        istringstream ss(line);
        vector<string> tokens{istream_iterator<string>{ss}, {}};

        for (int i = 9; i < tokens.size(); ++i) {
            const string &gt = tokens[i];
            if (gt.size() < 3) continue;

            char a = gt[0];
            char b = gt[2];
            haplotypes[2 * (i - 9)].push_back(a == '1' ? 1 : 0);
            haplotypes[2 * (i - 9) + 1].push_back(b == '1' ? 1 : 0);
        }
    }

    file.close();
    return haplotypes;
}
vector<vector<int>> read_vcf_to_hapmatrix(const string &filename) {
    ifstream infile(filename);
    string line;
    vector<vector<int>> haplotypes;
    bool header_found = false;
    int num_samples = 0;

    while (getline(infile, line)) {
        if (line.rfind("##", 0) == 0) continue;  // Skip metadata

        if (line.rfind("#CHROM", 0) == 0) {
            header_found = true;
            istringstream ss(line);
            vector<string> tokens{istream_iterator<string>{ss}, {}};
            num_samples = tokens.size() - 9;
            haplotypes.resize(num_samples * 2);
            continue;
        }

        if (!header_found) continue;

        istringstream ss(line);
        vector<string> tokens{istream_iterator<string>{ss}, {}};

        for (int i = 9; i < tokens.size(); ++i) {
            string gt = tokens[i];
            if (gt.size() < 3) continue;

            char a = gt[0];
            char b = gt[2];
            haplotypes[2 * (i - 9)].push_back(a == '1' ? 1 : 0);
            haplotypes[2 * (i - 9) + 1].push_back(b == '1' ? 1 : 0);
        }
    }

    infile.close();
    return haplotypes;
}

int main(int argc, char* argv[]) {
    // vector<vector<int> > X = {
    //     {0, 1, 0, 1, 0, 1},
    //     {1, 1, 0, 0, 0, 1},
    //     {1, 1, 1, 1, 1, 1},
    //     {0, 1, 1, 1, 1, 0},
    //     {0, 0, 0, 0, 0, 0},
    //     {1, 0, 0, 0, 1, 0},
    //     {1, 1, 0, 0, 0, 1},
    //     {0, 1, 0, 1, 1, 0}
    // };

    // vector<vector<int> > X = {
    //     {0, 1, 1, 1, 1, 0, 0, 0, 0, 0},
    //     {1, 1, 0, 1, 0, 0, 0, 1, 0, 1},
    //     {0, 1, 0, 1, 0, 0, 0, 1, 0, 1},
    //     {0, 0, 1, 1, 0, 1, 0, 1, 0, 0},
    //     {1, 0, 1, 1, 0, 0, 0, 0, 1, 0},
    //     {1, 0, 1, 1, 0, 0, 0, 0, 1, 1},
    //     {0, 1, 1, 1, 1, 0, 0, 0, 0, 1},
    //     {0, 1, 0, 0, 1, 0, 1, 0, 0, 1},
    //     {0, 1, 0, 0, 0, 0, 0, 0, 1, 0},
    //     {1, 0, 1, 1, 0, 0, 0, 0, 1, 0}
    // };

    // vector<vector<int> > X;
    // int width = 100000;
    // int height = 200;
    // for (int i = 0; i < height; ++i) {
    //     vector<int> tmp;
    //     for (int j = 0; j < width; ++j) {
    //         tmp.push_back(rand() % 2);
    //     } {
    //         X.push_back(tmp);
    //         printf("Init: %lu/%d\n", X.size(), height);
    //     }
    // }
    // printf("Init: done\n");

    // string f = "/home/pbthang/CLionProjects/untitled/chr6_5m.hap";
    // printf("Read file: %s\n",f.c_str());
    // vector<vector<int> > X = read_hap(f);

    //  string f = "/Users/NguyenLuan/pbwt-comparison/random_data/random_100x100000.vcf";
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <vcf_file_path>" << endl;
        return 1;
    }
    string f = argv[1];
     printf("Read file: %s\n",f.c_str());
    //  vector<vector<int> > X = read_hap(f);
    //  vector<vector<int> > X = read_hap_gz(f);
    // vector<vector<int> > X = read_vcf_gz_to_hapmatrix(f);
    vector<vector<int> > X = read_vcf_to_hapmatrix(f);



    const int retry = 50;
    const int L = 4;
    printf("Retry: %d\n",retry);
    printf("Matched length: %d\n",L);

        vector<int> ppa;
    vector<int> div;
    signed long int duration_total_build = 0;
        high_resolution_clock::time_point start,stop;
        vector<vector<vector<int> > > res;
        for (int j = 0; j < retry; ++j) {
            start = high_resolution_clock::now();
            build_prefix_and_divergence_arrays(X, ppa, div);
            stop = high_resolution_clock::now();
            duration_total_build += duration_cast<microseconds>(stop - start).count();
        }
        printf("1 thread(s), build time: %ld us\n", duration_total_build/retry);

        signed long int duration_total_match = 0;
        for (int j = 0; j < retry; ++j) {
            start = high_resolution_clock::now();
            vector<vector<int>> matches = report_long_matches(X, L);
            stop = high_resolution_clock::now();
            duration_total_match += duration_cast<microseconds>(stop - start).count();

        }
        printf("1 thread(s), match time: %ld us\n",  duration_total_match/retry);
        printf("1 thread(s), total time: %ld us\n\n",  duration_total_match/retry+duration_total_build/retry);

}
