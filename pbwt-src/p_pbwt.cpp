#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <omp.h>
#include "pbwt.h"
#include <thread>
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

vector<vector<vector<int> > > build_prefix_and_divergence_arrays(const vector<vector<int> > &X, int thread) {
    // Number of haplotypes
    int M = static_cast<int>(X.size());

    // Number of variable sites
    int N = static_cast<int>(X[0].size());
    vector<vector<int> > ppa_t;
    vector<vector<int> > div_t;
    vector<vector<int> > range;
    omp_set_dynamic(0);
    omp_set_num_threads(thread);
#pragma omp parallel
    {
        vector<int> ppa;
        vector<int> div;
        // Initialize positional prefix array and divergence array
        ppa.resize(M);
        div.resize(M, 0);

        // Fill initial ppa with indices
        for (int i = 0; i < M; ++i) {
            ppa[i] = i;
        }
        const int nthreads = omp_get_num_threads();
        const int ithread = omp_get_thread_num();
        int start = ithread * N / nthreads;
        int finish = (ithread + 1) * N / nthreads;
        // Iterate over variants
        for (int k = start; k < finish; ++k) {
            // Temporary vectors to store intermediates
            vector<int> a, b, d, e;
            int p = k + 1, q = k + 1;

            // Iterate over haplotypes in reverse prefix sorted order
            for (int i = 0; i < M; ++i) {
                int index = ppa[i];
                int match_start = div[i];
                const vector<int> &haplotype = X[index];
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
                    p = start;
                } else {
                    b.push_back(index);
                    e.push_back(q);
                    q = start;
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
#pragma omp for schedule(static) ordered

        for (int i = 0; i < omp_get_num_threads(); ++i) {
#pragma omp ordered
            {
                ppa_t.push_back(ppa);
                div_t.push_back(div);
                range.push_back(vector<int>{start, finish});
            }
        }
    }


    for (int slice = 1; slice <= range.size() - 1; ++slice) {
        // printf("\nStep %d:\n", slice);
        algorithm4(M, range[slice][1], range[slice][0], ppa_t[slice], div_t[slice]
                   , ppa_t[slice - 1], div_t[slice - 1]);
    }
    return vector<vector<vector<int> > >{ppa_t, div_t, range};
}


void algorithm4(const int M, const int k, const int b,
                vector<int> &sorted_k, vector<int> &start_k,
                const vector<int> &sorted_k_b, vector<int> &start_k_b
) {
    int index_k_b[M];
    int group_id = 0;
    for (int i = 0; i < M; ++i) {
        index_k_b[sorted_k_b[i]] = i;
    }
    for (int i = 0; i < M; ++i) {
        if (start_k[i] != b) {
            if (i - group_id > 1) {
                algorithm3(group_id, i, index_k_b, sorted_k, start_k, sorted_k_b, start_k_b);
            }
            group_id = i;
        }
    }
    if (M - group_id > 1) {
        algorithm3(group_id, M, index_k_b, sorted_k, start_k, sorted_k_b, start_k_b);
    }
}

void algorithm3(const int id_start, const int id_end, const int index_k_b[],
                vector<int> &sorted_k, vector<int> &start_k,
                const vector<int> &sorted_k_b, vector<int> &start_k_b
) {
    int arr[id_end - id_start];
    for (int i = id_start; i < id_end; ++i) {
        arr[i - id_start] = index_k_b[sorted_k[i]];
    }
    sort(arr, arr + (id_end - id_start));
    for (int i = id_start; i < id_end; ++i) {
        sorted_k[i] = sorted_k_b[arr[i - id_start]];
    }
    for (int i = id_start + 1; i < id_end; ++i) {
        int scan_start = arr[i - id_start - 1] + 1;
        int scan_stop = arr[i - id_start];
        auto new_val = max_element(start_k_b.begin() + scan_start, start_k_b.begin() + scan_stop + 1);
        start_k[i] = *new_val;
    }
}

vector<vector<int> > report_long_matches(const vector<vector<int> > &X, const int L, vector<vector<vector<int> > > res,
                                         int thread) {
    int M = static_cast<int>(X.size()); // Number of haplotypes
    int N = static_cast<int>(X[0].size()); // Number of variable sites

    vector<vector<int> > ppa_t = res[0];
    vector<vector<int> > div_t = res[1];
    vector<vector<int> > range = res[2];
    // Initialize positional prefix array and divergence array
    vector<int> ppa(M);
    iota(ppa.begin(), ppa.end(), 0); // Fill ppa with 0, 1, 2, ..., M-1
    vector<int> div(M, 0);

    ppa_t.insert(ppa_t.begin(), ppa);
    div_t.insert(div_t.begin(), div);

    vector<vector<int> > results; // To store results
    {
        omp_set_dynamic(0);
        omp_set_num_threads(thread);
#pragma omp parallel for private(ppa,div)
        for (int slice = 0; slice < range.size(); ++slice) {
            // Iterate over variants
            ppa = ppa_t[slice];
            div = div_t[slice];
            for (int k = range[slice][0]; k < range[slice][1]; ++k) {
                vector<int> a, b, d, e, ma, mb;
                int p = k + 1;
                int q = k + 1;

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
#pragma omp critical
                            results.push_back(result_entry);
                        }
                        ma.clear();
                        mb.clear();
                    }

                    // Current haplotype
                    const vector<int> &haplotype = X[index];
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
#pragma omp critical
                    results.push_back(result_entry);
                }

                // Construct the new arrays for k+1
                if (k < N - 1) {
                    vector<int> new_ppa, new_div;

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
        }
    }
    return results;
}

vector<vector<int> > read_hap(const string &filename) {
    vector<vector<int> > res;
    ifstream ifs(filename); // open the file
    string tempstr;
    while (getline(ifs, tempstr)) {
        stringstream lineStream(tempstr);
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

vector<vector<int>> read_vcf_to_hapmatrix(const string &filename) {
    ifstream infile(filename);
    string line;
    vector<vector<int>> haplotypes;
    bool header_found = false;
    int num_samples = 0;

    while (getline(infile, line)) {
        if (line.rfind("##", 0) == 0) continue;

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
            const string &gt = tokens[i];
            if (gt.size() < 3) continue;
            char a = gt[0], b = gt[2];
            haplotypes[2 * (i - 9)].push_back(a == '1' ? 1 : 0);
            haplotypes[2 * (i - 9) + 1].push_back(b == '1' ? 1 : 0);
        }
    }

    infile.close();
    return haplotypes;
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
        if (line.rfind("##", 0) == 0) continue;
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
            haplotypes[2 * (i - 9)].push_back(gt[0] == '1' ? 1 : 0);
            haplotypes[2 * (i - 9) + 1].push_back(gt[2] == '1' ? 1 : 0);
        }
    }

    file.close();
    return haplotypes;
}

int main(int argc, char* argv[]) {

     int width = 100;
     int height = 10000;

    // string f = "/Users/NguyenLuan/pbwt-comparison/random_data/random_100x100000.vcf";
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <vcf_file_path>" << endl;
        return 1;
    }
    string f = argv[1];
    printf("Read file: %s\n", f.c_str());
    vector<vector<int> > X = read_vcf_to_hapmatrix(f);
    
    if (X.empty() || X[0].empty()) {
        cerr << "Error: VCF data could not be loaded or is empty!\n";
        return 1;
    }
    const int retry = 1;
    const int L = 4;
    printf("Test size: %lu haps x %lu sites\n", X.size(), X[0].size());
    printf("Retry: %d\n", retry);
    printf("Matched length: %d\n", L);

    int nthreads = static_cast<int>(thread::hardware_concurrency());
    printf("Max thread:%d\n\n", nthreads);



        for (int i = 2; i <= nthreads; ++i) {
            // printf("Num thread:%d\n", i);
            signed long int duration_total_build = 0;
            signed long int duration_total_match = 0;
            high_resolution_clock::time_point start, stop;
            vector<vector<int> > matches;
            vector<vector<vector<int> > > res;
            for (int j = 0; j < retry; ++j) {
                start = high_resolution_clock::now();
                res = build_prefix_and_divergence_arrays(X, i);
                stop = high_resolution_clock::now();
                duration_total_build += duration_cast<microseconds>(stop - start).count();
                start = high_resolution_clock::now();
                matches = report_long_matches(X, L, res, i);
                stop = high_resolution_clock::now();
                duration_total_match += duration_cast<microseconds>(stop - start).count();

            }
            // printf("Build time: %ld us\n", duration_total_build / retry);
            // printf("Match time: %ld us\n",  duration_total_match / retry);
            // printf("Total time: %ld us\n\n",  duration_total_match / retry + duration_total_build / retry);
            
            printf("Num thread:%d, total time: %ld us\n", i, duration_total_match / retry + duration_total_build / retry);

    }
}

