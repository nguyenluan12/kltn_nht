
#define MAIN_H
#include <vector>
using namespace std;

vector<vector<vector<int>>> build_prefix_and_divergence_arrays(const vector<vector<int> > &X,int thread);

void algorithm4(int M, int k, int b,
                 vector<int> &sorted_k,  vector<int> &start_k,
                 const vector<int> &sorted_k_b,  vector<int> &start_k_b);

void algorithm3(int id_start, int id_end, const int index_k_b[],
                 vector<int> &sorted_k,  vector<int> &start_k,
                 const vector<int> &sorted_k_b,  vector<int> &start_k_b);

vector<vector<int> > report_long_matches(const vector<vector<int> > &X, const int L, vector<vector<vector<int>>> res,int thread);

#endif //MAIN_H
