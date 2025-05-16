// basic file operations
#include <iostream>
#include <fstream>
using namespace std;

int main() {
    ofstream myfile;
    // int hap=5004;
    // int sites=2370256;
    int hap=1000;
    int sites=100000;
    for (int i = 0; i < 4; ++i) {
        myfile.open(to_string(hap)+"hap x "+to_string(sites)+"site x"+to_string(i)+".hap",std::ios::out | std::ios::trunc);
        for (int j = 0; j < sites; ++j) {
            for (int k = 0; k < hap; ++k) {
                myfile << rand() % 2;
                myfile << " ";
            }
            myfile << "\n";
        }
        myfile.close();
    }
}
