#include <fstream>

/*
  Scales range limited set of random seeds to the full range of
  unsigned 64 bit integer, which the std rgen expects.

  Provide a unique set of 10000 integers between 0 and 42007935
  in a text file set1.txt. Format expected as you get it from
  random.org without comma delimiters and unordered.
*/

int main() {
    using namespace std;

    ifstream infile;
    ofstream outfile;
    infile.open("set1.txt");
    outfile.open("set1_out.txt");

    uint_fast64_t factor = 65537;
    factor *= 6700417;
    for (int i=0; i<10000; ++i) {
        uint_fast64_t seed;
        infile >> seed;
        seed *= factor;
        outfile << seed;
        outfile << " ";
    }
    infile.close();
    outfile.close();

    return 0;
}
