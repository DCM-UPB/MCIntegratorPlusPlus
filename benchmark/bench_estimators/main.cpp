#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "mci/Estimators.hpp"

#include "../../test/common/TestMCIFunctions.hpp"
#include "../common/MCIBenchmarks.hpp"

using namespace std;
using namespace mci;

void run_single_benchmark(const string &label, const double datax[],
                          BenchEstim estimatorType,
                          const int NMC, const int ndim, const int nruns)
{
    pair<double, double> result;
    const double time_scale = 1000000000.; //nanoseconds

    result = sample_benchmark_estimators(datax, estimatorType, NMC, ndim, nruns);
    cout << label << ":" << setw(max(1, 20 - static_cast<int>(label.length()))) << setfill(' ') << " " << result.first*time_scale/(NMC*ndim) << " +- " << result.second*time_scale/(NMC*ndim) << " nanoseconds" << endl;
}


int main()
{
    const bool flag_debug = false;
    const bool verbose = flag_debug; // some debug printout
    const int NMC = flag_debug ? 4096 : 8388608; // use power of two to allow MJBlocker
    const int nruns = 10;
    const int ntypes = 7;
    const BenchEstim estimatorTypes[ntypes] = {BenchEstim::Uncorr1D, BenchEstim::Block1D, BenchEstim::CorrFC1D,
                                               BenchEstim::UncorrND, BenchEstim::BlockND, BenchEstim::CorrFCND, BenchEstim::CorrMJND};
    const int ndims[3] = {1, 16, 128};
    const double stepSizes[3] = {1.59, 0.313, 0.1053};

    vector<string> labels{"noblock-1D", "500K-block-1D", "autoblock-FC-1D", "noblock-ND", "500K-block-ND", "autoblock-FC-ND", "autoblock-MJ-ND"};

    srand(1337); // consistent random seed

    cout << "=========================================================================================" << endl << endl;
    cout << "Benchmark results (time per sample and dimension):" << endl;

    // Estimator benchmark
    for (int i = 0; i < 3; ++i) { // go through ndims/stepSizes
        const int trueNMC = NMC/ndims[i];
        auto * datax = new double[trueNMC*ndims[i]];
        TestWalk1s testWalk(trueNMC, ndims[i], stepSizes[i], 1. /*all-particle steps*/);
        testWalk.generateWalk(datax);

        if (verbose) { // verbose stuff (better decrease NMC to 1000)
            cout << "Acceptance ratio was " << testWalk.getAcceptanceRate() << endl;
            cout << "datax:" << endl;
            for (int imc = 0; imc < trueNMC; ++imc) {
                cout << "imc " << imc << endl;
                cout << "x ";
                for (int idm = 0; idm < ndims[i]; ++idm) { cout << " " << datax[imc*ndims[i] + idm]; }
                cout << endl;
            }
        }

        // the steps sizes were tuned for 0.5+-0.01 acceptance rate (on 1337 seed)
        // so this serves as a little check that nothing was messed up in TestMCIFunctions
        assert(fabs(testWalk.getAcceptanceRate() - 0.5) < 0.01);

        for (int j = 0; j < ntypes; ++j) {
            if (!(ndims[i] > 1 && j < 3)) {
                run_single_benchmark("t/element ( ndim=" + to_string(ndims[i]) + ", " + labels[j] + " )", datax, estimatorTypes[j], trueNMC, ndims[i], nruns);
            }
        }

        delete[] datax;
    }
    cout << "=========================================================================================" << endl << endl << endl;

    return 0;
}
