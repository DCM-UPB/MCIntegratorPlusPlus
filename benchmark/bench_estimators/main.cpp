#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "mci/Estimators.hpp"

#include "../../test/common/TestMCIFunctions.hpp"
#include "MCIBenchmarks.hpp"

using namespace std;

void run_single_benchmark(const string &label, const double datax[],
                          const int estimatorType /*1 uncorr-1d, 2 block-1d, 3 corr-1d, 4 uncorr-nd, 5 block-nd, 6 corr-nd */,
                          const int NMC, const int ndim, const int nruns)
{
    pair<double, double> result;
    const double time_scale = 1000000000.; //nanoseconds

    result = sample_benchmark_estimators(datax, estimatorType, NMC, ndim, nruns);
    cout << label << ":" << setw(max(1, 20-static_cast<int>(label.length()))) << setfill(' ') << " " << result.first*time_scale/(NMC*ndim) << " +- " << result.second*time_scale/(NMC*ndim) << " nanoseconds" << endl;
}


int main ()
{
    const int NMC = 10000000;
    const int nruns = 10;
    const int estimatorTypes[6] = {1, 2, 3, 4, 5, 6};
    const int ndims[3] = {1, 10, 100};
    const double stepSizes[3] = {1.59, 0.404, 0.12};

    vector<string> labels {"noblock-1D", "20-block-1D", "autoblock-1D", "noblock-ND", "20-block-ND", "autoblock-ND"};

    srand(1337); // consistent random seed

    cout << "=========================================================================================" << endl << endl;
    cout << "Benchmark results (time per sample):" << endl;

    // Estimator benchmark
    for (int i=0; i<3; ++i) { // go through ndims/stepSizes
        const int trueNMC = NMC/ndims[i];
        auto * datax = new double[trueNMC*ndims[i]];
        TestWalk1s testWalk(trueNMC, ndims[i], stepSizes[i]);
        testWalk.generateWalk(datax);

        // the steps sizes were tuned for 0.5+-0.001 acceptance rate (on 1337 seed)
        // so this serves as a little check that nothing was messed up in TestMCIFunctions
        assert(fabs(testWalk.getAcceptanceRate()-0.5) < 0.001); // the steps sizes were tuned for 0.5+-0.001 acceptance rate (on 1337 seed)

        /* // verbose stuff (better decrease NMC to 1000)
           cout << "Acceptance ratio was " << testWalk.getAcceptanceRate() << endl;
           cout << "datax:" << endl;
           for (int imc=0; imc<trueNMC; ++imc) {
           cout << "imc " << imc << endl;
           cout << "x ";
           for (int idm=0; idm<ndims[i]; ++idm) { cout << " " << datax[imc*ndims[i] + idm]; }
           cout << endl;
           }
        */

        for (const int &etype : estimatorTypes) {
            if ( !(ndims[i]>1 && etype<4) ) {
                run_single_benchmark("t/element ( ndim=" + to_string(ndims[i]) + ", " + labels[etype-1] + " )", datax, etype, trueNMC, ndims[i], nruns);
            }
        }

        delete [] datax;
    }
    cout << "=========================================================================================" << endl << endl << endl;

    return 0;
}
