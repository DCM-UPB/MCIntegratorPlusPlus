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

bool isStepAccepted(double oldWFVal, double newWFVal)
{   // standard VMC acceptance criterion
    if (oldWFVal == 0) {
        return true;
    }
    if (newWFVal == 0) {
        return false;
    }
    const double threshold = (newWFVal*newWFVal)/(oldWFVal*oldWFVal);
    if (threshold >= 1.) {
        return true;
    }
    return ( rand()*(1.0 / RAND_MAX) <= threshold );
}

double calc1sOrbitalWFVal(const double * position, const int ndim)
{   // product of 1s orbitals in 1D
    double wfval = 0.;
    for (int i=0; i<ndim; ++i) {
        wfval += fabs(position[i]);
    }
    return exp(-wfval);
}

void generate1sOrbitalPosition(const double * oldPosition, double * newPosition, const int ndim)
{
    double oldWFVal = calc1sOrbitalWFVal(oldPosition, ndim);
    
    bool accepted = false;
    do {
        for (int i=0; i<ndim; ++i) {
            newPosition[i] = oldPosition[i] + 0.2*(rand()*(1.0 / RAND_MAX) - 0.5);
        }
        const double newWFVal = calc1sOrbitalWFVal(newPosition, ndim);
        accepted = isStepAccepted(oldWFVal, newWFVal);
    } while (!accepted);
}

void generate1sOrbitalWalk(double * datax, const int NMC, const int ndim)
{
    for (int j=0; j<ndim; ++j) { datax[j] = 0.; } // set first step to 0
    for (int i=1; i<NMC; ++i) {
        generate1sOrbitalPosition(datax+(i-1)*ndim, datax+i*ndim, ndim);
    }
}


void run_single_benchmark(const string &label, const int estimatorType /*1 uncorr, 2 block, 3 corr */,const double * datax, const int NMC, const int ndim, const int nruns)
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
    const int estimatorTypes[3] = {1, 2, 3};
    const int ndims[3] = {1, 10, 100};

    vector<string> labels {"noblock", "20-block", "autoblock"};

    srand(1337); // consistent random seed

    cout << "=========================================================================================" << endl << endl;
    cout << "Benchmark results (time per sample):" << endl;

    // Estimator benchmark
    for (const int ndim : ndims) {
        const int trueNMC = NMC/ndim;
        auto * datax = new double[trueNMC*ndim];
        generate1sOrbitalWalk(datax, trueNMC, ndim);

        for (const int etype : estimatorTypes) {
            run_single_benchmark("t/element ( ndim=" + to_string(ndim) + ", " + labels[etype-1] + " )", etype, datax, trueNMC, ndim, nruns);
        }

        delete [] datax;
    }
    cout << "=========================================================================================" << endl << endl << endl;

    return 0;
}
