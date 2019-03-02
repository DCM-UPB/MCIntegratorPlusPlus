#include "mci/Estimators.hpp"
#include "mci/MCISimpleAccumulator.hpp"
#include "mci/MCIBlockAccumulator.hpp"
#include "mci/MCIFullAccumulator.hpp"

#include <cassert>
#include <cmath>
#include <random>
#include <iostream>
#include <algorithm>
#include <string>

#include "../common/TestMCIFunctions.hpp"


using namespace std;

void reportAvgErr1D(const string &label, double avg, double err)
{
    cout << "- " << label << endl;
    cout << "     avg = " << avg << "     error = " << err << endl << endl;
}

void reportAvgErrND(const string &label, int ndim, double * avg, double * err)
{
    cout << "- " << label << endl;
    for (int i=0; i<ndim; ++i) {
        cout << "     avg" << i << " = " << avg[i] << "     error" << i << " = " << err[i] << endl;
    }
    cout << endl;
}

int main(){
    using namespace std;

    const double SMALL = 0.01;
    const double EXTRA_TINY = 0.0000000001; // 1e-10

    const int Nmc = 10000;
    const int nd = 2;
    const int ndata = Nmc*nd;

    const int nblocks=10; // blocks to use for fixed-block estimation

    // generate random walk
    double xND[ndata];
    bool accepted[Nmc]; // tells us which steps are new ones
    srand(1337); // seed standard random engine
    TestWalk1s testWalk(Nmc, nd, 1.0); // 2-particle walk in 1-dim 1s orbital
    testWalk.generateWalk(xND, accepted);
    cout << testWalk.getAcceptanceRate() << endl;

    // calculate reference average
    double refAvg[nd];
    std::fill(refAvg, refAvg+nd, 0.);
    for (int i=0; i<Nmc; ++i) {
        for (int j=0; j<nd; ++j) {
            refAvg[j] += xND[i*nd + j];
        }
    }
    for (int i=0; i<nd; ++i) { refAvg[i] /= Nmc; }

    cout << "Reference Average: " << endl << "avgND =";
    for (int i=0; i<nd; ++i) { cout << " " << refAvg[i]; }
    cout << endl << endl;

    // result arrays to pass
    double avgND[nd], errND[nd];
    // helper input for 1D
    double x1D[Nmc];

    // check 1D estimators
    for (int i=0; i<nd; ++i) {
        // copy subdata into continous array
        for (int j=0; j<Nmc; ++j) { x1D[j] = xND[j*nd + i]; }

        // perform check for correct averages
        double avg1D = 0., err1D = 0.;

        mci::OneDimUncorrelatedEstimator(Nmc, x1D, avg1D, err1D);
        reportAvgErr1D("UncorrelatedEstimator()", avg1D, err1D);
        assert(avg1D == refAvg[i]); // these should be absolutely identical

        mci::OneDimBlockEstimator(Nmc, x1D, nblocks, avg1D, err1D);
        reportAvgErr1D("BlockEstimator()", avg1D, err1D);
        assert(fabs(avg1D - refAvg[i]) < EXTRA_TINY ); // these should be nearly identical

        mci::OneDimCorrelatedEstimator(Nmc, x1D, avg1D, err1D);
        reportAvgErr1D("CorrelatedEstimator()", avg1D, err1D);
        assert(fabs(avg1D - refAvg[i]) < SMALL ); // this difference should be small
        // the following is currently true with the selected seed
        // and valid changes to the estimator are very unlikely
        // to break the check (just take a better seed then ;-) )
        assert(fabs(avg1D - refAvg[i]) < 3*err1D );
    }

    cout << endl << "Multidimensional version of Estimators" << endl << endl;

    mci::MultiDimUncorrelatedEstimator(Nmc, nd, xND, avgND, errND);
    reportAvgErrND("MultiDimUncorrelatedEstimator()", nd, avgND, errND);
    for (int i=0; i<nd; ++i) { assert(avgND[i] == refAvg[i]); } // these should be absolutely identical

    mci::MultiDimBlockEstimator(Nmc, nd, xND, nblocks, avgND, errND);
    reportAvgErrND("MultiDimBlockEstimator()", nd, avgND, errND);
    for (int i=0; i<nd; ++i) { assert(fabs(avgND[i] - refAvg[i]) < EXTRA_TINY ); } // these should be nearly identical

    mci::MultiDimCorrelatedEstimator(Nmc, nd, xND, avgND, errND);
    reportAvgErrND("MultiDimCorrelatedEstimator()", nd, avgND, errND);
    for (int i=0; i<nd; ++i) {
        assert(fabs(avgND[i] - refAvg[i]) < SMALL ); // these should be about identical
        assert(fabs(avgND[i] - refAvg[i]) < 3*errND[i] ); // like in the 1D case
    }

    cout << endl << "Checking that we get the same result" << endl << endl;
}
