#include "mci/Estimators.hpp"
#include "mci/MCIBlockAccumulator.hpp"
#include "mci/MCIFullAccumulator.hpp"
#include "mci/MCISimpleAccumulator.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <tuple>

#include "../common/TestMCIFunctions.hpp"


using namespace std;

void reportAvgErr1D(const string &label, double avg, double err)
{
    cout << "- " << label << endl;
    cout << "     avg = " << avg << "     error = " << err << endl << endl;
}

void reportAvgErrND(const string &label, int ndim, double avg[], double err[])
{
    cout << "- " << label << endl;
    for (int i=0; i<ndim; ++i) {
        cout << "     avg" << i << " = " << avg[i] << "     error" << i << " = " << err[i] << endl;
    }
    cout << endl;
}

void arrayAvgND(int N1, int N2, const double in[] /*N1*N2*/, double out[] /*N2*/)
{   // calculate average of length N2
    std::fill(out, out+N2, 0.);
    for (int i=0; i<N1; ++i) {
        for (int j=0; j<N2; ++j) {
            out[j] += in[i*N2 + j];
        }
    }
    for (int i=0; i<N2; ++i) { out[i] /= N1; }
}

void assertArraysEqual(int ndim, const double arr1[], const double arr2[], double tol = 0.)
{
    if (tol > 0.) {
        for (int i=0; i<ndim; ++i) { assert( fabs(arr1[i] - arr2[i]) < tol ); }
    } else {
        for (int i=0; i<ndim; ++i) { assert( arr1[i] == arr2[i] ); } // the check above may fail in this case
    }
}


void assertAccuAveragesEqual(MCIAccumulatorInterface * accu1, MCIAccumulatorInterface * accu2, double tol = 0.)
{   // check that averages of contained data are equal within tol
    const int nobs = accu1->getNObs();
    assert(nobs == accu2->getNObs());

    double avg1[nobs], avg2[nobs];
    arrayAvgND(accu1->getNStore(), nobs, accu1->getData(), avg1);
    arrayAvgND(accu2->getNStore(), nobs, accu2->getData(), avg2);
    assertArraysEqual(nobs, avg1, avg2, tol);
}

void assertAccuResetted(MCIAccumulatorInterface * accu)
{   // check that accu is in clean reset state (allocated/deallocated doesn't matter)
    assert(accu->getStepIndex() == 0);
    assert(accu->isClean());
    assert(!accu->isFinalized());
    for (int i=0; i<accu->getNData(); ++i) { assert(accu->getData()[i] == 0.); } 
}

void assertAccuDeallocated(MCIAccumulatorInterface * accu)
{   // check that the accu is in proper deallocated state
    assert(!accu->isAllocated());
    assert(accu->getNSteps() == 0);
    assert(accu->getNAccu() == 0);
    assert(accu->getNStore() == 0);
    assert(accu->getNData() == 0);

    assertAccuResetted(accu);
}

void assertAccuAllocated(MCIAccumulatorInterface * accu, int Nmc)
{   // check that the accu is in allocated state (not necessarily reset state)
    assert(accu->isAllocated());
    assert(accu->getNSteps() == Nmc);
    assert(accu->getNAccu() > 0);
    assert(accu->getNStore() > 0);
    assert(accu->getNData() > 0);
    assert(accu->getNData() == accu->getNStore() * accu->getNObs());
}

void assertAccuFinalized(MCIAccumulatorInterface * accu, int Nmc)
{   // check that the accu is properly finalized
    assert(accu->isAllocated());
    assert(!accu->isClean());
    assert(accu->isFinalized());
    assert(accu->getStepIndex() == Nmc);
}


void accumulateData(MCIAccumulatorInterface * accu, int Nmc, int ndim, const double datax[], const bool datacc[])
{   // simulated MC observable accumulation
    bool flags_xchanged[ndim];
    std::fill(flags_xchanged, flags_xchanged+ndim, true);
    for (int i=0; i<Nmc; ++i) {
        accu->accumulate(datax+i*ndim, datacc[i], flags_xchanged);
    }
    accu->finalize();
}

void checkAccumulator(MCIAccumulatorInterface * accu, int Nmc, int ndim, const double datax[], const bool datacc[],
                      double tol /* tolerance for avg */, bool verbose = false /* to enable printout */)
{
    // we expect walker-dim == obs-dim
    assert(accu->getNObs() == ndim);
    assert(accu->getNDim() == ndim);

    // verify that the accumulator is uninitialized
    assertAccuDeallocated(accu);

    // now allocate
    accu->allocate(Nmc);
    assertAccuAllocated(accu, Nmc); // allocated
    assertAccuResetted(accu); // but still clean

    // accumulate the data in pseudo MC loop
    double storedData[accu->getNData()]; // to store away obs data
    accumulateData(accu, Nmc, ndim, datax, datacc);
    assertAccuFinalized(accu, Nmc);

    // copy the stored data
    const double * const dataptr = accu->getData(); // we acquire a read-only pointer to data
    std::copy(dataptr, dataptr+accu->getNData(), storedData);

    // now do the same after reset
    accu->reset();
    assertAccuResetted(accu);
    accumulateData(accu, Nmc, ndim, datax, datacc);
    assertArraysEqual(accu->getNData(), storedData, accu->getData()); // check that we get the same result

    // now do the same after reallocation
    accu->deallocate();
    assertAccuDeallocated(accu);
    accu->allocate(Nmc);
    assertAccuAllocated(accu, Nmc);
    accu->allocate(Nmc); // do it twice on purpose
    accumulateData(accu, Nmc, ndim, datax, datacc);
    assertArraysEqual(accu->getNData(), storedData, accu->getData()); // check that we get the same result

    // finally check that average calculated from the data in
    // the accumulator matches the reference average within tol
    double refAvg[ndim];
    arrayAvgND(Nmc, ndim, datax, refAvg);

    if (accu->getNStore() > 1) {
        double avg[ndim];
        arrayAvgND(accu->getNStore(), ndim, accu->getData(), avg);
        for (int i=0; i<ndim; ++i) {
            assert( fabs(avg[i] - refAvg[i]) < tol );
            if (verbose) {
                cout << "avg" << i << " " << avg[i] << " refAvg" << i << " " << refAvg[i] << endl;
            }
        }
    }
    else {
        for (int i=0; i<ndim; ++i) {
            assert( fabs(accu->getData()[i] - refAvg[i]) < tol );
            if (verbose) {
                cout << "avg" << i << " " << accu->getData()[i] << " refAvg" << i << " " << refAvg[i] << endl;
            }
        }
    }
}

int main(){
    using namespace std;

    bool verbose = false;
    //verbose = true; // uncomment for output

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
    if (verbose) { cout << testWalk.getAcceptanceRate() << endl; }

    // calculate reference averages
    double refAvg[nd];
    arrayAvgND(Nmc, nd, xND, refAvg);

    if (verbose) {
        cout << "Reference Average: " << endl << "avgND =";
        for (double avg : refAvg) { cout << " " << avg; }
        cout << endl << endl;
    }


    // --- check 1D estimators ---
    if (verbose) { cout << endl << "1-dimensional versions of Estimators:" << endl << endl; }

    // helper input for 1D
    double x1D[Nmc];

    for (int i=0; i<nd; ++i) {
        // copy subdata into continous array
        for (int j=0; j<Nmc; ++j) { x1D[j] = xND[j*nd + i]; }

        // perform check for correct averages
        double avg1D = 0., err1D = 0.;

        mci::OneDimUncorrelatedEstimator(Nmc, x1D, avg1D, err1D);
        if (verbose) { reportAvgErr1D("UncorrelatedEstimator()", avg1D, err1D); }
        assert(avg1D == refAvg[i]); // these should be absolutely identical

        mci::OneDimBlockEstimator(Nmc, x1D, nblocks, avg1D, err1D);
        if (verbose) { reportAvgErr1D("BlockEstimator()", avg1D, err1D); }
        assert(fabs(avg1D - refAvg[i]) < EXTRA_TINY ); // these should be nearly identical

        mci::OneDimCorrelatedEstimator(Nmc, x1D, avg1D, err1D);
        if (verbose) { reportAvgErr1D("CorrelatedEstimator()", avg1D, err1D); }
        assert(fabs(avg1D - refAvg[i]) < SMALL ); // this difference should be small

        // the following is currently true with the selected seed
        // and valid changes to the estimator are very unlikely
        // to break the check (just take a better seed then ;-) )
        assert(fabs(avg1D - refAvg[i]) < 3*err1D );
    }


    // --- check ND estimators ---
    if (verbose) { cout << endl << "Multidimensional versions of Estimators:" << endl << endl; }

    // result arrays to pass
    double avgND[nd], errND[nd];

    mci::MultiDimUncorrelatedEstimator(Nmc, nd, xND, avgND, errND);
    if (verbose) { reportAvgErrND("MultiDimUncorrelatedEstimator()", nd, avgND, errND); }
    for (int i=0; i<nd; ++i) { assert(avgND[i] == refAvg[i]); } // these should be absolutely identical

    mci::MultiDimBlockEstimator(Nmc, nd, xND, nblocks, avgND, errND);
    if (verbose) { reportAvgErrND("MultiDimBlockEstimator()", nd, avgND, errND); }
    for (int i=0; i<nd; ++i) { assert(fabs(avgND[i] - refAvg[i]) < EXTRA_TINY ); } // these should be nearly identical

    mci::MultiDimCorrelatedEstimator(Nmc, nd, xND, avgND, errND);
    if (verbose) { reportAvgErrND("MultiDimCorrelatedEstimator()", nd, avgND, errND); }
    for (int i=0; i<nd; ++i) {
        assert(fabs(avgND[i] - refAvg[i]) < SMALL ); // this difference should be small
        assert(fabs(avgND[i] - refAvg[i]) < 3*errND[i] ); // like in the 1D case
    }


    // --- check accumulators ---
    if (verbose) { cout << endl << "Now using accumulator classes to store data:" << endl << endl; }
    auto * obsfun = new XND(nd); // n-dimensional position observable
    auto * simpleAccu = new MCISimpleAccumulator(obsfun, 1);
    auto * simpleAccuSkip2 = new MCISimpleAccumulator(obsfun, 2);
    auto * blockAccu = new MCIBlockAccumulator(obsfun, 1, 10);
    auto * blockAccuSkip2 = new MCIBlockAccumulator(obsfun, 2, 5);
    auto * fullAccu = new MCIFullAccumulator(obsfun, 1);
    auto * fullAccuSkip2 = new MCIFullAccumulator(obsfun, 2);

    vector<pair< MCIAccumulatorInterface *, string > > accuList;
    accuList.emplace_back(simpleAccu, "simpleAccu" );
    accuList.emplace_back(blockAccu, "blockAccu" );
    accuList.emplace_back(fullAccu, "fullAccu" );
    accuList.emplace_back(simpleAccuSkip2, "simpleAccuSkip2" );
    accuList.emplace_back(blockAccuSkip2, "blockAccuSkip2" );
    accuList.emplace_back(fullAccuSkip2, "fullAccuSkip2" );

    for (auto & accuTup : accuList) {
        if (verbose) { cout << endl << "Checking accumulator " << accuTup.second << " ..." << endl; }
        checkAccumulator(accuTup.first, Nmc, nd, xND, accepted, SMALL, verbose);
    }

    // check that values with same skip level are very equal
    assertAccuAveragesEqual(simpleAccu, blockAccu, EXTRA_TINY);
    assertAccuAveragesEqual(simpleAccu, fullAccu, EXTRA_TINY);

    assertAccuAveragesEqual(simpleAccuSkip2, blockAccuSkip2, EXTRA_TINY);
    assertAccuAveragesEqual(simpleAccuSkip2, fullAccuSkip2, EXTRA_TINY);

    for (auto & accuTup : accuList) {
        delete accuTup.first; // free memory
    }
    delete obsfun;
}
