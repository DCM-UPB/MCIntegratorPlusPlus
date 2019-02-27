#include <iostream>
#include <cmath>
#include <tuple>

#include "mci/MCIntegrator.hpp"
#include "mci/Estimators.hpp"
#include "Timer.hpp"

double benchmark_MCIntegrate(MCI * mci, const long &NMC) {
    Timer timer(1.);
    double average[mci->getNObs()];
    double error[mci->getNObs()];

    timer.reset();
    mci->integrate(NMC, average, error, false, false);
    return timer.elapsed();
}

std::pair<double, double> sample_benchmark_MCIntegrate(MCI * mci, const int nruns, const long &NMC) {
    double times[nruns];
    double mean = 0., err = 0.;

    for (int i=0; i<nruns; ++i) {
        times[i] = benchmark_MCIntegrate(mci, NMC);
        mean += times[i];
    }
    mean /= nruns;
    for (int i=0; i<nruns; ++i) err += pow(times[i]-mean, 2);
    err /= (nruns-1)*nruns; // variance of the mean
    err = sqrt(err); // standard error of the mean

    return std::pair<double, double>(mean, err);
}


double benchmark_estimators(const double * datax, const int estimatorType /* 1 uncorr, 2 block, 3 corr */, const long &NMC, const int ndim) {
    const int nblocks = 20;
    Timer timer(1.);
    double average[ndim];
    double error[ndim];

    timer.reset(); // we contain two if checks in the result, but shouldnt matter
    if (estimatorType == 1) {
        if (ndim == 1) { mci::UncorrelatedEstimator(NMC, datax, *average, *error); }
        else { mci::MultiDimUncorrelatedEstimator(NMC, ndim, datax, average, error); }
    }
    else if (estimatorType == 2) {
        if (ndim == 1) { mci::BlockEstimator(NMC, datax, nblocks, *average, *error); }
        else { mci::MultiDimBlockEstimator(NMC, ndim, datax, nblocks, average, error); }
    }
    else if (estimatorType == 3) {
        if (ndim == 1) { mci::CorrelatedEstimator(NMC, datax, *average, *error); }
        else { mci::MultiDimCorrelatedEstimator(NMC, ndim, datax, average, error); }
    }
    const double time = timer.elapsed();

    return time;
}

std::pair<double, double> sample_benchmark_estimators(const double * datax, const int estimatorType, const long &NMC, const int ndim, const int nruns) {
    double times[nruns];
    double mean = 0., err = 0.;

    for (int i=0; i<nruns; ++i) {
        times[i] = benchmark_estimators(datax, estimatorType, NMC, ndim);
        mean += times[i];
    }
    mean /= nruns;
    for (int i=0; i<nruns; ++i) err += pow(times[i]-mean, 2);
    err /= (nruns-1)*nruns; // variance of the mean
    err = sqrt(err); // standard error of the mean

    return std::pair<double, double>(mean, err);
}

