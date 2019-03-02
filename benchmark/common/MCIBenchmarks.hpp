#include <cmath>
#include <iostream>
#include <tuple>
#include <stdexcept>

#include "Timer.hpp"
#include "mci/Estimators.hpp"
#include "mci/MCIntegrator.hpp"

inline double benchmark_MCIntegrate(MCI * mci, const int NMC) {
    Timer timer(1.);
    double average[mci->getNObs()];
    double error[mci->getNObs()];
    std::fill(average, average+mci->getNObs(), 0.);
    std::fill(error, error+mci->getNObs(), 0.);

    timer.reset();
    mci->integrate(NMC, average, error, false, false);
    return timer.elapsed();
}

inline std::pair<double, double> sample_benchmark_MCIntegrate(MCI * mci, const int nruns, const int NMC) {
    double times[nruns];
    double mean = 0., err = 0.;

    for (int i=0; i<nruns; ++i) {
        times[i] = benchmark_MCIntegrate(mci, NMC);
        mean += times[i];
    }
    mean /= nruns;
    for (int i=0; i<nruns; ++i) { err += pow(times[i]-mean, 2);
}
    err /= (nruns-1)*nruns; // variance of the mean
    err = sqrt(err); // standard error of the mean

    return std::pair<double, double>(mean, err);
}


inline double benchmark_estimators(const double * datax,
                                   const int estimatorType /* 1 uncorr-1d, 2 block-1d, 3 corr-1d, 4 uncorr-nd, 5 block-nd, 6 corr-nd */,
                                   const int NMC, const int ndim)
{
    const int nblocks = 20;
    Timer timer(1.);
    double average[ndim];
    double error[ndim];
    std::fill(average, average+ndim, 0.);
    std::fill(error, error+ndim, 0.);

    switch(estimatorType)
        {
        case 1:
            if (ndim>1) { throw std::invalid_argument("OneDimUncorrelatedEstimator requested, but ndim>1"); }
            else {timer.reset(); mci::OneDimUncorrelatedEstimator(NMC, datax, average[0], error[0]); }
            break;

        case 2:
            if (ndim>1) { throw std::invalid_argument("OneDimBlockEstimator requested, but ndim>1"); }
            else {timer.reset(); mci::OneDimBlockEstimator(NMC, datax, nblocks, average[0], error[0]); }
            break;

        case 3:
            if (ndim>1) { throw std::invalid_argument("OneDimCorrelatedEstimator requested, but ndim>1"); }
            else {timer.reset(); mci::OneDimCorrelatedEstimator(NMC, datax, average[0], error[0]); }
            break;

        case 4:
            timer.reset();
            mci::MultiDimUncorrelatedEstimator(NMC, ndim, datax, average, error);
            break;

        case 5:
            timer.reset();
            mci::MultiDimBlockEstimator(NMC, ndim, datax, nblocks, average, error);
            break;

        case (6):
            timer.reset();
            mci::MultiDimCorrelatedEstimator(NMC, ndim, datax, average, error);
            break;

        default:
            throw std::invalid_argument("Invalid estimator typeID.");
        }

    return timer.elapsed();
}


inline std::pair<double, double> sample_benchmark_estimators(const double * datax, const int estimatorType, const int NMC, const int ndim, const int nruns) {
    double times[nruns];
    double mean = 0., err = 0.;

    for (int i=0; i<nruns; ++i) {
        times[i] = benchmark_estimators(datax, estimatorType, NMC, ndim);
        mean += times[i];
    }
    mean /= nruns;
    for (int i=0; i<nruns; ++i) { err += pow(times[i]-mean, 2); }
    err /= (nruns-1)*nruns; // variance of the mean
    err = sqrt(err); // standard error of the mean

    return std::pair<double, double>(mean, err);
}

