#ifndef MCI_MCIBENCHMARKS_HPP
#define MCI_MCIBENCHMARKS_HPP

#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <tuple>

#include "Timer.hpp"
#include "mci/Estimators.hpp"
#include "mci/MCIntegrator.hpp"

inline std::pair<double, double> sample_benchmark(const std::function< double () > &run_benchmark /*all parameters bound*/, const int nruns)
{
    double times[nruns];
    double mean = 0., err = 0.;

    for (int i=0; i<nruns; ++i) { // run given benchmark nruns times
        times[i] = run_benchmark();
        mean += times[i];
    }
    mean /= nruns;
    if (nruns>1) {
        for (int i=0; i<nruns; ++i) { err += pow(times[i]-mean, 2); }
        err /= (nruns-1.)*nruns; // variance of the mean
        err = sqrt(err); // standard error of the mean
    }

    return std::pair<double, double>(mean, err);
}

inline double benchmark_MCIntegrate(mci::MCI &mci, const int64_t NMC) {
    Timer timer(1.);
    double average[mci.getNObsDim()];
    double error[mci.getNObsDim()];
    std::fill(average, average+mci.getNObsDim(), 0.);
    std::fill(error, error+mci.getNObsDim(), 0.);

    timer.reset();
    mci.integrate(NMC, average, error, false, false);
    return timer.elapsed();
}

inline std::pair<double, double> sample_benchmark_MCIntegrate(mci::MCI & mci, const int nruns, const int64_t NMC) {
    return sample_benchmark([&] { return benchmark_MCIntegrate(mci, NMC); }, nruns);
}

// supported estimator configurations
enum class BenchEstim{Uncorr1D,
                           Block1D,
                           CorrFC1D,
                           UncorrND,
                           BlockND,
                           CorrFCND,
                           CorrMJND};

inline double benchmark_estimators(const double datax[],
                            BenchEstim estimatorType,
                            const int64_t NMC, const int ndim)
{
    const int blocksize = 16; // for fixed blocking estimator
    Timer timer(1.);
    double average[ndim];
    double error[ndim];
    std::fill(average, average+ndim, 0.);
    std::fill(error, error+ndim, 0.);

    switch(estimatorType)
        {
        case BenchEstim::Uncorr1D:
            if (ndim>1) { throw std::invalid_argument("OneDimUncorrelatedEstimator requested, but ndim>1"); }
            else {timer.reset(); mci::OneDimUncorrelatedEstimator(NMC, datax, average[0], error[0]); }
            break;

        case BenchEstim::Block1D:
            if (ndim>1) { throw std::invalid_argument("OneDimBlockEstimator requested, but ndim>1"); }
            else {timer.reset(); mci::OneDimBlockEstimator(NMC, datax, NMC/blocksize, average[0], error[0]); }
            break;

        case BenchEstim::CorrFC1D:
            if (ndim>1) { throw std::invalid_argument("OneDimFCBlockerEstimator requested, but ndim>1"); }
            else {timer.reset(); mci::OneDimFCBlockerEstimator(NMC, datax, average[0], error[0]); }
            break;

        case BenchEstim::UncorrND:
            timer.reset();
            mci::MultiDimUncorrelatedEstimator(NMC, ndim, datax, average, error);
            break;

        case BenchEstim::BlockND:
            timer.reset();
            mci::MultiDimBlockEstimator(NMC, ndim, datax, NMC/blocksize, average, error);
            break;

        case (BenchEstim::CorrFCND):
            timer.reset();
            mci::MultiDimFCBlockerEstimator(NMC, ndim, datax, average, error);
            break;

        case (BenchEstim::CorrMJND):
            timer.reset();
            mci::MJBlockerEstimator(NMC, ndim, datax, average, error);
            break;


        default:
            throw std::invalid_argument("Invalid enumerator.");
        }

    return timer.elapsed();
}


inline std::pair<double, double> sample_benchmark_estimators(const double datax[], BenchEstim estimatorType, const int64_t NMC, const int ndim, const int nruns) {
    return sample_benchmark([&] { return benchmark_estimators(datax, estimatorType, NMC, ndim); }, nruns);
}

#endif
