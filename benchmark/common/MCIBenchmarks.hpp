#include <iostream>
#include <cmath>
#include <tuple>

#include "mci/MCIntegrator.hpp"
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

    const std::pair<double, double> result(mean, err);
    return result;
}
