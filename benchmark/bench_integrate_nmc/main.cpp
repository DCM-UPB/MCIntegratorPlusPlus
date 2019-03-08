#include <iomanip>
#include <iostream>
#include <string>

#include "mci/MCIntegrator.hpp"

#include "../../test/common/TestMCIFunctions.hpp"
#include "../common/MCIBenchmarks.hpp"

using namespace std;
using namespace mci;

void run_single_benchmark(const string &label, MCI * mci, const int nruns, const int NMC) {
    pair<double, double> result;
    const double time_scale = 1000000.; //microseconds

    result = sample_benchmark_MCIntegrate(mci, nruns, NMC);
    cout << label << ":" << setw(max(1, 20-static_cast<int>(label.length()))) << setfill(' ') << " " << result.first/NMC*time_scale << " +- " << result.second/NMC*time_scale << " microseconds" << endl;
}

int main () {
    // debug settings
    //const int NMC[5] = {1000, 1000, 1000, 1000, 1000};
    //const int nruns[5] = {2, 2, 2, 2, 2};

    // benchmark settings
    const int NMC[5] = {10000000, 1000000, 100000, 10000, 1000};
    const int nruns[5] = {5, 50, 500, 5000, 50000};

    auto * pdf = new ThreeDimGaussianPDF();
    auto * obs1 = new XND(3);
    auto * obs2 = new XSquared();
    auto * obs3 = new XYZSquared();

    MCI mci(3);
    mci.setSeed(5649871);
    mci.addSamplingFunction(pdf);
    mci.addObservable(*obs1, 0); // use simple accu (i.e. no error), no skipping
    mci.addObservable(*obs2, 1); // use auto-blocking, no skipping
    mci.addObservable(*obs3, 5, 2); // fixed blocksize 5, skip every 2nd step (i.e. "effective" block size of 10)

    cout << mci.getNObsDim() << endl;
    double avg[mci.getNObsDim()],err[mci.getNObsDim()];
    const double mrt2step[3] = {1.85, 1.85, 1.85};
    mci.setMRT2Step(mrt2step);
    mci.integrate(5000, avg, err, true, true); // decorrelate

    cout << "=========================================================================================" << endl << endl;
    cout << "Benchmark results (time per sample):" << endl;

    // MCIntegrate benchmark
    for (int inmc=0; inmc<5; ++inmc) {
        run_single_benchmark("t/step (" + std::to_string(NMC[inmc]) + " steps)", &mci, nruns[inmc], NMC[inmc]);
    }
    cout << "=========================================================================================" << endl << endl << endl;

    delete obs3;
    delete obs2;
    delete obs1;
    delete pdf;

    return 0;
}
