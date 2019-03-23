#include <iomanip>
#include <iostream>
#include <string>

#include "mci/MCIntegrator.hpp"

#include "../../test/common/TestMCIFunctions.hpp"
#include "../common/MCIBenchmarks.hpp"

using namespace std;
using namespace mci;

void run_single_benchmark(const string &label, MCI &mci, const int nruns, const int NMC) {
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
    const int NMC[5] = {1000, 10000, 100000, 1000000, 10000000};
    const int nruns[5] = {50000, 5000, 500, 50, 5};

    Exp1DPDF pdf;
    X1D obs;

    MCI mci(1);
    mci.setSeed(1337);
    mci.addSamplingFunction(pdf);
    mci.addObservable(obs, 0, 1); // use simple accu (i.e. no error), no skipping

    double avg, err;
    const double mrt2step = 3.185;
    mci.setMRT2Step(&mrt2step);

    mci.integrate(100000, &avg, &err, false, false); // warmup&decorrelate
    cout << "avg " << avg <<  ", err " << err << endl;
    cout << "acceptance rate " << mci.getAcceptanceRate() << endl;

    cout << "=========================================================================================" << endl << endl;
    cout << "Benchmark results (time per sample):" << endl;

    // MCIntegrate benchmark
    for (int inmc=0; inmc<5; ++inmc) {
        run_single_benchmark("t/step (" + std::to_string(NMC[inmc]) + " steps)", mci, nruns[inmc], NMC[inmc]);
    }
    cout << "=========================================================================================" << endl << endl << endl;

    return 0;
}
