#include <iomanip>
#include <iostream>

#include "mci/MCIntegrator.hpp"

#include "../../test/common/TestMCIFunctions.hpp"
#include "../common/MCIBenchmarks.hpp"

using namespace std;
using namespace mci;

void run_single_benchmark(const string &label, MCI &mci, const int nruns, const int64_t NMC)
{
    pair<double, double> result;

    result = sample_benchmark_MCIntegrate(mci, nruns, NMC);
    cout << label << ":" << setw(max(1, 20 - static_cast<int>(label.length()))) << setfill(' ') << " " << result.first << " +- " << result.second << " seconds" << endl;
}

int main()
{
    // benchmark settings
    const int64_t NMC = 3000000000; // 3G steps, enough to overflow 32 bit int
    const int nruns = 1;

    Exp1DPDF pdf;
    X1D obs;

    MCI mci(1);
    mci.setSeed(1337);
    mci.addSamplingFunction(pdf);
    mci.addObservable(obs, 20, 1); // use block accu with blocksize20, no skipping (i.e. a bit over 1GB RAM)

    double avg, err;
    const double mrt2step = 3.185;
    mci.setMRT2Step(&mrt2step);

    mci.integrate(100000, &avg, &err, false, false); // warmup&decorrelate
    cout << "avg " << avg << ", err " << err << endl;
    cout << "acceptance rate " << mci.getAcceptanceRate() << endl;

    cout << "=========================================================================================" << endl << endl;
    cout << "Benchmark results (time per 3 Giga-Steps):" << endl;

    // MCIntegrate benchmark
    run_single_benchmark("t/run (" + std::to_string(NMC) + " steps)", mci, nruns, NMC);
    cout << "=========================================================================================" << endl << endl << endl;
    cout << "acceptance rate " << mci.getAcceptanceRate() << endl;

    return 0;
}
