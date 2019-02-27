#include <iostream>
#include <iomanip>
#include <string> 

#include "mci/MCIntegrator.hpp"

#include "MCIBenchmarks.hpp"
#include "../../test/common/TestMCIFunctions.hpp"

using namespace std;

void run_single_benchmark(const string &label, MCI * mci, const int nruns, const int NMC) {
    pair<double, double> result;
    const double time_scale = 1000000.; //microseconds

    result = sample_benchmark_MCIntegrate(mci, nruns, NMC);
    cout << label << ":" << setw(max(1, 20-(int)label.length())) << setfill(' ') << " " << result.first/NMC*time_scale << " +- " << result.second/NMC*time_scale << " microseconds" << endl;
}

int main (void) {
    const int NMC[5] = {10000000, 1000000, 100000, 10000, 1000};
    const int nruns[5] = {3, 30, 300, 3000, 30000};

    ThreeDimGaussianPDF * pdf = new ThreeDimGaussianPDF();
    XSquared * obs = new XSquared();

    MCI mci(3);
    mci.setSeed(5649871);
    mci.addSamplingFunction(pdf);
    mci.addObservable(obs);
    mci.setNBlocks(50); // don't use auto-blocking


    double avg[3],err[3];
    const double mrt2step[3] = {1.85, 1.85, 1.85};
    mci.setMRT2Step(mrt2step);
    mci.integrate(5000, avg, err, false, false); // decorrelate

    cout << "=========================================================================================" << endl << endl;
    cout << "Benchmark results (time per sample):" << endl;

    // MCIntegrate benchmark
    for (int inmc=0; inmc<5; ++inmc) {
        run_single_benchmark("t/step (" + std::to_string(NMC[inmc]) + " steps)", &mci, nruns[inmc], NMC[inmc]);
    }
    cout << "=========================================================================================" << endl << endl << endl;

    delete obs;
    delete pdf;

    return 0;
}

