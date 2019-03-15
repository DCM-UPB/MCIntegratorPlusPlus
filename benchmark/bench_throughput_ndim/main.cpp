#include <iomanip>
#include <iostream>
#include <string>

#include "mci/MCIntegrator.hpp"

#include "../../test/common/TestMCIFunctions.hpp"
#include "../common/MCIBenchmarks.hpp"

#include <memory>

using namespace std;
using namespace mci;

void run_single_benchmark(const string &label, MCI &mci, const int nruns, const int NMC) {
    pair<double, double> result;
    const double time_scale = 1000000.; //microseconds
    const double full_scale = (time_scale / NMC) / mci.getNDim(); // time per step per dim

    result = sample_benchmark_MCIntegrate(mci, nruns, NMC);
    cout << label << ":" << setw(max(1, 20-static_cast<int>(label.length()))) << setfill(' ') << " " << result.first*full_scale << " +- " << result.second*full_scale << " microseconds" << endl;
}

int main () {
    // debug settings
    //const int NMC[5] = {1000, 1000, 1000, 1000, 1000};
    //const int nruns[5] = {2, 2, 2, 2, 2};

    // benchmark settings
    const int NMC = 50000;
    const int ndims[6] = {1, 4, 16, 64, 256, 1024};
    const int nruns[6] = {5120, 1280, 320, 80, 20, 5};
    const double mrt2steps[6] = {3.0, 1.35, 0.6, 0.3, 0.15, 0.07};

    std::vector< std::unique_ptr<MCI> > mcis;
    for (int nd : ndims) {
        mcis.push_back( std::unique_ptr<MCI>(new MCI(nd)) );
    }

    for (int i=0; i<6; ++i) {
        MCI * mci = mcis[i].get();
        int nd = mci->getNDim();
        ExpNDPDF pdf(nd);
        XND obs(nd);

        mci->setSeed(1337);
        mci->addSamplingFunction(pdf);
        mci->addObservable(obs, 0, 1); // use simple accu (i.e. no error), no skipping

        double avg[nd], err[nd];
        for (int j=0; j<nd; ++j) { mci->setX(j, j%2==0 ? 0.1 : -0.05 ); }
        mci->setMRT2Step(mrt2steps[i]);

        mci->integrate(500000, avg, err, false, false); // warmup&decorrelate
        cout << "acceptance rate " << mci->getAcceptanceRate() << endl;
    }

    cout << "=========================================================================================" << endl << endl;
    cout << "Benchmark results (time per step and dimension):" << endl;

    // MCIntegrate benchmark
    for (int inmc=0; inmc<6; ++inmc) {
        run_single_benchmark("t/step (" + std::to_string(ndims[inmc]) + " dim)", *(mcis[inmc]), nruns[inmc], NMC);
    }
    cout << "=========================================================================================" << endl << endl << endl;

    return 0;
}
