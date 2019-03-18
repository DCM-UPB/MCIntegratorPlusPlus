#include <iomanip>
#include <iostream>
#include <memory>
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
    // benchmark settings
    const int nset = 6;
    const int NMC[nset] = {25000, 50000, 100000, 200000, 400000, 400000};
    const int ndims[nset] = {1, 4, 16, 64, 256, 1024};
    const int nruns[nset] = {5120, 1280, 320, 80, 20, 5};
    const double mrt2steps[nset] = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0};

    std::vector< std::unique_ptr<MCI> > mcis;
    for (int nd : ndims) {
        mcis.push_back( std::make_unique<MCI>(nd) );
    }

    for (int i=0; i<nset; ++i) {
        MCI * mci = mcis[i].get();
        int nd = mci->getNDim();
        ExpNDPDF pdf(nd);
        XND obs(nd); // use non-updateable XND, because the observable is not expensive enough for selective updates

        mci->setTrialMove(MoveType::Vec); // default (i.e. uniform) single-index move

        mci->setSeed(1337);
        mci->addSamplingFunction(pdf);
        mci->addObservable(obs, 0, 1); // use simple accu (i.e. no error), no skipping

        double avg[nd], err[nd];
        for (int j=0; j<nd; ++j) { mci->setX(j, j%2==0 ? 0.1 : -0.05 ); }
        mci->setMRT2Step(mrt2steps[i]);

        mci->setNdecorrelationSteps(500000);
        mci->integrate(0, avg, err, false, true); // warmup&decorrelate
        cout << "acceptance rate " << mci->getAcceptanceRate() << endl;
    }

    cout << "=========================================================================================" << endl << endl;
    cout << "Benchmark results (time per step and dimension):" << endl;

    // MCIntegrate benchmark
    for (int inmc=0; inmc<nset; ++inmc) {
        run_single_benchmark("t/step (" + std::to_string(ndims[inmc]) + " dim)", *(mcis[inmc]), nruns[inmc], NMC[inmc]);
    }
    cout << "=========================================================================================" << endl << endl << endl;

    return 0;
}
