#include "mci/MCIntegrator.hpp"

#include <cassert>

#include "../common/TestMCIFunctions.hpp"

using namespace std;
using namespace mci;

int main()
{
    const int NMC = 16384; // use a power of 2 this time, to let MCI use MJBlocker
    const double CORRECT_RESULT = 0.5;

    ThreeDimGaussianPDF pdf;
    XSquared obs;

    MCI mci(3);
    mci.setSeed(5649871);
    mci.addSamplingFunction(pdf);
    mci.addObservable(obs);

    // the integral should provide 0.5 as answer!

    double x[3]{5., -5., 10.}; // bad starting point
    double average;
    double error;

    // configure a fixed amount steps to use in pre-sampling
    mci.setNfindMRT2Iterations(20);
    mci.setNdecorrelationSteps(2000);

    // this integral should provide the right answer (after findMRT2step&decorrelation phase)
    mci.setX(x);
    mci.integrate(NMC, &average, &error, true, true);
    //std::cout << "average " << average << ", error " << error << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
    assert(fabs(average - CORRECT_RESULT) < 2.*error);

    // now, doing an integral without finding again the MRT2step and doing the initialDecorrelation will also result in a correct result
    mci.integrate(NMC, &average, &error, false, false);
    //std::cout << "average " << average << ", error " << error << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
    assert(fabs(average - CORRECT_RESULT) < 2.*error);

    // and using fixed blocking also gives the same result
    mci.clearObservables();
    mci.addObservable(obs, 16); // blocks of size 16
    mci.integrate(NMC, &average, &error, false, false);
    //std::cout << "average " << average << ", error " << error << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
    assert(fabs(average - CORRECT_RESULT) < 2.*error);

    // and half the block size with skipping every second step, should be similar again
    mci.clearObservables();
    mci.addObservable(obs, 8, 2); // blocksize 4, nskip 2 (i.e. "effective" blocksize of 16)
    mci.integrate(NMC, &average, &error, false, false);
    //std::cout << "average " << average << ", error " << error << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
    assert(fabs(average - CORRECT_RESULT) < 2.*error);

    // The previous two integrations implicitly used UncorrelatedEstimator (because blocksize>1).
    // CorrelatedEstimator (is MJBlocker, because NMC is power of 2) should yield similar result.

    // first set it by boolean arguments
    mci.clearObservables();
    mci.addObservable(obs, 8, 2, false /*flag_equil*/, true /*flag_correlated*/); // forcing correlated estimator, although the other settings would "imply" uncorrelated samples
    mci.integrate(NMC, &average, &error, false, false);
    //std::cout << "average " << average << ", error " << error << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
    assert(fabs(average - CORRECT_RESULT) < 2.*error);

    // now via enumerator
    mci.clearObservables();
    mci.addObservable(obs, 8, 2, false /*flag_equil*/, EstimatorType::Correlated);
    mci.integrate(NMC, &average, &error, false, false);
    //std::cout << "average " << average << ", error " << error << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
    assert(fabs(average - CORRECT_RESULT) < 2.*error);


    return 0;
}
