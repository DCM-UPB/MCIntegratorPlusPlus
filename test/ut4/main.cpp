#include "mci/MCIntegrator.hpp"
#include "mci/SamplingFunctionInterface.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

#include "../common/TestMCIFunctions.hpp"

using namespace std;
using namespace mci;

int main(){
    const int NMC = 10000;
    const double CORRECT_RESULT = 0.5;

    auto * pdf = new ThreeDimGaussianPDF();
    auto * obs = new XSquared();

    MCI mci(3);
    mci.setSeed(5649871);
    mci.addSamplingFunction(pdf);
    mci.addObservable(*obs);

    // the integral should provide 0.5 as answer!

    double x[3];
    x[0] = 5.; x[1] = -5.; x[2] = 10.;

    double average;
    double error;

    // configure a fixed amount steps to use in pre-sampling
    mci.setNfindMRT2Iterations(10); // 10 iterations
    mci.setNdecorrelationSteps(1000);

    // this integral will give a wrong answer! This is because the starting point is very bad and initialDecorrelation is skipped (as well as the MRT2step automatic setting)
    mci.setX(x);
    mci.integrate(NMC, &average, &error, false, false);
    //std::cout << "average " << average << ", error " << error << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
    assert( fabs(average-CORRECT_RESULT) > 2.*error );

    // this integral, instead, will provide the right answer
    mci.setX(x);
    mci.integrate(NMC, &average, &error, true, true);
    //std::cout << "average " << average << ", error " << error << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
    assert( fabs(average-CORRECT_RESULT) < 2.*error );

    // now, doing an integral without finding again the MRT2step and doing the initialDecorrelation will also result in a correct result
    mci.integrate(NMC, &average, &error, false, false);
    //std::cout << "average " << average << ", error " << error << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
    assert( fabs(average-CORRECT_RESULT) < 2.*error );

    // and using fixed blocking also gives the same result
    mci.clearObservables();
    mci.addObservable(*obs, 10); // blocks of size 10
    mci.integrate(NMC, &average, &error, false, false);
    //std::cout << "average " << average << ", error " << error << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
    assert( fabs(average-CORRECT_RESULT) < 2.*error );

    // and half the block size with skipping every second step, should be similar again
    mci.clearObservables();
    mci.addObservable(*obs, 5, 2); // blocksize 5, nskip 2 (i.e. "effective" blocksize of 10)
    mci.integrate(NMC, &average, &error, false, false);
    //std::cout << "average " << average << ", error " << error << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
    assert( fabs(average-CORRECT_RESULT) < 2.*error );

    delete pdf;
    delete obs;

    return 0;
}
