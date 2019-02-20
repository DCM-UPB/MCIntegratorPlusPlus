#include "mci/MCIntegrator.hpp"
#include "mci/MCISamplingFunctionInterface.hpp"

#include <iostream>
#include <math.h>
#include <assert.h>

#include "TestMCIFunctions.hpp"

using namespace std;


int main(){
    const long NMC = 10000;
    const double CORRECT_RESULT = 0.5;

    ThreeDimGaussianPDF * pdf = new ThreeDimGaussianPDF();
    XSquared * obs = new XSquared();

    MCI * mci = new MCI(3);
    mci->setSeed(5649871);
    mci->addSamplingFunction(pdf);
    mci->addObservable(obs);
    mci->setNBlocks(0); // use auto-blocking

    // the integral should provide 0.5 as answer!

    double * x = new double[3];
    x[0] = 5.; x[1] = -5.; x[2] = 10.;

    double * average = new double;
    double * error = new double;

    // configure the number of steps to use in pre-sampling
    mci->setNfindMRT2steps(10);
    mci->setNdecorrelationSteps(1000);

    // this integral will give a wrong answer! This is because the starting point is very bad and initialDecorrelation is skipped (as well as the MRT2step automatic setting)
    mci->setX(x);
    mci->integrate(NMC, average, error, false, false);
    assert( fabs(average[0]-CORRECT_RESULT) > 2.*error[0] );

    // this integral, instead, will provide the right answer
    mci->setX(x);
    mci->integrate(NMC, average, error, true, true);
    assert( fabs(average[0]-CORRECT_RESULT) < 2.*error[0] );

    // now, doing an integral without finding again the MRT2step and doing the initialDecorrelation will also result in a correct result
    mci->integrate(NMC, average, error, false, false);
    assert( fabs(average[0]-CORRECT_RESULT) < 2.*error[0] );

    // and using fixed blocking also gives the same result
    mci->setNBlocks(15);
    mci->integrate(NMC, average, error, false, false);
    assert( fabs(average[0]-CORRECT_RESULT) < 2.*error[0] );


    delete pdf;
    delete obs;
    delete mci;
    delete [] x;
    delete average;
    delete error;

    return 0;
}
