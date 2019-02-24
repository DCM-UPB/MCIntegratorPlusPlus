#include "mci/MCIntegrator.hpp"
#include "mci/MCISamplingFunctionInterface.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

#include "../common/TestMCIFunctions.hpp"

using namespace std;


int main(){
    const long NMC = 10000;
    const double CORRECT_RESULT = 0.5;

    auto * pdf = new ThreeDimGaussianPDF();
    auto * obs1d = new XSquared();
    auto * obs3d = new XYZSquared();

    MCI * mci = new MCI(3);
    mci->setSeed(5649871);
    mci->addSamplingFunction(pdf);
    mci->addObservable(obs1d);
    mci->addObservable(obs3d);
    mci->setNBlocks(0); // use auto-blocking

    // the integral should provide 0.5 as answer!

    auto * x = new double[3];
    x[0] = 5.; x[1] = -5.; x[2] = 10.;

    auto * average = new double[4];
    auto * error = new double[4];

    // this integral will give a wrong answer! This is because the starting point is very bad and initialDecorrelation is skipped (as well as the MRT2step automatic setting)
    mci->setX(x);
    mci->integrate(NMC, average, error, false, false);
    for (int i=0; i<mci->getNObsDim(); ++i) {
        assert( fabs(average[i]-CORRECT_RESULT) > 2.*error[i] );
    }

    // this integral, instead, will provide the right answer
    mci->setX(x);
    mci->integrate(NMC, average, error);
    for (int i=0; i<mci->getNObsDim(); ++i) {
        assert( fabs(average[i]-CORRECT_RESULT) < 2.*error[i] );
    }

    // now, doing an integral without finding again the MRT2step and doing the initialDecorrelation will also result in a correct result
    mci->integrate(NMC, average, error, false, false);
    for (int i=0; i<mci->getNObsDim(); ++i) {
        assert( fabs(average[i]-CORRECT_RESULT) < 2.*error[i] );
    }


    delete pdf;
    delete obs1d;
    delete obs3d;
    delete mci;
    delete [] x;
    delete [] average;
    delete [] error;

    return 0;
}
