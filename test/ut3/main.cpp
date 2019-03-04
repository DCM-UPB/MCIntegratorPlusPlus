#include "mci/MCIntegrator.hpp"
#include "mci/MCISamplingFunctionInterface.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

#include "../common/TestMCIFunctions.hpp"

using namespace std;


int main(){
    const int NMC = 10000;
    const double CORRECT_RESULT = 0.5;

    auto * pdf = new ThreeDimGaussianPDF();
    auto * obs1d = new XSquared();
    auto * obs3d = new XYZSquared();

    MCI * mci = new MCI(3);
    mci->setSeed(5649871);
    mci->addSamplingFunction(pdf);
    mci->addObservable(obs1d);
    mci->addObservable(obs3d);

    // the integral should provide 0.5 as answer!

    double x[3];
    x[0] = 5.; x[1] = -5.; x[2] = 10.;

    double average[4];
    double error[4];

    // this integral will give a wrong answer! This is because the starting point is very bad and initialDecorrelation is skipped (as well as the MRT2step automatic setting)
    mci->setX(x);
    mci->integrate(NMC, average, error, false, false);
    for (int i=0; i<mci->getNObsDim(); ++i) {
        //std::cout << "i " << i << ", average[i] " << average[i] << ", error[i] " << error[i] << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
        assert( fabs(average[i]-CORRECT_RESULT) > 2.*error[i] );
    }
    //std::cout << std::endl;

    // this integral, instead, will provide the right answer
    mci->setX(x);
    mci->integrate(NMC, average, error);
    for (int i=0; i<mci->getNObsDim(); ++i) {
        //std::cout << "i " << i << ", average[i] " << average[i] << ", error[i] " << error[i] << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
        assert( fabs(average[i]-CORRECT_RESULT) < 2.*error[i] );
    }
    //std::cout << std::endl;

    // now, doing an integral without finding again the MRT2step and doing the initialDecorrelation will also result in a correct result
    mci->integrate(NMC, average, error, false, false);
    for (int i=0; i<mci->getNObsDim(); ++i) {
        //std::cout << "i " << i << ", average[i] " << average[i] << ", error[i] " << error[i] << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
        assert( fabs(average[i]-CORRECT_RESULT) < 2.*error[i] );
    }
    //std::cout << std::endl;


    delete pdf;
    delete obs1d;
    delete obs3d;
    delete mci;

    return 0;
}
