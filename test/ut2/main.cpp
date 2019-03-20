#include "mci/MCIntegrator.hpp"
#include "mci/SamplingFunctionInterface.hpp"
#include "mci/OrthoPeriodicDomain.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

#include "../common/TestMCIFunctions.hpp"

using namespace std;
using namespace mci;

int main(){
    MCI mci(3);

    // first let's test some set/get-like methods

    // unbound domain (default)
    for (int i=0; i<3; ++i) { assert(mci.getX(i) == 0.); } // x should be 0 per default
    mci.moveX();
    bool flag_same = true;
    for (int i=0; i<3; ++i) { // check for changes
        if (mci.getX(i) != 0.) { flag_same = false; }
    }
    assert(!flag_same); // at least one index should have changed

    // pbc domain (set by manual construction)
    mci.setDomain( OrthoPeriodicDomain(3, -5., 5.) );
    mci.centerX(); // recenter X within domain
    for (int i=0; i<3; ++i) { assert(mci.getX(i) == 0.); }
    mci.newRandomX();
    for (int i=0; i<3; ++i) { // this is extremely unlikely to fail
        assert(mci.getX(i) != 0.);
    }

    mci.resetDomain(); // restore unbound domain


    // now we test actual integration
    const int NMC = 10000;
    const double CORRECT_RESULT = 0.5;

    // the integral of this pdf/obs combination should provide 0.5 as answer!
    ThreeDimGaussianPDF pdf;
    XSquared obs;

    mci.addSamplingFunction(pdf);
    mci.addObservable(obs);

    // choose very bad initial position
    double x[3];
    x[0] = 5.; x[1] = -5.; x[2] = 10.;

    // set the magic seed
    mci.setSeed(5649871);

    double average;
    double error;

    // this integral will give a wrong answer! This is because the starting point is very bad and initialDecorrelation is skipped (as well as the MRT2step automatic setting)
    mci.setX(x);
    mci.integrate(NMC, &average, &error, false, false);
    //std::cout << "average " << average << ", error " << error << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
    assert( fabs(average-CORRECT_RESULT) > 2.*error );

    // this integral, instead, will provide the right answer
    mci.setX(x);
    mci.integrate(NMC, &average, &error);
    //std::cout << "average " << average << ", error " << error << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
    assert( fabs(average-CORRECT_RESULT) < 2.*error );

    // now, doing an integral without finding again the MRT2step and doing the initialDecorrelation will also result in a correct result
    mci.integrate(NMC, &average, &error, false, false);
    //std::cout << "average " << average << ", error " << error << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
    assert( fabs(average-CORRECT_RESULT) < 2.*error );

    // now setting a fixed integration range
    mci.setIRange(-5., 5.); // implicitly sets pbc domain
    mci.integrate(NMC, &average, &error, false, false);
    //std::cout << "average " << average << ", error " << error << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
    assert( fabs(average-CORRECT_RESULT) < 2.*error );

    // and now also without sampling function
    mci.clearSamplingFunctions();
    mci.clearObservables();
    GaussXSquared obs_nopdf; // is XSquared multiplied by Gaussian
    mci.addObservable(obs_nopdf, 1, 1, false, false); // don't use correlated estimator (but would work as well)
    mci.integrate(NMC, &average, &error, false, false);
    //std::cout << "average " << average << ", error " << error << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
    assert( fabs(average-CORRECT_RESULT) < 2.*error );


    return 0;
}
