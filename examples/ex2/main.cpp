#include "mpi.h"
#include <cmath>
#include <cmath>
#include <fstream>
#include <iostream>

#include "mci/MCIntegrator.hpp"
#include "mci/MPIMCI.hpp"


// Observable functions
class Parabola: public MCIObservableFunctionInterface{
public:
    explicit Parabola(const int ndim): MCIObservableFunctionInterface(ndim, 1) {}

    void observableFunction(const double in[], double out[]) override{
        out[0] = 4.*in[0] - in[0]*in[0];
    }
};

class NormalizedParabola: public MCIObservableFunctionInterface{
public:
    explicit NormalizedParabola(const int ndim): MCIObservableFunctionInterface(ndim, 1) {}

    void observableFunction(const double in[], double out[]) override{
        out[0] = (4. - in[0]) * 5.;
        if (std::signbit(in[0])) { out[0] = -out[0];}
    }
};



// Sampling function
// the 48 is for normalization (even if not strictly necessary)
class NormalizedLine: public MCISamplingFunctionInterface{
public:
    explicit NormalizedLine(const int ndim): MCISamplingFunctionInterface(ndim, 1) {}

    void samplingFunction(const double in[], double protovalue[]) override{
        protovalue[0] = 0.2 * fabs(in[0]);
    }

    double getAcceptance(const double protoold[], const double protonew[]) override{
        return protonew[0] / protoold[0];
    }
};



int main() {
    using namespace std;

    int myrank = MPIMCI::init(); // first run custom MPI init

    if (myrank == 0) {
        // intro
        cout << "We want to compute the integral" << endl;
        cout << "    Integral[-1:3] dx (4x - x^2)" << endl << endl;
        cout << "Notice that we can re-write the integral as" << endl;
        cout << "    Integral[-1:3] dx (5*sign(x)*(4-x) * |x|/5)" << endl;
        cout << "where g(x)=|x|/5 is a candidate sampling function since it's positive on the domain [-1:3] and its integral on the domain is equal to 1." << endl << endl;

        cout << "We start by initializing MPI and setting the MCI:" << endl;
    }

    const int ndim = 1;
    MCI mci(ndim);

    if (myrank == 0) { cout << "ndim = " << mci.getNDim() << endl;}

    // set the integration range to [-1:3]
    mci.setIRange(-1., 3.);

    if (myrank == 0) { cout << "irange = [ " << mci.getLBound(0) << " ; " << mci.getUBound(0) << " ]" << endl;}

    // initial walker position
    double initpos[ndim];
    initpos[0] = -0.5;
    mci.setX(initpos);

    if (myrank == 0) { cout << "initial walker position = " << mci.getX(0) << endl;}

    // initial MRT2 step
    double step[ndim];
    step[0] = 0.5;
    mci.setMRT2Step(step);

    if (myrank == 0) { cout << "MRT2 step = " << mci.getMRT2Step(0) << endl;}

    // target acceptance rate
    mci.setTargetAcceptanceRate(0.7);

    if (myrank == 0) {
        cout << "Acceptance rate = " << mci.getTargetAcceptanceRate() << endl;
        cout << endl << endl;

        // first way of integrating
        cout << "We first compute the integral without setting a sampling function." << endl;
        cout << "    f(x) = (4x - x^2) " << endl;
        cout << "    g(x) = - " << endl << endl;
    }

    // observable
    MCIObservableFunctionInterface * obs = new Parabola(ndim);
    mci.addObservable(obs);

    if (myrank == 0) {
        cout << "Number of observables set = " << mci.getNObs() << endl;
        cout << "Dimension of observables set = " << mci.getNObsDim() << endl;
        // sampling function
        cout << "Number of sampling function set = " << mci.getNSampF() << endl;
    }

    // integrate
    const int Nmc = 1000000;
    double average[mci.getNObsDim()];
    double error[mci.getNObsDim()];

    // ! set fixed amount of findMRT2 and decorrelation steps    !
    // ! this is very important for efficient parallel execution !
    mci.setNfindMRT2Iterations(50);
    mci.setNdecorrelationSteps(5000);

    MPIMCI::integrate(&mci, Nmc, average, error);

    if (myrank == 0) {
        cout << "The integral gives as result = " << average[0] << "   +-   " << error[0] << endl;
        cout << "--------------------------------------------------------" << endl << endl;

        // second way of integrating
        cout << "Now we compute the integral using a sampling function." << endl;
        cout << "    f(x) = 5*sign(x)*(4-x) " << endl;
        cout << "    g(x) = |x|/5 " << endl << endl;
    }

    // observable
    delete obs;
    obs = new NormalizedParabola(ndim);
    mci.clearObservables();  // we first remove the old observable
    mci.addObservable(obs);

    // sampling function
    MCISamplingFunctionInterface * sf = new NormalizedLine(ndim);
    mci.addSamplingFunction(sf);

    if (myrank == 0) {
        cout << "Number of observables set = " << mci.getNObs() << endl;
        cout << "Dimension of observables set = " << mci.getNObsDim() << endl;
        cout << "Number of sampling function set = " << mci.getNSampF() << endl;
    }

    // integrate
    MPIMCI::integrate(&mci, Nmc, average, error);

    if (myrank == 0) {
        cout << "The integral gives as result = " << average[0] << "   +-   " << error[0] << endl;
        cout << "--------------------------------------------------------" << endl << endl;


        // final comments
        cout << "Using a sampling function in this case gives worse performance. In fact, the error bar is larger." << endl;
        cout << "This implies that the variance of the re-factored f(x) written for introducing a sampling function, is larger than the original f(x)." << endl;
    }

    // deallocate per-thread allocations
    delete sf;
    delete obs;

    // finalize MPI
    MPIMCI::finalize();

    // end
    return 0;
}
