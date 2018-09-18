#include "mpi.h"
#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>

#include "MCIntegrator.hpp"
#include "MPIMCI.hpp"


// Observable functions
class Parabola: public MCIObservableFunctionInterface{
public:
    Parabola(const int &ndim): MCIObservableFunctionInterface(ndim, 1) {}

    void observableFunction(const double * in, double * out){
        out[0] = 4.*in[0] - in[0]*in[0];
    }
};

class NormalizedParabola: public MCIObservableFunctionInterface{
public:
    NormalizedParabola(const int &ndim): MCIObservableFunctionInterface(ndim, 1) {}

    void observableFunction(const double * in, double * out){
        out[0] = (4. - in[0]) * 5.;
        if (std::signbit(in[0])) out[0] = -out[0];
    }
};



// Sampling function
// the 48 is for normalization (even if not strictly necessary)
class NormalizedLine: public MCISamplingFunctionInterface{
public:
    NormalizedLine(const int &ndim): MCISamplingFunctionInterface(ndim, 1) {}

    void samplingFunction(const double * in, double * protovalue){
        protovalue[0] = 0.2 * abs(in[0]);
    }

    double getAcceptance(const double * protoold, const double * protonew){
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
    MCI * mci = new MCI(ndim);

    if (myrank == 0) cout << "ndim = " << mci->getNDim() << endl;

    // set the integration range to [-1:3]
    double ** irange = new double*[ndim];
    irange[0] = new double[2];
    irange[0][0] = -1.;
    irange[0][1] = 3.;
    mci->setIRange(irange);

    if (myrank == 0) cout << "irange = [ " << mci->getIRange(0, 0) << " ; " << mci->getIRange(0, 1) << " ]" << endl;

    // initial walker position
    double * initpos = new double[ndim];
    initpos[0] = -0.5;
    mci->setX(initpos);

    if (myrank == 0) cout << "initial walker position = " << mci->getX(0) << endl;


    // initial MRT2 step
    double * step = new double[ndim];
    step[0] = 0.5;
    mci->setMRT2Step(step);

    if (myrank == 0) cout << "MRT2 step = " << mci->getMRT2Step(0) << endl;

    // target acceptance rate
    double * targetacceptrate = new double[1];
    targetacceptrate[0] = 0.7;
    mci->setTargetAcceptanceRate(targetacceptrate);

    if (myrank == 0) {
        cout << "Acceptance rate = " << mci->getTargetAcceptanceRate() << endl;
        cout << endl << endl;

        // first way of integrating
        cout << "We first compute the integral without setting a sampling function." << endl;
        cout << "    f(x) = (4x - x^2) " << endl;
        cout << "    g(x) = - " << endl << endl;
    }

    // observable
    MCIObservableFunctionInterface * obs = new Parabola(ndim);
    mci->addObservable(obs);

    if (myrank == 0) {
        cout << "Number of observables set = " << mci->getNObs() << endl;
        // sampling function
        cout << "Number of sampling function set = " << mci->getNSampF() << endl;
    }

    // integrate
    const long Nmc = 1000000;
    double * average;
    double * error;
    if (myrank == 0) { // allocate only for root
        average = new double[mci->getNObsDim()];
        error = new double[mci->getNObsDim()];
    }
    MPIMCI::integrate(mci, Nmc, average, error);

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
    mci->clearObservables();  // we first remove the old observable
    mci->addObservable(obs);

    if (myrank == 0) cout << "Number of observables set = " << mci->getNObs() << endl;


    // sampling function
    MCISamplingFunctionInterface * sf = new NormalizedLine(ndim);
    mci->addSamplingFunction(sf);

    if (myrank == 0) cout << "Number of sampling function set = " << mci->getNSampF() << endl;


    // integrate
    MPIMCI::integrate(mci, Nmc, average, error);

    // deallocate per-thread allocations
    delete sf;

    delete obs;

    delete[] targetacceptrate;

    delete[] step;

    delete[] initpos;

    delete[] irange[0];
    delete[] irange;

    MPIMCI::finalize(mci); // also deletes all mci

    // from now on just root thread
    if (myrank == 0) {
        cout << "The integral gives as result = " << average[0] << "   +-   " << error[0] << endl;
        cout << "--------------------------------------------------------" << endl << endl;


        // final comments
        cout << "Using a sampling function in this case gives worse performance. In fact, the error bar is larger." << endl;
        cout << "This implies that the variance of the re-factored f(x) written for introducing a sampling function, is larger than the original f(x)." << endl;

        // deallocate
        delete[] average;
        delete[] error;
    }

    // end
    return 0;
}
