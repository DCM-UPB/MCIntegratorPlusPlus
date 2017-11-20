#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>

#include "MCIntegrator.hpp"



// Observable functions
class Parabola: public MCIObservableFunctionInterface{
public:
    Parabola(const int &ndim): MCIObservableFunctionInterface(ndim, 1) {}
protected:
    void observableFunction(const double * in, double * out){
        out[0] = 4.*in[0] - in[0]*in[0];
    }
};

class NormalizedParabola: public MCIObservableFunctionInterface{
public:
    NormalizedParabola(const int &ndim): MCIObservableFunctionInterface(ndim, 1) {}
protected:
    void observableFunction(const double * in, double * out){
        out[0] = (4. - in[0]) * 5.;
        if (signbit(in[0])) out[0] = -out[0];
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
    
    double getAcceptance(){
        return this->getProtoNew(0) / this->getProtoOld(0);
    }
};



int main() {
    using namespace std;
    
    
    // intro
    cout << "We want to compute the integral" << endl;
    cout << "    Integral[-1:3] dx (4x - x^2)" << endl << endl;
    cout << "Notice that we can re-write the integral as" << endl;
    cout << "    Integral[-1:3] dx (5*sign(x)*(4-x) * |x|/5)" << endl;
    cout << "where g(x)=|x|/5 is a candidate sampling function since it's positive on the domain [-1:3] and its integral on the domain is equal to 1." << endl << endl;
    
    cout << "We start by setting the MCI:" << endl;
    
    

    // declare a 1-dimensional integrator
    const int ndim = 1;
    MCI * mci = new MCI(ndim);
    
    cout << "ndim = " << mci->getNDim() << endl;
    
    
    
    // set the integration range to [-1:3]
    double ** irange = new double*[ndim];
    irange[0] = new double[2];
    irange[0][0] = -1.;
    irange[0][1] = 3.;
    mci->setIRange(irange);
    
    cout << "irange = [ " << mci->getIRange(0, 0) << " ; " << mci->getIRange(0, 1) << " ]" << endl;
    
    
    
    // initial walker position
    double * initpos = new double[ndim];
    initpos[0] = -0.5;
    mci->setX(initpos);
    
    cout << "initial walker position = " << mci->getX(0) << endl;
    
    
    
    // initial MRT2 step
    double * step = new double[ndim];
    step[0] = 0.5;
    mci->setMRT2Step(step);
   
    cout << "MRT2 step = " << mci->getMRT2Step(0) << endl;
    
    
    
    // target acceptance rate
    double * targetacceptrate = new double[1];
    targetacceptrate[0] = 0.7;
    mci->setTargetAcceptanceRate(targetacceptrate);
    
    cout << "Acceptance rate = " << mci->getTargetAcceptanceRate() << endl;
    cout << endl << endl;
   
   
   
   
   
    // first way of integrating
    cout << "We first compute the integral without setting a sampling function." << endl;
    cout << "    f(x) = (4x - x^2) " << endl;
    cout << "    g(x) = - " << endl << endl;
    
    
    
    // observable
    MCIObservableFunctionInterface * obs = new Parabola(ndim);
    mci->addObservable(obs);
    
    cout << "Number of observables set = " << mci->getNObs() << endl;
    
    
    
    // sampling function
    cout << "Number of sampling function set = " << mci->getNSampF() << endl;
    
    
    
    // integrate
    const long Nmc = 1000000;
    double * average = new double[mci->getNObsDim()];
    double * error = new double[mci->getNObsDim()];
    mci->integrate(Nmc, average, error);
    
    cout << "The integral gives as result = " << average[0] << "   +-   " << error[0] << endl;
    cout << "--------------------------------------------------------" << endl << endl;
    
    
    
    
    
    
    // second way of integrating
    cout << "Now we compute the integral using a sampling function." << endl;
    cout << "    f(x) = 5*sign(x)*(4-x) " << endl;
    cout << "    g(x) = |x|/5 " << endl << endl;
    
   
   
   
    // observable
    obs = new NormalizedParabola(ndim);
    mci->clearObservables();  // we first remove the old observable
    mci->addObservable(obs);
    
    cout << "Number of observables set = " << mci->getNObs() << endl;
    
    
    
    // sampling function
    MCISamplingFunctionInterface * sf = new NormalizedLine(ndim);
    mci->addSamplingFunction(sf);

    cout << "Number of sampling function set = " << mci->getNSampF() << endl;
   
   
   
    // integrate
    mci->integrate(Nmc, average, error);
    
    cout << "The integral gives as result = " << average[0] << "   +-   " << error[0] << endl;
    cout << "--------------------------------------------------------" << endl << endl;
    
    
    
    // final comments
    cout << "Using a sampling function in this case gives worse performance. In fact, the error bar is larger." << endl;
    cout << "This implies that the variance of the re-factored f(x) written for introducing a sampling function, is larger than the original f(x)." << endl;
   
   
      
    // deallocate
    delete[] average;
    delete[] error;
    
    delete sf;
    
    delete obs;
    
    delete[] targetacceptrate;
    
    delete[] step;
    
    delete[] initpos;
   
    delete[] irange[0];
    delete[] irange;

    delete mci;
    


    // end
    return 0;
}
