#include "MCIntegrator.hpp"
#include "MCIObservableFunctionInterface.hpp"
#include "FeedForwardNeuralNetwork.hpp"

#include <iostream>

/*
Build an example for MCI++:

compute the integral

         \int_{-10}^{+10} dx x^2 nn(x)^2
(nn(x) is a normalized neural network) with MC and a non-MC method

do the same with

         \int_{-10}^{+10} dx \laplacian^2 nn(x)^2
*/





class NN2: public MCIObservableFunctionInterface{
private:
    FeedForwardNeuralNetwork * _ffnn;

public:
    NN2(FeedForwardNeuralNetwork * ffnn): MCIObservableFunctionInterface(1, 1){
        _ffnn = ffnn;
    }

    virtual void observableFunction(const double * in, double *out){
        _ffnn->setInput(1, in);
        _ffnn->FFPropagate();
        const double v = _ffnn->getOutput(1);
        out[0] = v*v;
    }

};



class NN2X2: public MCIObservableFunctionInterface{
private:
    FeedForwardNeuralNetwork * _ffnn;

public:
    NN2X2(FeedForwardNeuralNetwork * ffnn): MCIObservableFunctionInterface(1, 1){
        _ffnn = ffnn;
    }

    virtual void observableFunction(const double * in, double *out){
        _ffnn->setInput(1, in);
        _ffnn->FFPropagate();
        const double v = _ffnn->getOutput(1);
        out[0] = v*v*in[0]*in[0];
    }

};


class X2: public MCIObservableFunctionInterface{
private:
    FeedForwardNeuralNetwork * _ffnn;

public:
    X2(): MCIObservableFunctionInterface(1, 1){}

    virtual void observableFunction(const double * in, double *out){
        out[0] = in[0]*in[0];
    }

};




class NN2Sampling: public MCISamplingFunctionInterface{
private:
    FeedForwardNeuralNetwork * _ffnn;
public:
    NN2Sampling(FeedForwardNeuralNetwork * ffnn): MCISamplingFunctionInterface(1, 1){
        _ffnn = ffnn;
    }

    void samplingFunction(const double *in, double * protovalues){
        _ffnn->setInput(1, in);
        _ffnn->FFPropagate();
        const double v = _ffnn->getOutput(1);
        protovalues[0] = v*v;
    }

    double getAcceptance(){
        return getProtoNew(0)/getProtoOld(0);
    }
};


class NN2Lapl2: public MCIObservableFunctionInterface{
private:
    FeedForwardNeuralNetwork * _ffnn;

public:
    NN2Lapl2(FeedForwardNeuralNetwork * ffnn): MCIObservableFunctionInterface(1, 1){
        _ffnn = ffnn;
    }

    virtual void observableFunction(const double * in, double *out){
        _ffnn->setInput(1, in);
        _ffnn->FFPropagate();
        const double v = _ffnn->getOutput(1);
        const double d2 = _ffnn->getSecondDerivative(1, 0);
        out[0] = v*v*d2;
    }

};


class Lapl2: public MCIObservableFunctionInterface{
private:
    FeedForwardNeuralNetwork * _ffnn;

public:
    Lapl2(FeedForwardNeuralNetwork * ffnn): MCIObservableFunctionInterface(1, 1){
        _ffnn = ffnn;
    }

    virtual void observableFunction(const double * in, double *out){
        _ffnn->setInput(1, in);
        _ffnn->FFPropagate();
        const double d2 = _ffnn->getSecondDerivative(1, 0);
        out[0] = d2;
    }

};




int main(){
    using namespace std;

    const long NMC = 400000;
    const double DX = 0.0001;

    FeedForwardNeuralNetwork * ffnn = new FeedForwardNeuralNetwork(2, 10, 2);
    ffnn->connectFFNN();
    ffnn->addFirstDerivativeSubstrate();
    ffnn->addSecondDerivativeSubstrate();

    MCI * mci = new MCI(1);

    double ** irange = new double*[1];
    irange[0] = new double[2];
    irange[0][0] = -10.;
    irange[0][1] = +10.;

    mci->setIRange(irange);


    // compute   \int_{-10}^{+10} dx nn(x)^2    which is the normalization factor
    NN2 * nn2 = new NN2(ffnn);
    mci->addObservable(nn2);
    double * nn2_av = new double;
    double * nn2_er = new double;
    mci->integrate(NMC, nn2_av, nn2_er);
    cout << "Compute the normalization factor:";
    cout << "    normalization = " << *nn2_av << " +- " << *nn2_er << endl;




    // ---   \int_{-10}^{+10} dx x^2 nn(x)^2
    cout << endl << endl << "We now compute the integral" << endl;
    cout << "    int_{-10}^{+10} dx x^2 nn(x)^2" << endl << endl;

    // compute   \int_{-10}^{+10} dx x^2 nn(x)^2    MC without sampling
    NN2X2 * nn2x2 = new NN2X2(ffnn);
    mci->clearObservables();
    mci->addObservable(nn2x2);
    double * nn2x2_av = new double;
    double * nn2x2_er = new double;
    mci->integrate(NMC, nn2x2_av, nn2x2_er);
    mci->clearObservables();
    cout << "1. MC without sampling: ";
    cout << *nn2x2_av / *nn2_av << " +- " << *nn2x2_er / *nn2_av << endl << endl;


    // compute   \int_{-10}^{+10} dx x^2 nn(x)^2    MC with sampling
    X2 * x2 = new X2();
    mci->addObservable(x2);
    NN2Sampling * nn2_samp = new NN2Sampling(ffnn);
    mci->addSamplingFunction(nn2_samp);
    double * nn2x2_samp_av = new double;
    double * nn2x2_samp_er = new double;
    mci->integrate(NMC, nn2x2_samp_av, nn2x2_samp_er);
    mci->clearObservables();
    mci->clearSamplingFunctions();
    cout << "2. MC with sampling: ";
    cout << *nn2x2_samp_av << " +- " << *nn2x2_samp_er << endl << endl;


    // compute   \int_{-10}^{+10} dx x^2 nn(x)^2    direct integral
    double x = irange[0][0];
    double y = 0.;
    double integral = 0.;
    while (x < irange[0][1]){
        nn2x2->observableFunction(&x, &y);
        integral += y*DX;
        x += DX;
    }
    cout << "3. Direct integral: ";
    cout << "nn2x2_direct = " << integral / *nn2_av << endl << endl;



    // ---   \int_{-10}^{+10} dx \laplacian^2 nn(x)^2
    cout << endl << endl << "We now compute the integral" << endl;
    cout << "    int_{-10}^{+10} dx laplacian^2 nn(x)^2" << endl << endl;

    // compute   \int_{-10}^{+10} dx x^2 nn(x)^2    MC without sampling
    NN2Lapl2 * nn2lapl2 = new NN2Lapl2(ffnn);
    mci->addObservable(nn2lapl2);
    double * nn2lapl2_av = new double;
    double * nn2lapl2_er = new double;
    mci->integrate(NMC, nn2lapl2_av, nn2lapl2_er);
    mci->clearObservables();
    cout << "1. MC without sampling: ";
    cout << *nn2lapl2_av / *nn2_av << " +- " << *nn2lapl2_er / *nn2_av << endl << endl;

    // compute   \int_{-10}^{+10} dx x^2 nn(x)^2    MC with sampling
    Lapl2 * lapl2 = new Lapl2(ffnn);
    mci->addObservable(lapl2);
    mci->addSamplingFunction(nn2_samp);
    double * nn2lapl2_samp_av = new double;
    double * nn2lapl2_samp_er = new double;
    mci->integrate(NMC, nn2lapl2_samp_av, nn2lapl2_samp_er);
    mci->clearObservables();
    mci->clearSamplingFunctions();
    cout << "2. MC with sampling: ";
    cout << *nn2lapl2_samp_av << " +- " << *nn2lapl2_samp_er << endl << endl;

    // compute   \int_{-10}^{+10} dx x^2 nn(x)^2    direct integral
    x = irange[0][0];
    y = 0.;
    integral = 0.;
    while (x < irange[0][1]){
        nn2lapl2->observableFunction(&x, &y);
        integral += y*DX;
        x += DX;
    }
    cout << "3. Direct integral: ";
    cout << integral / *nn2_av << endl;





    delete nn2x2_er;
    delete nn2x2_av;
    delete nn2_er;
    delete nn2_av;
    delete nn2;
    delete[] irange[0];
    delete[] irange;
    delete mci;
    delete ffnn;

    return 0;
}
