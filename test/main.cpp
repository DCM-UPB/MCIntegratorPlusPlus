#include "mci/Estimators.hpp"
#include "mci/MCIntegrator.hpp"
#include <iostream>
#include <random>
#include <math.h>

class Constval: public MCIObservableFunctionInterface
{
public:
    Constval(const int &ndim): MCIObservableFunctionInterface(ndim, 1) {}

    void observableFunction(const double * in, double * out)
    {
        out[0] = 1.3;
    }
};

class Polynom: public MCIObservableFunctionInterface
{
public:
    Polynom(const int &ndim): MCIObservableFunctionInterface(ndim, 1) {}

    void observableFunction(const double * in, double * out)
    {
        out[0]=0.;
        for (int i=0; i<this->getNDim(); ++i)
            {
                out[0] += in[i];
            }
    }
};

class X2: public MCIObservableFunctionInterface
{
public:
    X2(const int &ndim): MCIObservableFunctionInterface(ndim,1) {}

    void observableFunction(const double * in, double * out)
    {
        out[0]=0.;
        for (int i=0; i<this->getNDim(); ++i)
            {
                out[0] += in[i]*in[i];
            }
    }
};

class Gauss: public MCISamplingFunctionInterface
{
public:
    Gauss(const int &ndim): MCISamplingFunctionInterface(ndim,1) {}

    void samplingFunction(const double * in, double * out)
    {
        out[0]=0.;
        for (int i=0; i<this->getNDim(); ++i)
            {
                out[0] += (in[i])*(in[i]);
            }
    }

    double getAcceptance(const double * protoold, const double * protonew)
    {
        return exp(-(protonew[0]-protoold[0]));
    }
};



int main(){
    using namespace std;

    default_random_engine rand_gen;
    uniform_real_distribution<double> rand_num(0.0,1.0);

    long N = 1000l;
    double *x = new double[N];
    for (long i=0; i<N; ++i){
        *(x+i)=3.5 + rand_num(rand_gen);
    }
    double * average = new double(0.);
    double * error = new double(0.);

    mci::UncorrelatedEstimator(N, x, average, error);
    cout << "- UncorrelatedEstimator()" << endl;
    cout << "     average = " << *average << "     error = " << *error << endl << endl;

    int nblocks=12;
    mci::BlockEstimator(N, x, nblocks, average, error);
    cout << "- BlockEstimator()" << endl;
    cout << "     average = " << *average << "     error = " << *error << endl << endl;

    mci::CorrelatedEstimator(N, x, average, error);
    cout << "- CorrelatedEstimator()" << endl;
    cout << "     average = " << *average << "     error = " << *error << endl << endl;

    delete average;
    delete error;
    delete[] x;



    cout << endl << "Multidimensional version of Estimators" << endl << endl;


    int nd=2;
    double **data = new double*[N];
    for (int i=0; i<N; ++i){ *(data+i) = new double[nd]; }

    for (long i=0; i<N; ++i){
        *(*(data+i)) = 2.5 + rand_num(rand_gen);
        *(*(data+i)+1) = -5.5 + rand_num(rand_gen);
    }
    average = new double[nd];
    error = new double[nd];

    mci::MultiDimUncorrelatedEstimator(N, nd, data, average, error);
    cout << "- UncorrelatedEstimator()" << endl;
    cout << "     average1 = " << *average << "     error1 = " << *error << endl;
    cout << "     average2 = " << *(average+1) << "     error2 = " << *(error+1) << endl << endl;

    mci::MultiDimBlockEstimator(N, nd, data, nblocks, average, error);
    cout << "- MultiDimBlockEstimator()" << endl;
    cout << "     average1 = " << *average << "     error1 = " << *error << endl;
    cout << "     average2 = " << *(average+1) << "     error2 = " << *(error+1) << endl << endl;

    mci::MultiDimCorrelatedEstimator(N, nd, data, average, error);
    cout << "- MultiDimCorrelatedEstimator()" << endl;
    cout << "     average1 = " << *average << "     error1 = " << *error << endl;
    cout << "     average2 = " << *(average+1) << "     error2 = " << *(error+1) << endl << endl;

    for (int j=0; j<N; ++j){ delete [] *(data+j); }
    delete[] data;
    delete[] average;
    delete[] error;



    cout << endl << "Monte Carlo integrator" << endl << endl;


    nd = 3;
    MCI mci(nd);
    cout << "Initialized a MCI object for an integration in a space with ndim=" << mci.getNDim() << endl;

    double **irange = new double*[nd];
    for (int i=0; i<nd; ++i){
        irange[i] = new double[2];
        irange[i][0] = 0.;
        irange[i][1] = 1.;
    }
    mci.setIRange(irange);
    cout << "Integration range: " << mci.getIRange(0,0) << "   " << mci.getIRange(0,1) << endl;
    for (int i=1; i<nd; ++i){
        cout << "                   " << mci.getIRange(i,0) << "   " << mci.getIRange(i,1) << endl;
    }

    double *r=new double[nd];
    for (int i=0; i<nd; ++i){ r[i]=0.5;}//-1.*i; }
    mci.setX(r);
    cout << "Set starting position X: ";
    for (int i=0; i<nd; ++i){ cout << mci.getX(i) << "   "; }
    cout << endl;

    Constval constval(nd);
    Polynom polynom(nd);
    mci.addObservable(&constval);
    mci.addObservable(&polynom);

    N=10000l;
    average = new double[2];
    error = new double[2];
    mci.integrate(N,average,error);

    cout << "Compute average: " << endl;
    cout << "Average 1 (Constval=1.3)         = " << average[0] << " +- " << error[0] << endl;
    cout << "Average 2 (Polynom=x+y+z -> 1.5) = " << average[1] << " +- " << error[1] << endl << endl << endl;

    delete[] r;
    for (int i=0; i<nd; ++i){delete [] irange[i];}
    delete [] irange;
    delete [] average;
    delete [] error;

    nd = 1;
    MCI mci_1d(1);
    irange = new double*[nd];
    for (int i=0; i<nd; ++i){
        irange[i] = new double[2];
        irange[i][0] = -1.;
        irange[i][1] = 1.;
    }
    mci_1d.setIRange(irange);
    X2 x2(nd);
    mci_1d.addObservable(&x2);
    average = new double[1];
    error = new double[1];
    N=10000;
    mci_1d.integrate(N,average,error);
    cout << "Integral of x^2 between -1 and +1 (expected 2./3.): " << endl;
    cout << "Integral = " << *average << " +- " << *error << endl << endl << endl;
    delete[] average;
    delete[] error;
    for (int i=0; i<nd; ++i){delete [] irange[i];}
    delete [] irange;


    nd = 1;
    MCI mci_1dgauss(1);
    Gauss gauss(nd);
    mci_1dgauss.addSamplingFunction(&gauss);
    mci_1dgauss.addObservable(&x2);
    mci_1dgauss.addObservable(&constval);
    mci_1dgauss.storeObservablesOnFile("observables.txt", 100);
    mci_1dgauss.storeWalkerPositionsOnFile("walker.txt", 100);
    average = new double[2];
    error = new double[2];
    N=10000;
    mci_1dgauss.integrate(N,average,error);
    cout << "Integral of x^2 between -1 and +1 sampling from a normalized gaussian (expected 1./2.): " << endl;
    cout << "Integral = " << *average << " +- " << *error << endl << endl << endl;
    delete[] average;
    delete[] error;

    //nd = 1;
    //MCI mci_1dgauss_complex(1);
    //mci_1dgauss_complex.setProtoSamplingFunction(&expgauss1, 2);
    //mci_1dgauss_complex.setProtoSamplingFunction(&expgauss2, 1);
    //mci_1dgauss_complex.setAcceptanceFunction(&gaussratio);
    //mci_1dgauss_complex.setObservable(&x2,1);
    //mci_1dgauss_complex.setObservable(&plain_x,1);
    //average = new double[2];
    //error = new double[2];
    //N=100000;
    //mci_1dgauss_complex.storeObservablesOnFile("observables.txt", 100);
    //mci_1dgauss_complex.storeWalkerPositionsOnFile("walker.txt", 100);
    //mci_1dgauss_complex.integrate(N,average,error);
    //cout << "Integral of x^2 between -1 and +1 sampling from a normalized gaussian (expected 1./2.): " << endl;
    //cout << "Integral = " << average[0] << " +- " << error[0] << endl ;
    //cout << "Integral of x between -1 and +1 sampling from a normalized gaussian (expected 0.): " << endl;
    //cout << "Integral = " << average[1] << " +- " << error[1] << endl ;
    //cout << " + Acceptance = " << mci_1dgauss_complex.getAcceptanceRate() << endl << endl << endl;
    //delete[] average;
    //delete[] error;

}
