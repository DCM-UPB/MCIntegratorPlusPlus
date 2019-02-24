#include "mci/Estimators.hpp"
#include "mci/MCIntegrator.hpp"
#include <cmath>
#include <iostream>
#include <random>

class Constval: public MCIObservableFunctionInterface
{
public:
    explicit Constval(const int &ndim): MCIObservableFunctionInterface(ndim, 1) {}

    void observableFunction(const double *  /*in*/, double * out) override
    {
        out[0] = 1.3;
    }
};

class Polynom: public MCIObservableFunctionInterface
{
public:
    explicit Polynom(const int &ndim): MCIObservableFunctionInterface(ndim, 1) {}

    void observableFunction(const double * in, double * out) override
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
    explicit X2(const int &ndim): MCIObservableFunctionInterface(ndim,1) {}

    void observableFunction(const double * in, double * out) override
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
    explicit Gauss(const int &ndim): MCISamplingFunctionInterface(ndim,1) {}

    void samplingFunction(const double * in, double * out) override
    {
        out[0]=0.;
        for (int i=0; i<this->getNDim(); ++i)
            {
                out[0] += (in[i])*(in[i]);
            }
    }

    double getAcceptance(const double * protoold, const double * protonew) override
    {
        return exp(-(protonew[0]-protoold[0]));
    }
};



int main(){
    using namespace std;

    default_random_engine rand_gen;
    uniform_real_distribution<double> rand_num(0.0,1.0);

    long N = 1000l;
    auto *x = new double[N];
    for (long i=0; i<N; ++i){
        *(x+i)=3.5 + rand_num(rand_gen);
    }
    double avg1D = 0.;
    double err1D = 0.;

    mci::UncorrelatedEstimator(N, x, &avg1D, &err1D);
    cout << "- UncorrelatedEstimator()" << endl;
    cout << "     avg1D = " << avg1D << "     error = " << err1D << endl << endl;

    int nblocks=12;
    mci::BlockEstimator(N, x, nblocks, &avg1D, &err1D);
    cout << "- BlockEstimator()" << endl;
    cout << "     average = " << avg1D << "     error = " << err1D << endl << endl;

    mci::CorrelatedEstimator(N, x, &avg1D, &err1D);
    cout << "- CorrelatedEstimator()" << endl;
    cout << "     average = " << avg1D << "     error = " << err1D << endl << endl;

    delete[] x;



    cout << endl << "Multidimensional version of Estimators" << endl << endl;


    int nd=2;
    auto **data = new double*[N];
    for (int i=0; i<N; ++i){ *(data+i) = new double[nd]; }

    for (long i=0; i<N; ++i){
        data[i][0] = 2.5 + rand_num(rand_gen);
        data[i][1] = -5.5 + rand_num(rand_gen);
    }
    double avg2D[nd];
    double err2D[nd];

    mci::MultiDimUncorrelatedEstimator(N, nd, data, avg2D, err2D);
    cout << "- UncorrelatedEstimator()" << endl;
    cout << "     average1 = " << avg2D[0] << "     error1 = " << err2D[0] << endl;
    cout << "     average2 = " << avg2D[1] << "     error2 = " << err2D[1] << endl << endl;

    mci::MultiDimBlockEstimator(N, nd, data, nblocks, avg2D, err2D);
    cout << "- MultiDimBlockEstimator()" << endl;
    cout << "     average1 = " << avg2D[0] << "     error1 = " << err2D[0] << endl;
    cout << "     average2 = " << avg2D[1] << "     error2 = " << err2D[1] << endl << endl;

    mci::MultiDimCorrelatedEstimator(N, nd, data, avg2D, err2D);
    cout << "- MultiDimCorrelatedEstimator()" << endl;
    cout << "     average1 = " << avg2D[0] << "     error1 = " << err2D[0] << endl;
    cout << "     average2 = " << avg2D[1] << "     error2 = " << err2D[1] << endl << endl;

    for (int j=0; j<N; ++j){ delete [] *(data+j); }
    delete[] data;

    cout << endl << "Monte Carlo integrator" << endl << endl;


    nd = 3;
    MCI mci(nd);
    cout << "Initialized a MCI object for an integration in a space with ndim=" << mci.getNDim() << endl;

    auto **irange = new double*[nd];
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

    double r[nd];
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
    mci.integrate(N,avg2D,err2D);

    cout << "Compute average: " << endl;
    cout << "Average 1 (Constval=1.3)         = " << avg2D[0] << " +- " << err2D[0] << endl;
    cout << "Average 2 (Polynom=x+y+z -> 1.5) = " << avg2D[1] << " +- " << err2D[1] << endl << endl << endl;

    for (int i=0; i<nd; ++i){delete [] irange[i];}
    delete [] irange;

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
    N=10000;
    mci_1d.integrate(N,&avg1D,&err1D);
    cout << "Integral of x^2 between -1 and +1 (expected 2./3.): " << endl;
    cout << "Integral = " << avg1D << " +- " << err1D << endl << endl << endl;
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
    N=10000;
    mci_1dgauss.integrate(N,avg2D,err2D);
    cout << "Integral of x^2 between -1 and +1 sampling from a normalized gaussian (expected 1./2.): " << endl;
    cout << "Integral = " << avg2D[0] << " +- " << err2D[0] << endl << endl << endl;
}
