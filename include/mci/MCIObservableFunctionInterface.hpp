#ifndef MCI_MCIOBSERVABLEFUNCTIONINTERFACE_HPP
#define MCI_MCIOBSERVABLEFUNCTIONINTERFACE_HPP

#include <algorithm>

class MCIObservableFunctionInterface
{
protected:
    const int _ndim;  //dimension of the input array (walker poistion)
    const int _nobs;  //number of values provided by the observable
    double * _obs; //array that stores the last observables computed


public:
    MCIObservableFunctionInterface(int ndim, int nobs): _ndim(ndim), _nobs(nobs)
    {
        _obs = new double[nobs];
        std::fill(_obs, _obs+nobs, 0.);
    }
    virtual ~MCIObservableFunctionInterface()
    {
        delete[] _obs;
    }

    int getNObs(){return _nobs;}
    int getNDim(){return _ndim;}
    double getValue(int i){ return _obs[i]; }
    const double * getValues(){ return _obs; }

    void computeValues(const double in[])
    {
        observableFunction(in,_obs);
    }



    // --- METHOD THAT MUST BE IMPLEMENTED

    // Compute the observable and store it inside out
    virtual void observableFunction(const double in[], double out[]) = 0;
    //                               ^input = walker positions  ^resulting observables
};



#endif
