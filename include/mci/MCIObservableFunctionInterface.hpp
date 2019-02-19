#ifndef MCI_OBSERVABLE_FUNCTION_INTERFACE
#define MCI_OBSERVABLE_FUNCTION_INTERFACE



class MCIObservableFunctionInterface
{
protected:
    int _ndim;  //dimension of the input array (walker poistion)
    int _nobs;  //number of values provided by the observable
    double * _obs; //array that stores the last observables computed


public:
    MCIObservableFunctionInterface(const int &ndim, const int &nobs)
    {
        _ndim=ndim; _nobs=nobs;
        _obs = new double[nobs];
        for (int i=0; i<nobs; ++i){_obs[i]=0.;}
    }
    virtual ~MCIObservableFunctionInterface()
    {
        delete[] _obs;
    }

    int getNObs(){return _nobs;}
    int getNDim(){return _ndim;}
    double getObservable(const int &i){ return _obs[i]; }

    void computeObservables(const double *in)
    {
        observableFunction(in,_obs);
    }



    // --- METHOD THAT MUST BE IMPLEMENTED

    // Compute the observable and store it inside out
    virtual void observableFunction(const double * in, double * out) = 0;
    //                               ^input = walker positions  ^resulting observables
};



#endif
