#ifndef MCI_SAMPLING_FUNCTION_INTERFACE
#define MCI_SAMPLING_FUNCTION_INTERFACE


class MCISamplingFunctionInterface
{
protected:
    int _ndim; //dimension of the input array (walker position)
    int _nproto; //number of proto sampling functions given as output
    double * _protonew; //array containing the new proto sampling functions
    double * _protoold; //array containing the old proto sampling functions

public:
    MCISamplingFunctionInterface(const int &ndim, const int &nproto)
    {
        _ndim=ndim;
        _protonew = 0; _protoold = 0;
        setNProto(nproto);
    }
    virtual ~MCISamplingFunctionInterface()
    {
        delete[] _protonew; delete[] _protoold;
    }


    // Setters
    void setNProto(const int &nproto){
        _nproto=nproto;
        if (_protonew != 0) delete[] _protonew;
        if (_protoold != 0) delete[] _protoold;
        _protonew = new double[_nproto]; _protoold = new double[_nproto];
        for (int i=0; i<_nproto; ++i){ _protonew[i]=0.; }
        for (int i=0; i<_nproto; ++i){ _protoold[i]=0.; }
    }


    // Getters
    int getNDim(){ return _ndim;}
    int getNProto(){ return _nproto;}

    // Utilities
    void newToOld()
    {   // pointer swap
        double * foo=_protonew;
        _protonew=_protoold;
        _protoold=foo;
    }

    void computeNewSamplingFunction(const double * in)
    {
        samplingFunction(in, _protonew);
    }

    double getAcceptance(){
        return getAcceptance(_protoold, _protonew);
    }


    // --- METHODS THAT MUST BE IMPLEMENTED
    // Function that MCI uses for the proto-sampling function. Computes _protonew
    virtual void samplingFunction(const double *in, double * protovalues) = 0;
    //                                      ^walker position  ^resulting proto-values

    // Acceptance function, that uses the new and old values of the proto sampling function
    virtual double getAcceptance(const double * protoold, const double * protonew) = 0;
};


#endif
