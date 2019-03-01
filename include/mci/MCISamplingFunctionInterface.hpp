#ifndef MCI_MCISAMPLINGFUNCTIONINTERFACE_HPP
#define MCI_MCISAMPLINGFUNCTIONINTERFACE_HPP

#include <algorithm>

class MCISamplingFunctionInterface
{
protected:
    const int _ndim; //dimension of the input array (walker position)
    int _nproto{}; //number of proto sampling functions given as output
    double * _protonew; //array containing the new proto sampling functions
    double * _protoold; //array containing the old proto sampling functions

public:
    MCISamplingFunctionInterface(const int &ndim, const int &nproto): _ndim(ndim)
    {
        _protonew = nullptr;
        _protoold = nullptr;
        setNProto(nproto);
    }
    virtual ~MCISamplingFunctionInterface()
    {
        delete[] _protoold;
        delete[] _protonew;

    }


    // Setters
    void setNProto(const int &nproto){
        _nproto=nproto;
         delete[] _protoold;
         delete[] _protonew;
        _protonew = new double[_nproto];
        _protoold = new double[_nproto];
        std::fill(_protonew, _protonew+_nproto, 0.);
        std::fill(_protoold, _protoold+_nproto, 0.);
    }


    // Getters
    int getNDim(){ return _ndim;}
    int getNProto(){ return _nproto;}

    // Utilities
    void newToOld()
    {   // pointer swap
        double * const foo = _protonew;
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
