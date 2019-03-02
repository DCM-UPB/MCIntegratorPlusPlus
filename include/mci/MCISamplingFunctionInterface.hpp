#ifndef MCI_MCISAMPLINGFUNCTIONINTERFACE_HPP
#define MCI_MCISAMPLINGFUNCTIONINTERFACE_HPP

class MCISamplingFunctionInterface
{
protected:
    const int _ndim; //dimension of the input array (walker position)
    int _nproto; //number of proto sampling functions given as output
    double * _protonew; //array containing the new proto sampling functions
    double * _protoold; //array containing the old proto sampling functions

public:
    MCISamplingFunctionInterface(int ndim, int nproto);

    virtual ~MCISamplingFunctionInterface();

    // Setters
    void setNProto(int nproto);

    // Getters
    int getNDim(){ return _ndim;}
    int getNProto(){ return _nproto;}

    // Utilities
    void newToOld(); // swap old and new protovalues

    void computeNewSamplingFunction(const double * in) { samplingFunction(in, _protonew); }

    double getAcceptance() { return getAcceptance(_protoold, _protonew); }


    // --- METHODS THAT MUST BE IMPLEMENTED
    // Function that MCI uses for the proto-sampling function. Computes _protonew
    virtual void samplingFunction(const double *in, double * protovalues) = 0;
    //                                      ^walker position  ^resulting proto-values

    // Acceptance function, that uses the new and old values of the proto sampling function
    virtual double getAcceptance(const double * protoold, const double * protonew) = 0;
};


#endif
