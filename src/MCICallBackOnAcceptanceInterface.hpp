#ifndef CALL_BACK_ON_ACCEPTANCE_INTERFACE
#define CALL_BACK_ON_ACCEPTANCE_INTERFACE


class MCICallBackOnAcceptanceInterface
{
protected:
    int _ndim; //dimension of the input array (walker position)

public:
    MCICallBackOnAcceptanceInterface(const int &ndim){
        _ndim=ndim;
    }
    virtual ~MCICallBackOnAcceptanceInterface(){}

    // Getters
    int getNDim(){ return _ndim;}


    // --- METHODS THAT MUST BE IMPLEMENTED

    // Call-back function, called if a move is accepted by the MCI Integrator
    virtual void callBackFunction(const double *in) = 0;
    //                                   ^walker position

};


#endif
