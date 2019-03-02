#ifndef MCI_MCICALLBACKONACCEPTANCEINTERFACE_HPP
#define MCI_MCICALLBACKONACCEPTANCEINTERFACE_HPP


class MCICallBackOnAcceptanceInterface
{
protected:
    const int _ndim; //dimension of the input array (walker position)

public:
    explicit MCICallBackOnAcceptanceInterface(int ndim): _ndim(ndim) {}
    virtual ~MCICallBackOnAcceptanceInterface()= default;

    // Getters
    int getNDim(){return _ndim;}


    // --- METHODS THAT MUST BE IMPLEMENTED

    // Call-back function, called if a move is accepted by the MCI Integrator
    virtual void callBackFunction(const double *in, bool flag_observation) = 0;
    //                                   ^walker position
};


#endif
