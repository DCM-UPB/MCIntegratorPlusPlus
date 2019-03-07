#ifndef MCI_CALLBACKONACCEPTANCEINTERFACE_HPP
#define MCI_CALLBACKONACCEPTANCEINTERFACE_HPP

namespace mci
{
    class CallBackOnAcceptanceInterface
    {
    protected:
        const int _ndim; //dimension of the input array (walker position)

    public:
        explicit CallBackOnAcceptanceInterface(int ndim): _ndim(ndim) {}
        virtual ~CallBackOnAcceptanceInterface()= default;

        // Getters
        int getNDim(){return _ndim;}


        // --- METHODS THAT MUST BE IMPLEMENTED

        // Call-back function, called if a move is accepted by the MCI Integrator
        virtual void callBackFunction(const double in[], bool flag_observation) = 0;
        //                                   ^walker position
    };
}

#endif
