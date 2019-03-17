#ifndef MCI_CALLBACKONMOVEINTERFACE_HPP
#define MCI_CALLBACKONMOVEINTERFACE_HPP

#include "mci/Clonable.hpp"
#include "mci/WalkerState.hpp"

namespace mci
{
    // Base class for adding callbacks to MCI. When added to MCI, they will be invoked directly
    // after every walker move & acceptance/rejection. You will get passed the WalkerState struct,
    // containing old positions, new positions (which might be rejected), which and how many indices
    // are changed, and finally if the step was accepted.
    //
    // Derive from this and implement the virtual callBackFunction(...) method.
    // You also need to provide a protected _clone method returning a raw pointer of
    // type CallBackOnMoveInterface, pointing to an object of your type MyCallback, e.g.:
    //
    // class MyCallback: public CallBackOnMoveInterface {
    // protected:
    //     CallBackOnMoveInterface * _clone() const final {
    //         return new MyCallback(...); // create a cloned version here
    //     }
    // public:
    //     void callBackFunction(...) final;
    //     ...
    // };
    //
    // Your class will have a public clone() method returning std::unique_ptr<CallBackOnMoveInterface> .
    class CallBackOnMoveInterface: public Clonable<CallBackOnMoveInterface>
    {
    protected:
        const int _ndim; //dimension of the input array (walker position)

    public:
        explicit CallBackOnMoveInterface(int ndim): _ndim(ndim) {}

        // Getters
        int getNDim() const {return _ndim;}


        // --- METHODS THAT MUST BE IMPLEMENTED

        // Call-back function, called by the MCI Integrator at every step
        virtual void callBackFunction(const WalkerState &wlkstate) = 0;
        //                                              ^walker state struct
    };
}  // namespace mci

#endif
