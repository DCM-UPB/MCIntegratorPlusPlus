#ifndef MCI_CALLBACKONMOVEINTERFACE_HPP
#define MCI_CALLBACKONMOVEINTERFACE_HPP

namespace mci
{
    // Base class for adding callbacks to MCI. When added to MCI, they will be invoked on every walker move.
    //
    // Derive from this and implement the virtual callBackFunction(...) method.
    // You also need to provide a protected _clone method returning a raw pointer of
    // type CallBackOnMoveInterface, pointing to an object of your type MyCallback, e.g.:
    //
    // class MyCallback: public CallBackOnMoveInterface {
    // protected:
    //     CallBackOnMoveInterface * _clone() const override {
    //         return new MyCallback(...); // create a cloned version here
    //     }
    // public:
    //     void samplingFunction(...) overwrite;
    //     double getAcceptance(...) const overwrite;
    //     ...
    // };
    //
    // Your class will have a public clone() method returning std::unique_ptr<CallBackOnMoveInterface> .
    // If you want/need it, also create a non-overriding clone() method returning a pointer of type MyCallback.
    class CallBackOnMoveInterface: public Clonable<CallBackOnMoveInterface>
    {
    protected:
        const int _ndim; //dimension of the input array (walker position)

    public:
        explicit CallBackOnMoveInterface(int ndim): _ndim(ndim) {}
        virtual ~CallBackOnMoveInterface()= default;

        // Getters
        int getNDim() const {return _ndim;}


        // --- METHODS THAT MUST BE IMPLEMENTED

        // Call-back function, called if a move is accepted by the MCI Integrator
        virtual void callBackFunction(const double in[], bool flag_observation) = 0;
        //                                   ^walker position
    };
}  // namespace mci

#endif
