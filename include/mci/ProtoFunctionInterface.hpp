#ifndef MCI_PROTOFUNCTIONINTERFACE_HPP
#define MCI_PROTOFUNCTIONINTERFACE_HPP

namespace mci
{
    // Base class for all proto functions
    //
    // Proto Functions are functions to be evaluated for every new proposed trial move within
    // a accept/reject sampling scheme. In that scenario, they can often make use of temporary
    // values that we call proto values, to optimize updates for single/few-particle moves. A
    // proto function stores a "constant" set of temporaries (_protoold) for the previously
    // accepted position and a variable set (_protonew) for the newly proposed position.
    // Whenever the new position is accepted (externally), the method newToOld() get's called to
    // copy new values to old.
    //
    // Child classes have to implement the virtual protoFunction(...) method, which takes a set
    // of inputs and calculates the proto values. The derived interfaces decide how to make actual
    // use of the values and also how to handle selective updating.
    // If you have own proto-value like data and for some reason don't want to store them in the
    // protovalue arrays, please implement the protected method _newToOld and copy your new data
    // to your old data. This makes sure the old values are initialized at the first step and
    // copied on newToOld. Update your custom data in protoFunction and in the derived interface's
    // selective updating methods.
    class ProtoFunctionInterface
    {
    protected:
        const int _ndim; // dimension of the input array (walker position)
        int _nproto; // number of proto values calculated in protoFunction
        double * _protoold; // ptr to block containing the new old values
        double * _protonew; // ptr to block containing the new proto values

        // internal setters
        void setNProto(int nproto); // you may freely choose the amount of values you need

        // Overwrite this if you have own data to copy on acceptance/rejection.
        // It will be called in the public newToOld()/oldToNew() methods.
        virtual void _newToOld() {};
        virtual void _oldToNew() {};

    public:
        ProtoFunctionInterface(int ndim, int nproto);
        virtual ~ProtoFunctionInterface();

        // Getters
        int getNDim() const { return _ndim;}
        int getNProto() const { return _nproto;}

        // --- Main operational methods

        // initializer for proto values
        void initializeProtoValues(const double in[]);

        // copy new to old protov, call _newToOld()/_oldToNew()
        void newToOld(); // called on acceptance
        void oldToNew(); // called on rejection

        // --- METHOD THAT MUST BE IMPLEMENTED

        // Function that MCI uses to calculate your proto-function values
        // Calculate and store them in the protovalues array as they are
        // expected from the derived interface's methods.
        virtual void protoFunction(const double in[], double protovalues[]) = 0; // e.g. the summands of an exponent: exp(-sum(protovalues))
        //                             ^walker position  ^resulting proto-values
    };

}  // namespace mci

#endif
