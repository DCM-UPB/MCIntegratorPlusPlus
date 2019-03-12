#ifndef MCI_SAMPLINGFUNCTIONINTERFACE_HPP
#define MCI_SAMPLINGFUNCTIONINTERFACE_HPP

#include "mci/Clonable.hpp"

namespace mci
{
    // Base class for MC sampling functions (probability distribution functions)
    //
    // Derive from this and implement the virtual samplingFunction(...) and getAcceptance (..) methods.
    // You also need to provide a protected _clone method returning a raw pointer of
    // type SamplingFunctionInterface, pointing to an object of your type MyPDF, e.g.:
    //
    // class MyPDF: public SamplingFunctionInterface {
    // protected:
    //     SamplingFunctionInterface * _clone() const override {
    //         return new MyPDF(...); // create a cloned version here
    //     }
    // public:
    //     void samplingFunction(...) overwrite;
    //     double getAcceptance(...) const overwrite;
    //     ...
    // };
    //
    // Your class will have a public clone() method returning std::unique_ptr<SamplingFunctionInterface> .
    // If you want/need it, also create a non-overriding clone() method returning a pointer of type MyPDF.
    class SamplingFunctionInterface: public Clonable<SamplingFunctionInterface>
    {
    protected:
        const int _ndim; //dimension of the input array (walker position)
        int _nproto; //number of proto sampling functions given as output
        double * _protonew; //array containing the new proto sampling functions
        double * _protoold; //array containing the old proto sampling functions

        // Setters
        void setNProto(int nproto); // you may freely choose the amount of values you need

    public:
        SamplingFunctionInterface(int ndim, int nproto);

        virtual ~SamplingFunctionInterface();

        // Getters
        int getNDim() const { return _ndim;}
        int getNProto() const { return _nproto;}

        // Main operational methods
        void computeNewSamplingFunction(const double in[]) { samplingFunction(in, _protonew); }
        void computeNewSamplingFunction(const double xold[], const double xnew[], int nchanged, const int changedIdx[]) { samplingFunction(xold, xnew, nchanged, changedIdx, _protonew); }
        double getAcceptance() const { return getAcceptance(_protoold, _protonew); }

        // overwrite this if you have own data to swap, but remember
        // to call SamplingFunctionInterface::newToTold() in your newToOld()
        virtual void newToOld(); // swap old and new protovalues


        // --- METHODS THAT MUST BE IMPLEMENTED
        // Function that MCI uses for the proto-sampling function. Computes _protonew
        virtual void samplingFunction(const double in[], double protovalues[]) = 0;
        //                                      ^walker position  ^resulting proto-values

        // Acceptance function, that uses the new and old values of the proto sampling function
        virtual double getAcceptance(const double protoold[], const double protonew[]) const = 0;

        // --- OPTIONALLY ALSO OVERWRITE THIS (to optimize for single/few particle moves)
        // Compute the proto values, given both previous and current walker positions (xold/xnew), and additionally the
        // array changedIdx containing the indices of the nchanged(!) elements that differ between xold and xnew. The
        // indices in changedIdx are guaranteed to be in ascending order.
        // This means:
        //     a) you never have to store the previous walker position in your child class and
        //     b) you can use the indices in changedIdx to provide efficient recalculation of your protovalues
        // If full recalculation is almost always more efficient in your case, you may also choose not to overwrite this method.
        virtual void samplingFunction(const double xold[], const double xnew[], int /*nchanged*/, const int[] /*changedIdx*/, double protov[]) {
            samplingFunction(xnew, protov); // default to "calculate all"
        }
    };
}  // namespace mci

#endif
