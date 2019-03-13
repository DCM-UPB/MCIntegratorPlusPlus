#ifndef MCI_SAMPLINGFUNCTIONINTERFACE_HPP
#define MCI_SAMPLINGFUNCTIONINTERFACE_HPP

#include "mci/Clonable.hpp"

namespace mci
{
    // Base class for MC sampling functions (probability distribution functions)
    //
    // Derive from this and implement the virtual protoFunction(...) and acceptanceFunction(..)
    // methods. You also need to provide a protected _clone method returning a raw pointer of
    // type SamplingFunctionInterface, pointing to an object of your type MyPDF, e.g.:
    //
    // class MyPDF: public SamplingFunctionInterface {
    // protected:
    //     SamplingFunctionInterface * _clone() const override {
    //         return new MyPDF(...); // create a cloned version here
    //     }
    // public:
    //     void protoFunction(...) overwrite;
    //     double acceptanceFunction(...) const overwrite;
    //     ...
    // };
    //
    // Your class will have a public clone() method returning std::unique_ptr<SamplingFunctionInterface> .
    // If you want/need it, also create a non-overriding clone() method returning a pointer of type MyPDF.
    class SamplingFunctionInterface: public Clonable<SamplingFunctionInterface>
    {
    protected:
        const int _ndim; // dimension of the input array (walker position)
        int _nproto; // number of proto sampling functions given as output
        double * _protonew; // array containing the new proto values (temporaries to compute scalar function value)
        double * _protoold; // array containing the new old values

        // internal setters
        void setNProto(int nproto); // you may freely choose the amount of values you need

        // Overwrite this if you have own data to copy on acceptance
        // acceptance. It will be called in base's public newToOld()
        virtual void _newToOld() {};

    public:
        SamplingFunctionInterface(int ndim, int nproto);
        virtual ~SamplingFunctionInterface();

        // Getters
        int getNDim() const { return _ndim;}
        int getNProto() const { return _nproto;}


        // --- Main operational methods

        // initializer for old
        void computeOldProtoValues(const double in[]) { protoFunction(in, _protoold); }

        // compute full protonew and return acceptance
        double computeAcceptance(const double in[]);

        // update protonew and return acceptance, given the nchanged indices changedIdx, that differ between xold and xnew
        double computeAcceptance(const double xold[], const double xnew[], int nchanged, const int changedIdx[]);

        // copy new to old protov (we need to copy, not swap, to allow elementary updates), call _newToOld()
        void newToOld();


        // --- METHODS THAT MUST BE IMPLEMENTED

        // Function that MCI uses to calculate your proto-sampling function values.
        // Calculate them as they are expected from your acceptanceFunction().
        virtual void protoFunction(const double in[], double protovalues[]) = 0;
        //                             ^walker position  ^resulting proto-values

        // Acceptance function, that uses the old and new proto sampling function values
        // If your actual sampling function is an exponential like exp(-sum(pv)), you would compute
        // something like exp( -sum(protonew)+sum(protoold) ) here, i.e. division of exponentials.
        virtual double acceptanceFunction(const double protoold[], const double protonew[]) const = 0;


        // --- OPTIONALLY ALSO OVERWRITE THIS (to optimize for single/few particle moves)
        // Return step acceptance AND update(!) protonew elements, given both previous and current
        // walker positions (xold/xnew), and additionally the array changedIdx containing the indices
        // of the nchanged elements that differ between xold and xnew. The indices in changedIdx are
        // guaranteed to be in ascending order.
        // This means:
        //     a) you never have to store the previous walker position in your child class and
        //     b) you can use the indices in changedIdx to provide efficient recalculation of your protovalues
        // Remember that in this method you should only update the protov[] elements that need to change due
        // to the nchanged input indices in changedIdx.
        // If full recalculation is more efficient in your case, you may also choose not to overwrite this method.
        virtual double updatedAcceptance(const double[] /*xold*/, const double xnew[],
                                           int /*nchanged*/, const int[] /*changedIdx*/,
                                           const double protoold[], double protonew[] /* update this! */)
        {
            // default to "calculate all"
            this->protoFunction(xnew, protonew);
            return this->acceptanceFunction(protoold, protonew);
        }
    };

}  // namespace mci

#endif
