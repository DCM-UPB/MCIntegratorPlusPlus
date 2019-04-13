#ifndef MCI_SAMPLINGFUNCTIONINTERFACE_HPP
#define MCI_SAMPLINGFUNCTIONINTERFACE_HPP

#include "mci/Clonable.hpp"
#include "mci/ProtoFunctionInterface.hpp"
#include "mci/WalkerState.hpp"

namespace mci
{
// Base class for MC sampling functions (probability distribution functions)
//
// Derive from this and implement the protoFunction(...) (from ProtoFunctionInterface)
// and acceptanceFunction(..) methods. You also need to provide a protected _clone method
// returning a raw pointer of type SamplingFunctionInterface, pointing to an object of
// your type MyPDF, e.g.:
//
// class MyPDF: public SamplingFunctionInterface {
// protected:
//     SamplingFunctionInterface * _clone() const final {
//         return new MyPDF(...); // create a cloned version here
//     }
// public:
//     void protoFunction(...) final;
//     double samplingFunction(...) const final;
//     double acceptanceFunction(...) const final;
//     ...
// };
//
// Now your class has a public clone() method returning std::unique_ptr<SamplingFunctionInterface> .
//
// IMPORTANT: Formally, we require normalized PDFs, i.e. the integral of the PDF over the integration
// domain should be 1. However, in our program any (normalizable) PDF behaves as if it was normalized.
// This is usually very desirable behavior, but may lead to confusion when this "auto-normalization" is
// not expected or desired. So remember: We assume the PDF is either normalized by you or you want the
// integral with the normalized version anyway.
class SamplingFunctionInterface: public ProtoFunctionInterface, public Clonable<SamplingFunctionInterface>
{
protected:
    SamplingFunctionInterface(int ndim, int nproto): ProtoFunctionInterface(ndim, nproto) {}

public:
    // --- Main operational methods

    // return value of old sampling function
    double getOldSamplingFunction() const { return this->samplingFunction(_protoold); }

    // update protonew and return acceptance, given the Walkerstate, which
    // contains the changed indices changedIdx, that differ between xold and xnew
    double computeAcceptance(const WalkerState &wlk)
    {
        if (wlk.nchanged < _ndim) {
            return this->updatedAcceptance(wlk, _protoold, _protonew);
        }
        // all elements have changed
        this->protoFunction(wlk.xnew, _protonew);
        return this->acceptanceFunction(_protoold, _protonew);
    }


    // --- METHODS THAT MUST BE IMPLEMENTED

    // Remember to implement the protoFunction from ProtoFunctionInterface!
    // In the "e.g." comments below it is assumed that your sampling function is of the form
    // exp(-sum(protovalues)), where protovalues are the values computed in your protoFunction().

    // Function that can be used to calculate the true value of your sampling function,
    // from the given set of proto values. Note that this function is not used by the main
    // parts of MCI, but for example certain trial moves might require this method.
    virtual double samplingFunction(const double protovalues[]) const = 0; // e.g. exp(-sum(protovalues))

    // Acceptance function, that uses the old and new proto sampling function values to compute
    // the acceptance quotient to use in the metropolis criterion.
    virtual double acceptanceFunction(const double protoold[], const double protonew[]) const = 0; // e.g. exp(-sum(protonew)+sum(protoold))


    // --- OPTIONALLY ALSO OVERWRITE THIS (to optimize for single/few particle moves)
    // Return step acceptance AND update(!) protonew elements, given the WalkerSate, which
    // contains previous and current walker positions (xold/xnew), and additionally the array
    // changedIdx containing the indices of the nchanged elements that differ between xold and xnew.
    // The indices in changedIdx are guaranteed to be in ascending order.
    // This means:
    //     a) you never have to store the previous walker position in your child class and
    //     b) you can use the indices in changedIdx to provide efficient recalculation of your protovalues
    // Remember that in this method you should only update the protov[] elements that need to change due
    // to the nchanged input indices in changedIdx.
    // If full recalculation is more efficient in your case, you may also choose not to overwrite this method.
    virtual double updatedAcceptance(const WalkerState &wlk,
                                     const double protoold[], double protonew[] /* update this! */)
    {
        // default to "calculate all"
        this->protoFunction(wlk.xnew, protonew);
        return this->acceptanceFunction(protoold, protonew);
    }
};
}  // namespace mci

#endif
