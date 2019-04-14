#ifndef MCI_OBSERVABLEFUNCTIONINTERFACE_HPP
#define MCI_OBSERVABLEFUNCTIONINTERFACE_HPP

#include "mci/Clonable.hpp"
#include "mci/SamplingFunctionContainer.hpp"

#include <algorithm>

namespace mci
{
// Base class for MC observables
//
// Derive from this and implement the observableFunction() methods (see below).
// You also need to provide a protected _clone method returning a raw pointer of
// type ObservableFunctionInterface, pointing to an object of your type MyObservable, e.g.:
//
// class MyObservable: public ObservableFunctionInterface {
// protected:
//     ObservableFunctionInterface * _clone() const final {
//         return new MyObservable(...); // create a cloned version here
//     }
// public:
//     void observableFunction(...) final;
//     ...
// };
//
// Your class will have a public clone() method returning std::unique_ptr<ObservableFunctionInterface> .
//
// NOTE: If you are using single-particle moves and your observable allows efficient recomputation, given
// knowledge about changed input indices since last computation, you may pass isUpdateable=true to the
// constructor and override updatedObservable(..).
//
class ObservableFunctionInterface: public Clonable<ObservableFunctionInterface>
{
protected:
    const int _ndim;  //dimension of the input array (walker position)
    const int _nobs;  //number of values provided by the observable
    const bool _flag_updateable; // does the obs want to enable use of the updatedObservable() method

    ObservableFunctionInterface(int ndim, int nobs, bool isUpdateable): _ndim(ndim), _nobs(nobs), _flag_updateable(isUpdateable) {}

public:
    // getters
    int getNObs() const { return _nobs; }
    int getNDim() const { return _ndim; }
    bool isUpdateable() const { return _flag_updateable; }

    // --- METHOD THAT MUST BE IMPLEMENTED
    // Compute all observable elements and store them in out.
    virtual void observableFunction(const double in[], const SamplingFunctionContainer &pdfcont, double out[]) = 0;
    //                               ^input = walker positions            ^you may read the pdfs        ^resulting observables


    // --- YOU MAY ALSO OVERRIDE THIS (only used if flag_updateable set to true)
    // Compute the observable, given ndim flags indicating which inputs have changed since last observable calculation.
    // This means you don't need to store the previous walker position to decide how to efficiently calculate
    // the new observable, unless you really need the old positions in your calculation. The last observable values are
    // passed via the out array, so again you can use these values without storing internal "old" values, unless you explictly
    // need to do so.
    // You may use the nchanged argument to decide whether a full recalculation or flag-based recalculation is more efficient.
    // If full recalculation is almost always more efficient in your case, you may also choose not to override this method and
    // pass isUpdateable=false to the constructor.
    virtual void updatedObservable(const double in[], int/*nchanged*/, const bool/*flags_xchanged[ndim]*/[], const SamplingFunctionContainer &pdfcont, double out[]) { this->observableFunction(in, pdfcont, out); }
    //                             ^input = walker positions  ^how many inputs changed  ^which indices are new           ^allows to read from pdfs      ^resulting observables (passed containing old obs, so you may make use of those)
};
}  // namespace mci


#endif
