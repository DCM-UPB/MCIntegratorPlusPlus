#ifndef MCI_OBSERVABLEFUNCTIONINTERFACE_HPP
#define MCI_OBSERVABLEFUNCTIONINTERFACE_HPP

#include "mci/Clonable.hpp"

#include <algorithm>
#include <vector>
namespace mci
{
    // Base class for MC observables
    //
    // Derive from this and implement one or both observableFunction() methods (see below).
    // You also need to provide a protected _clone method returning a raw pointer of
    // type ObservableFunctionInterface, pointing to an object of your type MyObservable, e.g.:
    //
    // class MyObservable: public ObservableFunctionInterface {
    // protected:
    //     ObservableFunctionInterface * _clone() const override {
    //         return new MyObservable(...); // create a cloned version here
    //     }
    // public:
    //     void observableFunction(...) overwrite;
    //     ...
    // };
    //
    // Your class will have a public clone() method returning std::unique_ptr<ObservableFunctionInterface> .
    // If you want/need it, also create a non-overriding clone() method returning a pointer of type MyObservable.
    class ObservableFunctionInterface: public Clonable<ObservableFunctionInterface>
    {
    protected:
        const int _ndim;  //dimension of the input array (walker position)
        const int _nobs;  //number of values provided by the observable

    public:
        ObservableFunctionInterface(int ndim, int nobs): _ndim(ndim), _nobs(nobs) {}

        // getters
        int getNObs() const { return _nobs; }
        int getNDim() const { return _ndim; }

        // --- METHOD THAT MUST BE IMPLEMENTED
        // Compute all observable elements and store them in out.
        virtual void observableFunction(const double in[], double out[]) = 0;
        //                               ^input = walker positions  ^resulting observables

        // --- OPTIONALLY ALSO OVERWRITE THIS
        // Compute the observable, given ndim flags indicating which inputs have changed since last observable calculation.
        // This means you don't need to store the previous walker position to decide how to efficiently calculate
        // the new observable, unless you really need the old positions in your calculation. The last observable values are
        // passed via the out array, so again you can use these values without storing internal "old" values, unless you explictly
        // need to do so.
        // You may use the nchanged argument to decide whether a full recalculation or flag-based recalculation is more efficient.
        // If full recalculation is almost always more efficient in your case, simply don't overwrite this method.
        virtual void observableFunction(const double in[], const int /*nchanged*/, const bool[]/*flags[ndim]*/, double out[]){ observableFunction(in, out); }
        //                               ^input = walker positions  ^how many inputs changed  ^which indices are new   ^resulting observables (passed containing old obs, so you can make use of those)
    };
}  // namespace mci


#endif
