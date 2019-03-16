#ifndef MCI_UPDATEABLEOBSERVABLEINTERFACE_HPP
#define MCI_UPDATEABLEOBSERVABLEINTERFACE_HPP

#include "mci/ObservableFunctionInterface.hpp"

namespace mci
{
    // Base class for MC observables with partial update method
    //
    // Derive from this if your observable allows efficient recomputation,
    // given knowledge about changed input indices since last computation,
    // in the form of a boolean array with one bool per input index (true if changed).
    // Besides the methods from ObservableFunctionInterface, you also need to
    // implement the mentioned partial update method.
    // So typically your class would look like the following:
    //
    // class MyObservable: public UpdateableObservableInterface {
    // protected:
    //     ObservableFunctionInterface * _clone() const final {
    //         // returns ptr of type ObservableFunctionInterface
    //         return new MyObservable(...); // create a cloned version here
    //     }
    // public:
    //     void observableFunction(in, out) final; // full update
    //     void observableUpdate(in, nchanged, changeFlags, out) final; // partial update, see below
    //     ...
    // };
    class UpdateableObservableInterface: public ObservableFunctionInterface
    {
    public:
        UpdateableObservableInterface(int ndim, int nobs): ObservableFunctionInterface(ndim ,nobs) {}

        // --- YOU MUST IMPLEMENT THIS
        // Compute the observable, given ndim flags indicating which inputs have changed since last observable calculation.
        // This means you don't need to store the previous walker position to decide how to efficiently calculate
        // the new observable, unless you really need the old positions in your calculation. The last observable values are
        // passed via the out array, so again you can use these values without storing internal "old" values, unless you explictly
        // need to do so.
        // You may use the nchanged argument to decide whether a full recalculation or flag-based recalculation is more efficient.
        // If full recalculation is almost always more efficient in your case, you may also choose not to overwrite this method.
        virtual void updatedObservable(const double in[], int nchanged, const bool flags_xchanged[]/*[ndim]*/, double out[]) = 0;
        //                             ^input = walker positions  ^how many inputs changed  ^which indices are new  ^resulting observables (passed containing old obs, so you may make use of those)
    };

}  // namespace mci


#endif
