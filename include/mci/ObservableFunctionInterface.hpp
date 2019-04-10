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
    virtual void observableFunction(const double in[], const SamplingFunctionContainer &pdfcont, double out[]) = 0;
    //                               ^input = walker positions            ^you may read the pdfs        ^resulting observables
};
}  // namespace mci


#endif
