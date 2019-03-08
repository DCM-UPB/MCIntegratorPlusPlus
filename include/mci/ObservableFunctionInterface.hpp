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
    // type ObservableFunctionInterface, for example:
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
        const int _ndim;  //dimension of the input array (walker poistion)
        const int _nobs;  //number of values provided by the observable
        double * const _obs; //array that stores the last observables computed

    public:
        ObservableFunctionInterface(int ndim, int nobs): _ndim(ndim), _nobs(nobs), _obs(new double[nobs])
        {
            std::fill(_obs, _obs+nobs, 0.);
        }

        ~ObservableFunctionInterface() override{ delete [] _obs; }

        int getNObs() const { return _nobs; }
        int getNDim() const { return _ndim; }
        double getValue(int i) const { return _obs[i]; }
        const double * getValues() const { return _obs; }

        void computeValues(const double in[])
        {
            observableFunction(in, _obs);
        }

        void computeValues(const double in[], const bool flags_xchanged[] /*flag_i true if x_i changed*/)
        {
            observableFunction(in, flags_xchanged, _obs);
        }


        // --- METHOD THAT MUST BE IMPLEMENTED
        // Compute the observable and store it in out
        virtual void observableFunction(const double in[], double out[]) = 0;
        //                               ^input = walker positions  ^resulting observables

        // --- OPTIONALLY ALSO OVERWRITE THIS
        // Compute the observable, given knowledge of what changed
        virtual void observableFunction(const double in[], const bool /*flags*/[], double out[]){ observableFunction(in, out); }
        //                               ^input = walker positions   ^which input is new   ^resulting observables
    };
}  // namespace mci


#endif
