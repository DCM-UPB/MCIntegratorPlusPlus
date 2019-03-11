#ifndef MCI_OBSERVABLECONTAINER_HPP
#define MCI_OBSERVABLECONTAINER_HPP

#include "mci/AccumulatorInterface.hpp"
#include "mci/ObservableFunctionInterface.hpp"

#include <fstream>
#include <functional>
#include <memory>
#include <vector>

namespace mci
{
    class ObservableContainer
    {
    private:
        int _nobsdim {0}; // stores total dimension of contained observables

        // Accumulators
        std::vector< std::unique_ptr<AccumulatorInterface> > _accus;

        // Estimator functions used to obtain result of MC integration
        std::vector< std::function< void (double [] /*avg*/, double [] /*error*/) > > _estims; // corresponding accumulators are already bound

        // flags
        std::vector< unsigned char /*avoid special vector<bool>..*/ > _flags_equil; // equilibrate this observable when using automatic decorrelation?

    public:
        explicit ObservableContainer() = default;
        ~ObservableContainer() = default;

        // simple getters
        int getNObs() const { return _accus.size(); }
        int getNObsDim() const { return _nobsdim; }

        const ObservableFunctionInterface & getObservableFunction(int i) const { return _accus[i]->getObservableFunction(); }
        bool getFlagEquil(int i) const {return (_flags_equil[i]>0);}

        // operational methods
        // add accumulator&estimator for an observable
        void addObservable(std::unique_ptr<AccumulatorInterface> accumulator, // we acquire ownerhsip
                           const std::function< void (int /*nstored*/, int /*nobs*/, const double [] /*data*/, double [] /*avg*/, double [] /*error*/) > &estimator,
                           bool needsEquil);

        void allocate(int Nmc); // allocate data memory
        void accumulate(const double x[], int nchanged, const int changedIdx[]); // process accumulation for position x, with nchanged indices changedIdx
        void printObsValues(std::ofstream &file) const; // write last observables values to filestream
        void finalize(); // used after sampling to apply all necessary data normalization
        void estimate(double average[], double error[]) const; // eval estimators on finalized data and return average/error
        void reset(); // obtain clean state, but keep allocation
        void deallocate(); // free data memory
        void clear(); // clear everything
    };
}  // namespace mci

#endif
