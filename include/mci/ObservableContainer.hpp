#ifndef MCI_OBSERVABLECONTAINER_HPP
#define MCI_OBSERVABLECONTAINER_HPP

#include "mci/AccumulatorInterface.hpp"
#include "mci/ObservableFunctionInterface.hpp"
#include "mci/WalkerState.hpp"

#include <cstdint>
#include <fstream>
#include <functional>
#include <memory>
#include <vector>

namespace mci
{
    class ObservableContainer
    {
        struct ObservableContainerElement
        {
            // Estimator function used to obtain result of MC integration
            std::function< void (double [] /*avg*/, double [] /*error*/) > estim; // corresponding accumulator is already bound

            // Accumulator
            std::unique_ptr<AccumulatorInterface> accu;

            // flags
            bool flag_equil{}; // equilibrate this observable when using automatic decorrelation?
        };

    private:
        // vector with container elements
        std::vector< ObservableContainerElement > _cont;
        int _nobsdim {0}; // stores total dimension of contained observables

    public:
        explicit ObservableContainer() = default;
        ~ObservableContainer() = default;

        // simple getters
        int size() const { return _cont.size(); }
        int getNObs() const { return this->size(); }
        int getNObsDim() const { return _nobsdim; }

        bool empty() const { return _cont.empty(); }
        bool hasObs() const { return !this->empty(); }

        ObservableFunctionInterface & getObservableFunction(int i) const { return _cont[i].accu->getObservableFunction(); }
        bool getFlagEquil(int i) const { return _cont[i].flag_equil; }

        // operational methods
        // add accumulator&estimator for an observable
        void addObservable(std::unique_ptr<AccumulatorInterface> accumulator, // we acquire ownerhsip
                           const std::function< void (int64_t /*nstored*/, int /*nobs*/, const double [] /*data*/, double [] /*avg*/, double [] /*error*/) > &estimator,
                           bool needsEquil);

        void allocate(int64_t Nmc); // allocate data memory
        void accumulate(const WalkerState &wlk); // process accumulation for new step, described by WalkerState
        void printObsValues(std::ofstream &file) const; // write last observables values to filestream
        void finalize(); // used after sampling to apply all necessary data normalization
        void estimate(double average[], double error[]) const; // eval estimators on finalized data and return average/error
        void reset(); // obtain clean state, but keep allocation
        void deallocate(); // free data memory
        std::unique_ptr<AccumulatorInterface> pop_back(); // remove and return last accu
        void clear(); // clear everything
    };
}  // namespace mci

#endif
