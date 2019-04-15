#ifndef MCI_OBSERVABLECONTAINER_HPP
#define MCI_OBSERVABLECONTAINER_HPP

#include "mci/ObservableFunctionInterface.hpp"
#include "mci/DependentObservableInterface.hpp"
#include "mci/Factories.hpp"
#include "mci/WalkerState.hpp"
#include "mci/SamplingFunctionContainer.hpp"

#include <cstdint>
#include <fstream>
#include <functional>
#include <memory>
#include <vector>

namespace mci
{

// Internally used container for objects related to observable accumulation
class ObservableContainer
{
    struct ObservableContainerElement
    {
        // Observable
        std::unique_ptr<ObservableFunctionInterface> obs;

        // obs ptr dynamic_casted to DependentObservableInterface (if possible, else nullptr)
        DependentObservableInterface * depobs;

        // Accumulator
        std::unique_ptr<AccumulatorInterface> accu;

        // Estimator function used to obtain result of MC integration
        std::function<void(double [] /*avg*/, double [] /*error*/)> estim; // corresponding accumulator is already bound

        // flags
        bool flag_equil{}; // equilibrate this observable when using automatic decorrelation?
    };

private:
    // vector with container elements
    std::vector<ObservableContainerElement> _cont;
    int _nobsdim{0}; // stores total dimension of contained observables

public:
    explicit ObservableContainer() = default;
    ~ObservableContainer() = default;

    // simple getters
    int size() const { return static_cast<int>(_cont.size()); }
    int getNObs() const { return this->size(); }
    int getNObsDim() const { return _nobsdim; }

    bool empty() const { return _cont.empty(); }
    bool hasObs() const { return !this->empty(); }

    ObservableFunctionInterface &getObservableFunction(int i) const { return *(_cont[i].obs); }
    const AccumulatorInterface &getAccumulator(int i) const { return *(_cont[i].accu); }
    bool getFlagEquil(int i) const { return _cont[i].flag_equil; }

    // operational methods
    // add observable (+internally accumulator&estimator)
    void addObservable(std::unique_ptr<ObservableFunctionInterface> obs /*we acquire ownership*/,
                       int blocksize, int nskip, bool needsEquil, EstimatorType estimType);

    void allocate(int64_t Nmc, const SamplingFunctionContainer &pdfcont); // allocate data memory and register dependencies
    void accumulate(const WalkerState &wlk); // process accumulation for new step, described by WalkerState
    void printObsValues(std::ofstream &file) const; // write last observables values to filestream
    void finalize(); // used after sampling to apply all necessary data normalization
    void estimate(double average[], double error[]) const; // eval estimators on finalized data and return average/error
    void reset(); // obtain clean state, but keep allocation
    void deallocate(); // free data memory
    std::unique_ptr<ObservableFunctionInterface> pop_back(); // remove and return last obs
    void clear(); // clear everything
};
}  // namespace mci

#endif
