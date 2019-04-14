#ifndef MCI_ACCUMULATORINTERFACE_HPP
#define MCI_ACCUMULATORINTERFACE_HPP

#include "mci/ObservableFunctionInterface.hpp"
#include "mci/WalkerState.hpp"
#include "mci/SamplingFunctionContainer.hpp"

#include <cstdint>
#include <memory>

namespace mci
{
// Interface class to handle accumulation of observable data
//
// Accumulators completely wrap around an exclusively owned Observable-
// FunctionInterface and are the "communication partner" for MCI during
// sampling. The derived classes of this interface implement different
// storage/accumulation techniques. Accumulators are typically contained
// within an Observable-Container (see ObservableContainer.hpp), where they
// are strictly paired with corresponding average/error estimation functions.
class AccumulatorInterface
{
protected:
    std::unique_ptr<ObservableFunctionInterface> _obs; // "unique" pointer to the passed observable function (we own it)
    const bool _flag_updobs; // is the passed observable supporting selective updating?

    const int _nobs; // number of values returned by the observable function
    const int _xndim; // dimension of walker positions/flags that get passed on accumulate
    const int _nskip; // evaluate observable only on every nskip-th step

    // fixed-size allocations
    double * const _obs_values; // observable's last values (length _nobs)
    bool * const _flags_xchanged; // remembers which x have changed since last obs evaluation (length _xndim)

    // variables
    int64_t _nsteps; // total number of sampling steps (set on allocate() to planned number of calls to accumulateObservables)
    double * _data; // childs use this to store data

    int _nchanged{}; // counter of how many x have changed since last obs evaluation
    int64_t _stepidx{}; // running step index
    int _skipidx{}; // to determine when to skip accumulation
    bool _flag_final{}; // was finalized called (without throwing error) ?

    // base methods
    void _init(); // used in construct/reset
    void _processOld(const WalkerState &wlk, const SamplingFunctionContainer &pdfcont); // used in accumulate() when observables need no computation
    void _processFull(const WalkerState &wlk, const SamplingFunctionContainer &pdfcont); // used else when obs not updateable
    void _processSelective(const WalkerState &wlk, const SamplingFunctionContainer &pdfcont); // and this is used otherwise

    // TO BE IMPLEMENTED BY CHILD
    virtual void _allocate() = 0; // allocate _data for a MC run of nsteps length ( expect deallocated state )
    virtual void _accumulate() = 0; // store new observable data ( sanity checks are done already )
    virtual void _finalize() = 0; // if necessary, apply normalization ( do nothing when deallocated )
    virtual void _reset() = 0; // reset data / child's members ( must work in deallocated state )
    virtual void _deallocate() = 0; // delete _data allocation ( reset will be called already )

    // Constructor
    AccumulatorInterface(std::unique_ptr<ObservableFunctionInterface> obs, int nskip);

public:
    virtual ~AccumulatorInterface();
    std::unique_ptr<ObservableFunctionInterface> removeObs(); // remove and return the contained obs. NOTE: After calling this, the accumulator should be immediately deleted.

    ObservableFunctionInterface &getObservableFunction() const { return *_obs; } // acquire raw read-only ref

    // Getters
    int getNObs() const { return _nobs; } // dimension of observable
    int getNDim() const { return _xndim; } // dimension of walkers

    int getNSkip() const { return _nskip; }
    int64_t getNSteps() const { return _nsteps; }
    int64_t getNAccu() const
    {
        return (_nsteps > 0) ? 1 + (_nsteps - 1)/_nskip : 0; // actual number of steps to accumulate
    }
    int64_t getNData() const { return this->getNStore()*_nobs; } // total length of allocated data

    int64_t getStepIndex() const { return _stepidx; }
    bool isAllocated() const { return (_nsteps > 0); }
    bool isClean() const { return (_stepidx == 0); }
    bool isFinalized() const { return _flag_final; }
    bool isUpdateable() const { return _flag_updobs; }

    // get data
    const double * getData() const { return _data; } // direct read-only access to internal data pointer
    const double * getObsValues() const { return _obs_values; } // read-only pointer to last calculated observable data
    double getObsValue(int i) const { return _obs_values[i]; } // element-wise access to last values

    // TO BE IMPLEMENTED BY CHILD
    virtual int64_t getNStore() const = 0; // get number of allocated data elements with _nobs length each


    // methods to call externally, in the following pattern:
    // allocate -> nsteps * accumulate -> finalize -> getData ( -> reset -> accumulate ...) -> delete/deallocate

    // call this before a MC run of nsteps length
    void allocate(int64_t nsteps); // will deallocate any existing allocation

    // externally call this on every MC step
    void accumulate(const WalkerState &wlk /*step info*/, const SamplingFunctionContainer &pdfcont /*pdf info*/); // process step described by WalkerState

    // finalize (e.g. normalize) stored data
    void finalize(); // will throw if called prematurely, but does nothing if deallocated or used repeatedly

    // reset for new accumulation
    void reset();

    // deallocate memory
    void deallocate();
};
}  // namespace mci

#endif
