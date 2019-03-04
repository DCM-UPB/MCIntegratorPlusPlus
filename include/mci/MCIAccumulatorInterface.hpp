#ifndef MCI_MCIACCUMULATORINTERFACE_HPP
#define MCI_MCIACCUMULATORINTERFACE_HPP

#include "mci/MCIObservableFunctionInterface.hpp"

// Class to handle accumulation of data for a given observable function
class MCIAccumulatorInterface
{
protected:
    MCIObservableFunctionInterface * const _obs; // pointer to the related observable function
    const int _nobs; // number of values returned by the observable function

    const int _nskip; // evaluate observable only on every nskip-th step
    int _nsteps; // total number of sampling steps (expected number of calls to accumulateObservables)

    int _stepidx; // running step index
    int _accuidx; // counting number of non-skipped accumulations
    int _skipidx; // to determine when to skip accumulation
    bool _flag_eval; // observable evaluation is needed before next accumulation
    bool _flag_final; // was finalized called (without throwing error) ?

    double * _data; // childs use this to store data

    // TO BE IMPLEMENTED BY CHILD
    virtual void _allocate() = 0; // allocate _data for a MC run of nsteps length ( expect deallocated state )
    virtual void _accumulate() = 0; // store new observable data ( sanity checks are done already )
    virtual void _finalize() = 0; // if necessary, apply normalization ( do nothing when deallocated )
    virtual void _reset() = 0; // reset data / child's members ( must work in deallocated state )
    virtual void _deallocate() = 0; // delete _data allocation ( reset will be called already )

public:
    MCIAccumulatorInterface(MCIObservableFunctionInterface * obs, int nskip);

    virtual ~MCIAccumulatorInterface() = default;

    MCIObservableFunctionInterface * getObservableFunction(){ return _obs; }

    // Getters
    int getNObs(){ return _nobs; }
    int getNSkip(){ return _nskip; }
    int getNSteps(){ return _nsteps; }
    int getNAccu(){ return (_nsteps>0) ? 1 + (_nsteps-1)/_nskip : 0; } // actual number of steps to accumulate
    int getNData(){ return this->getNStore()*_nobs; } // total length of allocated data

    int getStepIndex(){ return _stepidx; }
    bool isAllocated(){ return (_nsteps>0); }
    bool isClean(){ return (_stepidx == 0); }
    bool isFinalized(){ return _flag_final; }

    // TO BE IMPLEMENTED BY CHILD
    virtual int getNStore() = 0; // get number of allocated data elements with _nobs length each


    // methods to call externally, in the following pattern:
    // allocate -> nsteps * accumulate -> finalize -> getData ( -> reset -> accumulate ...) -> delete/deallocate

    // call this before a MC run of nsteps length
    void allocate(int nsteps); // will deallocate an existing allocation

    // externally call this on every MC step
    void accumulate(const double * in, bool flagacc); // will throw if not allocated (enough)

    // finalize (e.g. normalize) stored data
    void finalize(); // will throw if called prematurely, but does nothing if deallocated or used repeatedly

    // get data
    const double * getData(){ return _data; } // direct read-only access to internal data pointer

    // reset for new accumulation
    void reset();

    // deallocate memory
    void deallocate();
};


#endif
