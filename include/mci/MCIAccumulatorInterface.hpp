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

    double * _data; // childs use this to store data

    // --- storage method to be implemented
    virtual void _allocate() = 0; // allocate memory for a MC run of nsteps length
    virtual void _accumulate() = 0; // store new observable data
    virtual void _reset() = 0; // reset data/counters for new accumulation (but keep allocation)
    virtual void _deallocate() = 0; // free data memory

public:
    MCIAccumulatorInterface(MCIObservableFunctionInterface * obs, int nskip);

    virtual ~MCIAccumulatorInterface() = default;

    MCIObservableFunctionInterface * getObservableFunction(){ return _obs; }

    int getNObs(){ return _nobs; }
    int getNSkip(){ return _nskip; }

    int getNSteps(){ return _nsteps; }
    int getNAccu(){ return (_nsteps>0) ? 1 + (_nsteps-1)/_nskip : 0; } // actual number of steps to accumulate

    int getNData(){ return this->getNStore()*_nobs; } // total length of allocated data
    const double * getData(){ return _data; } // direct read-only access to data

    // externally call this before a MC run of nsteps length
    void allocate(int nsteps);

    // externally call this on every MC step
    void accumulate(const double * in, bool flagacc);

    // reset for new accumulation
    void reset();

    // deallocate memory
    void deallocate();

    // methods to be implemented
    virtual int getNStore() = 0; // get number of allocated data elements with _nobs length each
    virtual void finalize() = 0; // apply post-processing (usually normalization) to stored data
};


#endif
