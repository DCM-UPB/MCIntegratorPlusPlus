#ifndef MCI_MCIFULLACCUMULATOR_HPP
#define MCI_MCIFULLACCUMULATOR_HPP

#include "mci/MCIAccumulatorInterface.hpp"

// Class to handle accumulation of observables, when storing every single sample is desired
// Typically you want this for automatic blocking techniques after the sampling run
//
// NOTE: If the dimension of the observable is very large, consider using MCIBlockAccumulator
// or MCISimpleAccumulator (if no error is required) instead, because the memory requirements
// of the MCIFullAccumulator may become very large with a large number of MC steps.
//
class MCIFullAccumulator: public MCIAccumulatorInterface
{
protected:
    int _nstore; // number of allocated storage elements with _nobs length each
    int _storeidx; // storage index offset for next write

    // --- storage method to be implemented
    void _allocate() override;
    void _accumulate() override;
    void _finalize() override {} // nothing to do
    void _reset() override;
    void _deallocate() override;

public:
    MCIFullAccumulator(MCIObservableFunctionInterface * obs, int nskip):
        MCIAccumulatorInterface(obs, nskip), _nstore(0), _storeidx(0)
    {}

    ~MCIFullAccumulator() override { this->_deallocate(); }

    int getNStore() override { return _nstore; }
};


#endif
