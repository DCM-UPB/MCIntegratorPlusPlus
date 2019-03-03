#ifndef MCI_MCISIMPLEACCUMULATOR_HPP
#define MCI_MCISIMPLEACCUMULATOR_HPP

#include "mci/MCIAccumulatorInterface.hpp"

// Class to handle accumulation of observables that don't require an error (typically gradients)
// This mean that all we have to do is to sum up the individual samples
class MCISimpleAccumulator: public MCIAccumulatorInterface
{
protected:
    bool _flag_alloc; // to determine proper nstored/ndata

    // --- storage method to be implemented
    void _allocate() override;
    void _accumulate() override;
    void _finalize() override;
    void _reset() override;
    void _deallocate() override;

public:
    MCISimpleAccumulator(MCIObservableFunctionInterface * obs, int nskip):
        MCIAccumulatorInterface(obs, nskip), _flag_alloc(false)
    {}

    ~MCISimpleAccumulator() override { this->_deallocate(); }

    int getNStore() override {
        return _flag_alloc ? 1 : 0; // we don't store old values
    }
};


#endif
