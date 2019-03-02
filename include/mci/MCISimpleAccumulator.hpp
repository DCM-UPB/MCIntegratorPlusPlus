#ifndef MCI_MCISIMPLEACCUMULATOR_HPP
#define MCI_MCISIMPLEACCUMULATOR_HPP

#include "mci/MCIAccumulatorInterface.hpp"

#include <algorithm>

// Class to handle accumulation of observables that don't require an error (typically gradients)
// This mean that all we have to do is to sum up the individual samples
class MCISimpleAccumulator
{
protected:
    // --- storage method to be implemented
    void _allocate() override
    {
        _data = new double[_nobs];
        std::fill(_data, _data+_nobs, 0.);
    }

    void _accumulate() override
    {
        for (int i=0; i<_nobs; ++i) {
            _data[i] += _obs->getObs(i);
        }
    }

    void _reset() override
    {
        std::fill(_data, _data+_nobs, 0.);
    }

    void _deallocate() override
    {
        delete [] _data;
    }

public:
    MCISimpleAccumulator(MCIObservableFunctionInterface * obs, int nskip):
        MCIAccumulatorInterface(obs, nskip) {}

    ~MCISimpleAccumulator() override
    {
        this->_deallocate();
    }


    int getNStored() override {
        return 1; // we don't store old values
    }

    void finalize() override
    {
        const double normf = 1./this->getNAccu();
        for (int i=0; i<_nobs; ++i) { _data[i] *= normf; }
    }
};


#endif
