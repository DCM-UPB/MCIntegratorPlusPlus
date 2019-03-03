#ifndef MCI_MCISIMPLEACCUMULATOR_HPP
#define MCI_MCISIMPLEACCUMULATOR_HPP

#include "mci/MCIAccumulatorInterface.hpp"

#include <algorithm>

// Class to handle accumulation of observables that don't require an error (typically gradients)
// This mean that all we have to do is to sum up the individual samples
class MCISimpleAccumulator: public MCIAccumulatorInterface
{
protected:
    bool _flag_alloc; // to determine proper nstored/ndata

    // --- storage method to be implemented
    void _allocate() override
    {
        _data = new double[_nobs];
        std::fill(_data, _data+_nobs, 0.);
        _flag_alloc = true;
    }

    void _accumulate() override
    {   // no checks here for performance
        for (int i=0; i<_nobs; ++i) {
            _data[i] += _obs->getObservable(i);
        }
    }

    void _reset() override
    {   // reset must not fail on deallocated state
        if (_flag_alloc) { std::fill(_data, _data+_nobs, 0.); }
    }

    void _deallocate() override
    {
        delete [] _data;
        _data = nullptr;
        _flag_alloc = false;
    }

public:
    MCISimpleAccumulator(MCIObservableFunctionInterface * obs, int nskip):
        MCIAccumulatorInterface(obs, nskip), _flag_alloc(false) {}

    ~MCISimpleAccumulator() override
    {
        this->_deallocate();
    }


    int getNStore() override {
        return _flag_alloc ? 1 : 0; // we don't store old values
    }

    void finalize() override
    {
        if (_flag_alloc) {
            const double normf = 1./this->getNAccu();
            for (int i=0; i<_nobs; ++i) { _data[i] *= normf; }
        }
    }
};


#endif
