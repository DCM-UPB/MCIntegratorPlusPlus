#ifndef MCI_MCIFULLACCUMULATOR_HPP
#define MCI_MCIFULLACCUMULATOR_HPP

#include "mci/MCIAccumulatorInterface.hpp"

#include <algorithm>

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
    // we store these for fast access
    int _nstore;
    int _ndata;

    int _storeidx; // storage index offset for next write

    // --- storage method to be implemented
    void _allocate() override
    {
        _nstore = this->getNAccu();
        _data = new double[this->getNData()]; // _nstore * _nobs layout
        std::fill(_data, _data+this->getNData(), 0.);
    }

    void _accumulate() override
    {   // no checks here for performance
        for (int i=0; i<_nobs; ++i) {
            _data[_storeidx + i] += _obs->getObservable(i);
        }
        _storeidx += _nobs;
    }

    void _reset() override
    {
        _storeidx = 0;
        std::fill(_data, _data+this->getNData(), 0.);
    }

    void _deallocate() override
    {
        delete [] _data;
        _data = nullptr;
        _nstore = 0;
    }

public:
    MCIFullAccumulator(MCIObservableFunctionInterface * obs, int nskip):
        MCIAccumulatorInterface(obs, nskip), _nstore(0), _storeidx(0)
    {}

    ~MCIFullAccumulator() override { this->_deallocate(); }

    int getNStore() override { return _nstore; }

    void finalize() override {} // nothing to do
};


#endif
