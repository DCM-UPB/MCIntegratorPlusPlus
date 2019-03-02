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
class MCIFullAccumulator
{
protected:
    // we store these for fast access
    int _nstored;
    int _ndata;

    int _storeidx; // storage index offset for next write

    // --- storage method to be implemented
    void _allocate() override
    {
        _nstored = this->getNAccu();
        _ndata = this->getNData();
        _data = new double[_ndata]; // _nstored * _nobs layout
        std::fill(_data, _data+_ndata, 0.);
    }

    void _accumulate() override
    {
        for (int i=0; i<_nobs; ++i) {
            _data[_storeidx + i] += _obs->getObs(i);
        }
        _storeidx += _nobs;
    }

    void _reset() override
    {
        _storeidx = 0;
        std::fill(_data, _data+_ndata, 0.);
    }

    void _deallocate() override
    {
        delete [] _data;
        _nstored = 0;
        _ndata = 0;
    }

public:
    MCIFullAccumulator(MCIObservableFunctionInterface * obs, int nskip):
        MCIAccumulatorInterface(obs, nskip), _nstored(0), _ndata(0), _storeidx(0)
    {}

    ~MCIFullAccumulator() override { this->_deallocate(); }


    int getNStored() override { return _nstored; }

    void finalize() override {} // nothing to do
};


#endif
