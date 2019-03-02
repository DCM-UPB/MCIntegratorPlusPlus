#ifndef MCI_MCIBLOCKACCUMULATOR_HPP
#define MCI_MCIBLOCKACCUMULATOR_HPP

#include "mci/MCIAccumulatorInterface.hpp"

#include <algorithm>
#include <stdexcept>

// Class to handle accumulation of observables, when averaging samples in blocks of fixed size
// is desired. Typical use case is if you know how large the blocks have to be for uncorrelated samples
// and want to avoid the memory&CPU overhead of using automatic blocking.
//
class MCIBlockAccumulator: public MCIAccumulatorInterface
{
protected:
    const int _blocksize;

    // we store these for fast access
    int _nblocks;
    int _ndata;

    int _bidx; // counter to determine when block is finished
    int _storeidx; // storage index offset for next write

    // --- storage method to be implemented
    void _allocate() override
    {
        if (this->getNAccu() % _blocksize != 0) {
            throw std::invalid_argument("[MCIBlockAccumulator::allocate] Requested number of accumulations is not a multiple of the requested block size.");
        }
        _nblocks = this->getNAccu() / _blocksize;
        _ndata = this->getNData();
        _data = new double[_ndata]; // _nstored * _nobs layout
        std::fill(_data, _data+_ndata, 0.);
    }

    void _accumulate() override
    {
        for (int i=0; i<_nobs; ++i) {
            _data[_storeidx + i] += _obs->getObservable(i);
        }

        if (++_bidx == _blocksize) {
            _bidx = 0;
            _storeidx += _nobs; // move to next block
        }
    }

    void _reset() override
    {
        _bidx = 0;
        _storeidx = 0;
        std::fill(_data, _data+_ndata, 0.);
    }

    void _deallocate() override
    {
        delete [] _data;
        _nblocks = 0;
        _ndata = 0;
    }

public:
    MCIBlockAccumulator(MCIObservableFunctionInterface * obs, int nskip, int blocksize):
        MCIAccumulatorInterface(obs, nskip), _blocksize(blocksize), _nblocks(0), _ndata(0), _bidx(0), _storeidx(0)
    {
        if (_blocksize < 1) { throw std::invalid_argument("[MCIBlockAccumulator] Requested blocksize was < 1 ."); }
    }

    ~MCIBlockAccumulator() override { this->_deallocate(); }

    int getBlockSize() { return _blocksize; }
    int getNStored() override { return _nblocks; }

    void finalize() override
    {
        const double normf = 1./_blocksize;
        for (int i=0; i<_nblocks; ++i) {
            for (int j=0; j<_nobs; ++j) {
                _data[i*_nobs + j] *= normf;
            }
        }
    }
};


#endif