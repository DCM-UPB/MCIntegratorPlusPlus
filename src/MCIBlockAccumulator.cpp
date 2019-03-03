#include "mci/MCIBlockAccumulator.hpp"

#include <algorithm>
#include <stdexcept>

void MCIBlockAccumulator::_allocate()
{
    if (this->getNAccu() % _blocksize != 0) {
        throw std::invalid_argument("[MCIBlockAccumulator::allocate] Requested number of accumulations is not a multiple of the requested block size.");
    }
    _nblocks = this->getNAccu() / _blocksize;
    _data = new double[this->getNData()]; // _nstore * _nobs layout
    std::fill(_data, _data+this->getNData(), 0.);
}


void MCIBlockAccumulator::_accumulate()
{
    for (int i=0; i<_nobs; ++i) {
        _data[_storeidx + i] += _obs->getObservable(i);
    }

    if (++_bidx == _blocksize) {
        _bidx = 0;
        _storeidx += _nobs; // move to next block
    }
}


void MCIBlockAccumulator::_finalize()
{
    const double normf = 1./_blocksize;
    for (int i=0; i<_nblocks; ++i) {
        for (int j=0; j<_nobs; ++j) {
            _data[i*_nobs + j] *= normf;
        }
    }
}


void MCIBlockAccumulator::_reset()
{
    _bidx = 0;
    _storeidx = 0;
    std::fill(_data, _data+this->getNData(), 0.);
}


void MCIBlockAccumulator::_deallocate()
{
    delete [] _data;
    _data = nullptr;
    _nblocks = 0;
}
