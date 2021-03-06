#include "mci/BlockAccumulator.hpp"

namespace mci
{

void BlockAccumulator::_allocate()
{
    if (this->getNAccu() < _blocksize) {
        throw std::invalid_argument("[BlockAccumulator::allocate] Requested number of accumulations is smaller than the requested block size.");
    }
    if (this->getNAccu()%_blocksize != 0) {
        throw std::invalid_argument("[BlockAccumulator::allocate] Requested number of accumulations is not a multiple of the requested block size.");
    }
    _nblocks = this->getNAccu()/_blocksize;
    _data = new double[this->getNData()]; // _nstore * _nobs layout
    std::fill(_data, _data + this->getNData(), 0.);
}


void BlockAccumulator::_accumulate()
{
    for (int i = 0; i < _nobs; ++i) {
        _data[_storeidx + i] += _obs_values[i];
    }

    if (++_bidx == _blocksize) {
        _bidx = 0;
        _storeidx += _nobs; // move to next block
    }
}


void BlockAccumulator::_finalize()
{
    const double normf = 1./_blocksize;
    for (int64_t i = 0; i < _nblocks; ++i) {
        for (int j = 0; j < _nobs; ++j) {
            _data[i*_nobs + j] *= normf;
        }
    }
}


void BlockAccumulator::_reset()
{
    _bidx = 0;
    _storeidx = 0;
    std::fill(_data, _data + this->getNData(), 0.);
}


void BlockAccumulator::_deallocate()
{
    delete[] _data;
    _data = nullptr;
    _nblocks = 0;
}
}  // namespace mci
