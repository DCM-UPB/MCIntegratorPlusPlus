#include "mci/FullAccumulator.hpp"

namespace mci
{

void FullAccumulator::_allocate()
{
    _nstore = this->getNAccu();
    _data = new double[this->getNData()]; // _nstore * _nobs layout
    std::fill(_data, _data + this->getNData(), 0.); // not strictly necessary
}


void FullAccumulator::_accumulate()
{
    std::copy(_obs_values, _obs_values + _nobs, _data + _storeidx);
    _storeidx += _nobs;
}


void FullAccumulator::_reset()
{
    _storeidx = 0;
    std::fill(_data, _data + this->getNData(), 0.);
}


void FullAccumulator::_deallocate()
{
    delete[] _data;
    _data = nullptr;
    _nstore = 0;
}
}  // namespace mci
