#include "mci/MCIFullAccumulator.hpp"

#include <algorithm>

void MCIFullAccumulator::_allocate()
{
    _nstore = this->getNAccu();
    _data = new double[this->getNData()]; // _nstore * _nobs layout
    std::fill(_data, _data+this->getNData(), 0.); // not strictly necessary
}


void MCIFullAccumulator::_accumulate()
{
    std::copy(_obs_values, _obs_values+_nobs, _data+_storeidx);
    _storeidx += _nobs;
}


void MCIFullAccumulator::_reset()
{
    _storeidx = 0;
    std::fill(_data, _data+this->getNData(), 0.);
}


void MCIFullAccumulator::_deallocate()
{
    delete [] _data;
    _data = nullptr;
    _nstore = 0;
}
