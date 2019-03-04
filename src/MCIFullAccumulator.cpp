#include "mci/MCIFullAccumulator.hpp"

#include <algorithm>


void MCIFullAccumulator::_allocate()
{
    _nstore = this->getNAccu();
    _data = new double[this->getNData()]; // _nstore * _nobs layout
    std::fill(_data, _data+this->getNData(), 0.);
}


void MCIFullAccumulator::_accumulate()
{
    for (int i=0; i<_nobs; ++i) {
        _data[_storeidx + i] += _obs->getValue(i);
    }
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
