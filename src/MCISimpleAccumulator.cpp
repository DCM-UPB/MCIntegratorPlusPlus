#include "mci/MCISimpleAccumulator.hpp"

#include <algorithm>


void MCISimpleAccumulator::_allocate()
{
    _data = new double[_nobs];
    std::fill(_data, _data+_nobs, 0.);
    _flag_alloc = true;
}


void MCISimpleAccumulator::_accumulate()
{
    for (int i=0; i<_nobs; ++i) {
        _data[i] += _obs->getObservable(i);
    }
}


void MCISimpleAccumulator::_finalize()
{   // do nothing on deallocated state
    if (_flag_alloc) {
        const double normf = 1./this->getNAccu();
        for (int i=0; i<_nobs; ++i) { _data[i] *= normf; }
    }
}


void MCISimpleAccumulator::_reset()
{   // reset must not fail on deallocated state
    if (_flag_alloc) { std::fill(_data, _data+_nobs, 0.); }
}


void MCISimpleAccumulator::_deallocate()
{
    delete [] _data;
    _data = nullptr;
    _flag_alloc = false;
}
