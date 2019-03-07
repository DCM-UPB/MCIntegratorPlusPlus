#include "mci/SimpleAccumulator.hpp"

#include <algorithm>

namespace mci
{

    void SimpleAccumulator::_allocate()
    {
        _data = new double[_nobs];
        std::fill(_data, _data+_nobs, 0.);
        _flag_alloc = true;
    }


    void SimpleAccumulator::_accumulate()
    {
        for (int i=0; i<_nobs; ++i) {
            _data[i] += _obs->getValue(i);
        }
    }


    void SimpleAccumulator::_finalize()
    {   // do nothing on deallocated state
        if (_flag_alloc) {
            const double normf = 1./this->getNAccu();
            for (int i=0; i<_nobs; ++i) { _data[i] *= normf; }
        }
    }


    void SimpleAccumulator::_reset()
    {   // reset must not fail on deallocated state
        if (_flag_alloc) { std::fill(_data, _data+_nobs, 0.); }
    }


    void SimpleAccumulator::_deallocate()
    {
        delete [] _data;
        _data = nullptr;
        _flag_alloc = false;
    }

}
