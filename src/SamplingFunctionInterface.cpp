#include "mci/SamplingFunctionInterface.hpp"

#include <algorithm>

namespace mci
{

    SamplingFunctionInterface::SamplingFunctionInterface(const int ndim, const int nproto):
        _ndim(ndim), _nproto(nproto)
    {
        _protonew = new double[_nproto];
        _protoold = new double[_nproto];
        std::fill(_protonew, _protonew+_nproto, 0.);
        std::fill(_protoold, _protoold+_nproto, 0.);
    }

    SamplingFunctionInterface::~SamplingFunctionInterface()
    {
        delete[] _protoold;
        delete[] _protonew;
    }

    void SamplingFunctionInterface::setNProto(const int nproto)
    {
        delete[] _protoold;
        delete[] _protonew;
        _protonew = new double[_nproto];
        _protoold = new double[_nproto];
        std::fill(_protonew, _protonew+_nproto, 0.);
        std::fill(_protoold, _protoold+_nproto, 0.);
        _nproto=nproto;
    }

    void SamplingFunctionInterface::newToOld()
    {   // pointer swap
        std::swap(_protonew, _protoold);
    }

}  // namespace mci
