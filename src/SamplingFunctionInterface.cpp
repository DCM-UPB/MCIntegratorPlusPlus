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

    double SamplingFunctionInterface::computeAcceptance(const double in[])
    {
        this->protoFunction(in, _protonew);
        return this->acceptanceFunction(_protoold, _protonew);
    }

    double SamplingFunctionInterface::computeAcceptance(const double xold[], const double xnew[], int nchanged, const int changedIdx[])
    {
        if (nchanged < _ndim) {
            return this->updatedAcceptance(xold, xnew, nchanged, changedIdx, _protoold, _protonew);
        } else { // all elements have changed
            this->protoFunction(xnew, _protonew);
            return this->acceptanceFunction(_protoold, _protonew);
        }
    }


    void SamplingFunctionInterface::newToOld()
    {   // copy new values to old
        std::copy(_protonew, _protonew+_nproto, _protoold);
        this->_newToOld();
    }

}  // namespace mci
