#include "mci/ProtoFunctionInterface.hpp"

#include <algorithm>

namespace mci
{

    ProtoFunctionInterface::ProtoFunctionInterface(const int ndim, const int nproto):
        _ndim(ndim), _nproto(0), _protonew(nullptr), _protoold(nullptr)
    {
        this->setNProto(nproto);
    }

    ProtoFunctionInterface::~ProtoFunctionInterface()
    {
        delete[] _protoold;
        delete[] _protonew;
    }

    void ProtoFunctionInterface::setNProto(const int nproto)
    {
        delete[] _protoold;
        delete[] _protonew;
        if (nproto>0) {
            _protonew = new double[_nproto];
            _protoold = new double[_nproto];
            std::fill(_protonew, _protonew+_nproto, 0.);
            std::fill(_protoold, _protoold+_nproto, 0.);
            _nproto=nproto;
        }
        else {
            _protonew = nullptr;
            _protoold = nullptr;
            _nproto = 0;
        }
    }

    void ProtoFunctionInterface::computeOldProtoValues(const double in[]) {
        this->protoFunction(in, _protonew);
        this->newToOld();
    }

    void ProtoFunctionInterface::newToOld()
    {   // copy new values to old
        std::copy(_protonew, _protonew+_nproto, _protoold);
        this->_newToOld();
    }

}  // namespace mci
