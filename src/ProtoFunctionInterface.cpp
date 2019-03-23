#include "mci/ProtoFunctionInterface.hpp"

#include <algorithm>
#include <stdexcept>

namespace mci
{

    ProtoFunctionInterface::ProtoFunctionInterface(const int ndim, const int nproto):
        _ndim(ndim), _nproto(0), _protoold(nullptr), _protonew(nullptr)
    {
        if (ndim < 1) { throw std::invalid_argument("[ProtoFunctionInterface] Number of dimensions must be at least 1."); }
        this->setNProto(nproto);
    }

    ProtoFunctionInterface::~ProtoFunctionInterface()
    {
        delete [] _protonew;
        delete [] _protoold;
    }

    void ProtoFunctionInterface::setNProto(const int nproto)
    {
        delete [] _protonew;
        delete [] _protoold;
        if (nproto>0) {
            _protoold = new double[nproto];
            _protonew = new double[nproto];
            std::fill(_protoold, _protoold+nproto, 0.);
            std::fill(_protonew, _protonew+nproto, 0.);
            _nproto=nproto;
        }
        else {
            _protoold = nullptr;
            _protonew = nullptr;
            _nproto = 0;
        }
    }

    void ProtoFunctionInterface::initializeProtoValues(const WalkerState &wlk) {
        this->protoFunction(wlk.xold, _protonew);
        this->newToOld(wlk);
    }

    void ProtoFunctionInterface::newToOld(const WalkerState &wlk)
    {   // copy new values to old
        std::copy(_protonew, _protonew+_nproto, _protoold);
        this->_newToOld(wlk);
    }

    void ProtoFunctionInterface::oldToNew()
    {   // copy old values to new
        std::copy(_protoold, _protoold+_nproto, _protonew);
        this->_oldToNew();
    }

}  // namespace mci
