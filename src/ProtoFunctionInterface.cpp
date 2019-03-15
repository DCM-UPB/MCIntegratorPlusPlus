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
        delete[] _protoold; // this one was used to allocate
    }

    void ProtoFunctionInterface::setNProto(const int nproto)
    {
        delete[] _protoold;
        if (nproto>0) {

            _protoold = new double[2*nproto];
            _protonew = _protoold + nproto; // ptr to second half
            std::fill(_protoold, _protoold+2*nproto, 0.);
            _nproto=nproto;
        }
        else {
            _protonew = nullptr;
            _protoold = nullptr;
            _nproto = 0;
        }
    }

    void ProtoFunctionInterface::initializeProtoValues(const double in[]) {
        this->protoFunction(in, _protonew);
        this->newToOld();
    }

    void ProtoFunctionInterface::newToOld()
    {   // copy new values to old
        std::copy(_protonew, _protonew+_nproto, _protoold);
        this->_newToOld();
    }

    void ProtoFunctionInterface::oldToNew()
    {   // copy old values to new
        std::copy(_protoold, _protoold+_nproto, _protonew);
        this->_oldToNew();
    }

}  // namespace mci
