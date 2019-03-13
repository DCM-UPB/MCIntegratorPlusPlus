#include "mci/SamplingFunctionInterface.hpp"

#include <algorithm>

namespace mci
{
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

}  // namespace mci
