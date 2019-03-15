#include "mci/SamplingFunctionInterface.hpp"

#include <algorithm>

namespace mci
{
    double SamplingFunctionInterface::computeAcceptance(const double in[])
    {
        this->protoFunction(in, _protonew);
        return this->acceptanceFunction(_protoold, _protonew);
    }

    double SamplingFunctionInterface::computeAcceptance(const WalkerState &wlkstate)
    {
        if (wlkstate.nchanged < _ndim) {
            return this->updatedAcceptance(wlkstate, _protoold, _protonew);
        }
        // all elements have changed
        this->protoFunction(wlkstate.xnew, _protonew);
        return this->acceptanceFunction(_protoold, _protonew);
    }

}  // namespace mci
