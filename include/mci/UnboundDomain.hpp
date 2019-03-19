#ifndef MCI_UNBOUNDDOMAIN_HPP
#define MCI_UNBOUNDDOMAIN_HPP

#include "mci/DomainInterface.hpp"

#include <limits>
#include <algorithm>

namespace mci
{
    // Domain without boundaries, i.e. R^N
    struct UnboundDomain: public DomainInterface
    {
    protected:
        DomainInterface * _clone() const final {
            return new UnboundDomain(ndim);
        }

    public:
        UnboundDomain(int n_dim): DomainInterface(n_dim) {}

        // return somewhat safe value that should never be reached anyway
        double getMaxStepSize() const final { return 0.1*std::numeric_limits<double>::max(); }

        // these are trivial
        void getCenter(double centerX[]) const final { std::fill(centerX, centerX+ndim, 0.); }
        double getVolume() const final { return 0.; } // infinite volume -> convention is 0
        void applyDomain(double x[]) const final {} // do nothing
        void applyDomain(WalkerState &wlk) const final {} // still nothing
    };

} // namespace mci


#endif
