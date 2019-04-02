#ifndef MCI_UNBOUNDDOMAIN_HPP
#define MCI_UNBOUNDDOMAIN_HPP

#include "mci/DomainInterface.hpp"

#include <algorithm>
#include <limits>

namespace mci
{
// Domain without boundaries, i.e. R^N
struct UnboundDomain final: public DomainInterface
{
protected:
    DomainInterface * _clone() const final
    {
        return new UnboundDomain(ndim);
    }

public:
    explicit UnboundDomain(int n_dim): DomainInterface(n_dim) {}
    ~UnboundDomain() final = default;

    // most are trivial
    void applyDomain(double x[]) const final {} // do nothing
    void applyDomain(WalkerState &wlk) const final {} // still nothing

    void scaleToDomain(double normX[]) const final
    {
        for (int i = 0; i < ndim; ++i) {
            normX[i] = -domain_conv::infinity + normX[i]*domain_conv::infinityX2; // we use the infinity conventions
        }
    }

    void getSizes(double dimSizes[]) const final { std::fill(dimSizes, dimSizes + ndim, domain_conv::infinityX2); }
    double getVolume() const final { return domain_conv::infiniteVol; } // infinite volume convention
};
} // namespace mci


#endif
