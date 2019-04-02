#ifndef MCI_ORTHOPERIODICDOMAIN_HPP
#define MCI_ORTHOPERIODICDOMAIN_HPP

#include "mci/DomainInterface.hpp"

#include <limits>

namespace mci
{
// Domain enforcing orthorhombic periodic boundary conditions.
// The boundaries are defined by the two arrays _lbounds and _ubounds.
struct OrthoPeriodicDomain final: public DomainInterface
{
public:
    double * const lbounds; // lower boundaries
    double * const ubounds; // upper boundaries

protected:
    DomainInterface * _clone() const final
    {
        return new OrthoPeriodicDomain(ndim, lbounds, ubounds);
    }

    void _checkBounds() const; // make sure the set bounds are reasonable

public:
    explicit OrthoPeriodicDomain(int n_dim, // use the infinity conventions
                                 double l_bound = -domain_conv::infinity,
                                 double u_bound = domain_conv::infinity);

    OrthoPeriodicDomain(int n_dim, const double l_bounds[], const double u_bounds[]); // use arrays to set bounds
    ~OrthoPeriodicDomain() final;

    // apply PBC to full x
    void applyDomain(double x[]) const final;

    // apply PBC to updated walkerstate
    void applyDomain(WalkerState &wlk) const final;

    // transform normX in (0,1)^N to true box coordinates
    void scaleToDomain(double normX[]) const final;

    // fill with ubound - lbound
    void getSizes(double dimSizes[]) const final;

    // volume is product of dimension lengths
    double getVolume() const final;
};
} // namespace mci


#endif
