#ifndef MCI_ORTHOPERIODICDOMAIN_HPP
#define MCI_ORTHOPERIODICDOMAIN_HPP

#include "mci/DomainInterface.hpp"

#include <limits>

namespace mci
{
    // Domain enforcing orthorhombic periodic boundary conditions.
    // The boundaries are defined by the two arrays _lbounds and _ubounds.
    struct OrthoPeriodicDomain: public DomainInterface
    {
    public:
        double * const lbounds; // lower boundaries
        double * const ubounds; // upper boundaries

    protected:
        DomainInterface * _clone() const final {
            return new OrthoPeriodicDomain(ndim, lbounds, ubounds);
        }

        void _checkBounds() const; // make sure the set bounds are reasonable

    public:
        OrthoPeriodicDomain(int n_dim, // use some huge default boundaries, but far less than the boundaries of doubles
                            double l_bound = -std::numeric_limits<float>::max(),
                            double u_bound =  std::numeric_limits<float>::max());

        OrthoPeriodicDomain(int n_dim, const double l_bounds[], const double u_bounds[]); // use arrays to set bounds

        ~OrthoPeriodicDomain() override;


        // return smallest box dimension
        double getMaxStepSize() const final;

        // center is in the middle between boundaries
        void getCenter(double centerX[]) const final;

        // volume is product of dimension lengths
        double getVolume() const final;

        // apply PBC to full x
        void applyDomain(double x[]) const final;

        // apply PBC to updated walkerstate
        void applyDomain(WalkerState &wlk) const final;
    };

} // namespace mci


#endif
