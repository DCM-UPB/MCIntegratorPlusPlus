#ifndef MCI_DOMAININTERFACE_HPP
#define MCI_DOMAININTERFACE_HPP

#include "mci/Clonable.hpp"
#include "mci/WalkerState.hpp"

#include <stdexcept>

namespace mci
{
namespace domain_conv
{ // conventions on what to use when you want infinite boundaries
static constexpr double infinity = std::numeric_limits<float>::max(); // 3.402823e+38
static constexpr double infinityX2 = infinity + infinity; // length from -infinity to +infinity
static constexpr double infiniteVol = 0.; // convention is infinite box -> vol = 0
} // namespace domain_conv

// Derived classes of DomainInterface define integration domains for
// MCI to integrate over. Note that distance calculation in the form
// of minimum image convention is intentionally left out of the inter-
// face, because it has no relevance for MCI itself. If you need such
// functionality, create a new interface derived from DomainInterface,
// add your extra methods there and derive your final Domain classes
// from that new interface. You can still pass these to MCI.
struct DomainInterface: public Clonable<DomainInterface>
{
public:
    const int ndim; // total dimension to integrate over (usually spacedim*nparticles)

    explicit DomainInterface(int n_dim): ndim(n_dim)
    {
        if (ndim < 1) {
            throw std::invalid_argument("[DomainInterface] Number of dimensions must be at least 1.");
        }
    }

    // check if the domain can be used without sampling function
    bool isFinite() const { return (this->getVolume() != domain_conv::infiniteVol); } // see above namespace with conventions

    // this may be useful
    void getCenter(double centerX[]) const
    {
        std::fill(centerX, centerX + ndim, 0.5);
        this->scaleToDomain(centerX);
    }

    // METHODS TO BE IMPLEMENTED
    // NOTE: all exchanged arrays have length ndim
    virtual void applyDomain(double x[]) const = 0; // apply domain to passed vector
    virtual void applyDomain(WalkerState &wlk) const = 0; // apply domain (selectively) to walker state's xnew
    virtual void scaleToDomain(double normX[]/*inout*/) const = 0; // transform positions normX (in (0,1)^N) to actual positions in domain
    virtual void getSizes(double dimSizes[]) const = 0; // lengths of the ndim dimensions (if infinite, use domain_conv::infinityX2)
    virtual double getVolume() const = 0; // volume of the domain (if infinite, return domain_conv::infiniteVol)
};
} // namespace mci


#endif
