#ifndef MCI_DOMAININTERFACE_HPP
#define MCI_DOMAININTERFACE_HPP

#include "mci/Clonable.hpp"
#include "mci/WalkerState.hpp"

#include <stdexcept>

namespace mci
{
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

        DomainInterface(int n_dim): ndim(n_dim) {
            if (ndim < 1) {
                throw std::invalid_argument("[DomainInterface] Number of dimensions must be at least 1.");
            }
        }
        virtual ~DomainInterface() = default;

        bool isFinite() const { return (this->getVolume() != 0.); } // convention is infinite box -> vol = 0

        // METHODS TO BE IMPLEMENTED
        // NOTE: all exchanged arrays have length ndim
        virtual double getMaxStepSize() const = 0; // return the maximal reasonable scalar step size (e.g. length of shortest dimension)
        virtual void getCenter(double centerX[]) const = 0; // writes center of domain into centerX
        virtual double getVolume() const = 0; // volume of the domain (if infinite, return 0)
        virtual void applyDomain(double x[]) const = 0; // apply domain to passed vector
        virtual void applyDomain(WalkerState &wlk) const = 0; // apply domain (selectively) to walker state's xnew
    };

} // namespace mci


#endif
