#include "mci/OrthoPeriodicDomain.hpp"

#include <algorithm>
#include <stdexcept>

namespace mci
{

    void OrthoPeriodicDomain::_checkBounds() const
    {
        for (int i=0; i<ndim; ++i) {
            if (ubounds[i]<=lbounds[i]) {
                throw std::invalid_argument("[OrthoPeriodicDomain::checkBounds] All upper bounds must be truly greater than their corresponding lower bounds.");
            }
        }
    }

    OrthoPeriodicDomain::OrthoPeriodicDomain(const int n_dim, const double l_bound, const double u_bound):
        DomainInterface(n_dim), lbounds(new double[n_dim]), ubounds(new double[n_dim])
    {
        std::fill(lbounds, lbounds+ndim, l_bound);
        std::fill(ubounds, ubounds+ndim, u_bound);
        this->_checkBounds();
    }

    OrthoPeriodicDomain::OrthoPeriodicDomain(const int n_dim, const double l_bounds[], const double u_bounds[]):
        DomainInterface(n_dim), lbounds(new double[n_dim]), ubounds(new double[n_dim])
    {
        std::copy(l_bounds, l_bounds+ndim, lbounds);
        std::copy(u_bounds, u_bounds+ndim, ubounds);
        this->_checkBounds();
    }

    OrthoPeriodicDomain::~OrthoPeriodicDomain()
    {
        delete [] ubounds;
        delete [] lbounds;
    }


    void OrthoPeriodicDomain::applyDomain(double x[]) const
    {
        for (int i=0; i<ndim; ++i) {
            while ( x[i] < lbounds[i] ) {
                x[i] += ubounds[i] - lbounds[i];
            }
            while ( x[i] > ubounds[i] ) {
                x[i] -= ubounds[i] - lbounds[i];
            }
        }
    }

    void OrthoPeriodicDomain::applyDomain(WalkerState &wlk) const
    {
        for (int i=0; i<wlk.nchanged; ++i) {
            const int idx = wlk.changedIdx[i];
            while ( wlk.xnew[idx] < lbounds[idx] ) {
                wlk.xnew[idx] += ubounds[idx] - lbounds[idx];
            }
            while ( wlk.xnew[idx] > ubounds[idx] ) {
                wlk.xnew[idx] -= ubounds[idx] - lbounds[idx];
            }
        }
    }

    void OrthoPeriodicDomain::scaleToDomain(double normX[]) const
    {
        for (int i=0; i<ndim; ++i) {
            normX[i] = lbounds[i] + normX[i] * (ubounds[i]-lbounds[i]);
        }
    }

    void OrthoPeriodicDomain::getSizes(double dimSizes[]) const
    {
        for (int i=0; i<ndim; ++i) {
            dimSizes[i] = ubounds[i] - lbounds[i];
        }
    }

    double OrthoPeriodicDomain::getVolume() const
    {
        double vol=1.;
        for (int i=0; i<ndim; ++i) {
            vol *= ( ubounds[i] - lbounds[i] );
        }
        return vol;
    }

} // namespace mci
