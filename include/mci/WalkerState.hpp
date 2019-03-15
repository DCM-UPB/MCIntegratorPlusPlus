#ifndef MCI_WALKERSTATE_HPP
#define MCI_WALKERSTATE_HPP

#include <algorithm>
#include <numeric>

namespace mci
{

    // Holds together information about the current walker move,
    // which needs to be passed around a lot. Passing the struct
    // instead of multiple parameters reduces function call overhead.
    struct WalkerState
    {
        const int ndim;
        double * const xold;
        double * const xnew;
        int nchanged;
        int * const changedIdx;

        explicit WalkerState(int n_dim):
            ndim(n_dim), xold(new double[ndim]), xnew(new double[ndim]),
            nchanged(ndim), changedIdx(new int[ndim])
        {
            std::fill(xold, xold+ndim, 0.);
            this->initialize();
        }

        ~WalkerState() {
            delete [] changedIdx;
            delete [] xnew;
            delete [] xold;
        }

        void initialize() {
            // prepare sampling run
            std::copy(xold, xold+ndim, xnew);
            nchanged = ndim;
            std::iota(changedIdx, changedIdx+ndim, 0); // fill 0..ndim-1
        }

        void newToOld() { std::copy(xnew, xnew+ndim, xold); } // on acceptance
        void oldToNew() { std::copy(xold, xold+ndim, xnew); } // on rejection
    };

} // namespace mci

#endif
