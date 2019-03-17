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
        double * const xold; // ptr to old positions
        double * const xnew; // ptr to new positions

        int nchanged{}; // number of differing indices between xold and xnew
        /* NOTE: If nchanged=ndim, changedIdx is allowed to be invalid! You must check for that case (i.e. all-particle moves)!!) */
        int * const changedIdx; // first nchanged elements are the differing indices, in order
        bool accepted{}; // is the step accepted?

        explicit WalkerState(int n_dim): // initialize
            ndim(n_dim), xold(new double[ndim]),
            xnew(new double[ndim]), changedIdx(new int[ndim])
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
            // prepare sampling run (initial state is "accepted")
            std::copy(xold, xold+ndim, xnew);
            nchanged = ndim;
            std::iota(changedIdx, changedIdx+ndim, 0); // fill 0..ndim-1
            accepted = true;
        }

        void newToOld() { std::copy(xnew, xnew+ndim, xold); } // on acceptance
        void oldToNew() { std::copy(xold, xold+ndim, xnew); } // on rejection
    };

} // namespace mci

#endif
