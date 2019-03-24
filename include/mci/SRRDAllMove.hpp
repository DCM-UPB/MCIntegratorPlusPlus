#ifndef MCI_SRRDALLMOVE_HPP
#define MCI_SRRDALLMOVE_HPP

#include "mci/TypedMoveInterface.hpp"

#include <algorithm>

namespace mci
{
// Generic all-particle move that can be instantiated for any type SSRD for which
// a specialization of the createSymNormalRRD function exists (e.g. in TrialMoveInterface.hpp).
// In any case, SRRD type must be callable and, when called with passed stdlib
// random generator, return a random double from a symmetric distribution around 0.
// By including this header you automatically "use" instantiations for types that
// have builtin support (see end of file).
//
// About the all-particle move:
// Probably the simplest kind of move, which requires the least effort in your
// sampling functions/observables, since you do not need to provide single-particle
// update methods for efficiency.
// However, if the number of particles is large, all-particle moves suffer from
// small step sizes and therefore lead to suboptimal auto-correlation times.
// Note: If you pass ntypes>1, follow the rules of TypedMoveInterface .
//
template <class SRRD /*symmetric, real-valued random distribution that works like standard library dists*/>
class SRRDAllMove: public TypedMoveInterface
{
private:
    SRRD _rd; // real-valued random distribution for move

    TrialMoveInterface * _clone() const final
    {
        return new SRRDAllMove(_ndim, _ntypes, _typeEnds, _stepSizes, &_rd);
    }

    // not used, make final for that extra performance
    void _newToOld() final {}
    void _oldToNew() final {}

public:
    // Full constructor with scalar step init and optionally passed pre-made random dist
    SRRDAllMove(int ndim, int ntypes, const int typeEnds[] /*len ntypes*/, double initStepSize /*scalar init*/, const SRRD * rdist = nullptr):
            TypedMoveInterface(ndim, 0, ntypes, typeEnds, initStepSize),
            _rd((rdist != nullptr) ? *rdist : createSymNormalRRD<SRRD>() /*fall-back*/ ) {}

    // Full constructor, with array step init and optionally passed pre-made random dist
    SRRDAllMove(int ndim, int ntypes, const int typeEnds[], const double initStepSizes[] /*len ntypes*/, const SRRD * rdist = nullptr):
            SRRDAllMove(ndim, ntypes, typeEnds, 0., rdist) // reuse above constructor
    {
        std::copy(initStepSizes, initStepSizes + _ntypes, _stepSizes); // put the proper values in
    }

    // ntype=1 constructor (i.e. scalar size), with optionally passed pre-made random dist
    SRRDAllMove(int ndim, double initStepSize, const SRRD * rdist = nullptr):
            SRRDAllMove(ndim, 1, nullptr, initStepSize, rdist) // it is safe to use the constructor like this
    {}

    // Method required for auto-calibration
    double getChangeRate() const final { return 1.; } // chance for a single index to change is 1 (because they all change)


    void protoFunction(const double/*in*/[], double/*protovalues*/[]) final {} // not needed

    double trialMove(WalkerState &wlk, const double/*protoold*/[], double/*protonew*/[]) final
    {
        // do step
        int xidx = 0;
        for (int tidx = 0; tidx < _ntypes; ++tidx) {
            while (xidx < _typeEnds[tidx]) {
                wlk.xnew[xidx] += _stepSizes[tidx]*_rd(*(_rgen));
                ++xidx;
            }
        }
        wlk.nchanged = _ndim; // if we changed all, we don't need to fill changedIdx

        return 1.; // symmetric distribution -> no move acceptance factor
    }
};

// Instantiations for applicable standard-library distributions
using UniformAllMove = SRRDAllMove<std::uniform_real_distribution<double>>;
using GaussianAllMove = SRRDAllMove<std::normal_distribution<double>>;
using StudentAllMove = SRRDAllMove<std::student_t_distribution<double>>;
using CauchyAllMove = SRRDAllMove<std::cauchy_distribution<double>>;

// the following ones use the symmetrized wrapper
using ExponentialAllMove = SRRDAllMove<SymmetrizedPRRD < std::exponential_distribution<double> > >;
using GammaAllMove = SRRDAllMove<SymmetrizedPRRD < std::gamma_distribution<double> > >;
using WeibullAllMove = SRRDAllMove<SymmetrizedPRRD < std::weibull_distribution<double> > >;
using LognormalAllMove = SRRDAllMove<SymmetrizedPRRD < std::lognormal_distribution<double> > >;
using ChisqAllMove = SRRDAllMove<SymmetrizedPRRD < std::chi_squared_distribution<double> > >;
using FisherAllMove = SRRDAllMove<SymmetrizedPRRD < std::fisher_f_distribution<double> > >;
} // namespace mci

#endif
