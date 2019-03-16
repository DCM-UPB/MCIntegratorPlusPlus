#ifndef MCI_SRRDALLMOVE_HPP
#define MCI_SRRDALLMOVE_HPP

#include "mci/TypedMoveInterface.hpp"

#include <algorithm>

namespace mci
{
    // Generic all-particle move (instantiations for standard distributions at the end of file)
    //
    // Probably the simplest kind of move, which requires the least effort in your
    // sampling functions/observables, since you do not need to provide single-particle
    // update methods for efficiency.
    // However, if the number of particles is large, all-particle moves suffer from
    // small step sizes and therefore lead to suboptimal auto-correlation times.
    // Note: If you pass ntypes>1, follow the rules of TypedMoveInterface .
    //
    template < class SRRD /*symmetric, real-valued random distribution that works like standard library dists*/>
    class SRRDAllMove: public TypedMoveInterface
    {
    private:
        SRRD _rd; // real-valued random distribution for move

        TrialMoveInterface * _clone() const final {
            return new SRRDAllMove(_ndim, _ntypes, _typeEnds, _stepSizes, &_rd);
        }

    public:
        // Full constructor with scalar step init and optionally passed pre-made random dist
        SRRDAllMove(int ndim, int ntypes, const int typeEnds[] /*len ntypes*/, double initStepSize /*scalar init*/, const SRRD * rdist = nullptr):
            TypedMoveInterface(ndim, 0, ntypes, typeEnds, initStepSize),
            _rd( (rdist!=nullptr)? *rdist : createSymNormalRRD<SRRD>() /*fall-back*/ )
        {}

        // Full constructor, with array step init and optionally passed pre-made random dist
        SRRDAllMove(int ndim, int ntypes, const int typeEnds[], const double initStepSizes[] /*len ntypes*/, const SRRD * rdist = nullptr):
            SRRDAllMove(ndim, ntypes, typeEnds, 0., rdist) // reuse above constructor
        {
            std::copy(initStepSizes, initStepSizes+_ntypes, _stepSizes); // put the proper values in
        }

        // ntype=1 constructor (i.e. scalar size), with optionally passed pre-made random dist
        SRRDAllMove(int ndim, double initStepSize, const SRRD * rdist = nullptr):
            SRRDAllMove(ndim, 1, nullptr, initStepSize, rdist) // it is safe to use the constructor like this
        {}

        // Methods required for auto-calibration
        double getChangeRate() const final { return 1.; } // chance for a single index to change is 1 (because they all change)
        void getUsedStepSizes(const WalkerState&/*wlkstate*/, int &nusedSizes, int usedSizeIdx[]) const final
        { // we always use all step sizes
            nusedSizes = _ntypes;
            std::iota(usedSizeIdx, usedSizeIdx+_ntypes, 0); // fill 0..._ntypes-1
        }

        void protoFunction(const double/*in*/[], double/*protovalues*/[]) final {} // not needed

        void onAcceptance(const SamplingFunctionContainer&/*pdfcont*/, double/*protoold*/[]) final {} // not needed

        double trialMove(WalkerState &wlkstate, const double/*protoold*/[], double/*protonew*/[]) final
        {
            // do step
            int xidx = 0;
            for (int tidx=0; tidx<_ntypes; ++tidx) {
                while (xidx < _typeEnds[tidx]) {
                    wlkstate.xnew[xidx] += _stepSizes[tidx] * _rd( *(_rgen) );
                    ++xidx;
                }
            }
            wlkstate.nchanged = _ndim; // if we changed all, we don't need to fill changedIdx

            return 1.; // uniform -> no move acceptance factor
        }
    };

    // Instantiations for applicable standard-library distributions
    using UniformAllMove = SRRDAllMove<std::uniform_real_distribution<double>>;
    using GaussianAllMove = SRRDAllMove<std::normal_distribution<double>>;
    using StudentAllMove = SRRDAllMove<std::student_t_distribution<double>>;
    using CauchyAllMove = SRRDAllMove<std::cauchy_distribution<double>>;

} // namespace mci

#endif