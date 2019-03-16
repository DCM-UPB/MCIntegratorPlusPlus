#ifndef MCI_SRRDVECMOVE_HPP
#define MCI_SRRDVECMOVE_HPP

#include "mci/TypedMoveInterface.hpp"

#include <algorithm>
#include <stdexcept>

namespace mci
{
    // Generic single-vector move (instantiations for standard distributions at the end of file)
    //
    // Use this if your input x consists of vectors of length _veclen and you want to move a
    // random one of them on each step. Just use veclen=1, if you want single-index moves.
    // Note: If you use ntypes>1, follow the instructions of the TypedMoveInterface. Your
    // vectors must be aligned with the specified typeEnds (else constructor will throw).
    //
    template < class SRRD /*symmetric, real-valued random distribution that works like standard library dists*/>
    class SRRDVecMove: public TypedMoveInterface
    {
    private:
        const int _nvecs{}; // how many vector/particles are considered (calculated as ndim/veclen)
        const int _veclen{}; // how many indices does one particle/vector have? (i.e. space dimension)
        std::uniform_int_distribution<int> _rdidx; // uniform integer distribution to choose vector index
        SRRD _rdmov; // symmetric double-typed distribution to move vector

    protected:
        TrialMoveInterface * _clone() const final {
            return new SRRDVecMove(_nvecs, _veclen, _ntypes, _typeEnds, _stepSizes);
        }

    public:
        // Full constructor with scalar step init and optionally passed pre-made random dist
        SRRDVecMove(int nvecs, int veclen, int ntypes, const int typeEnds[] /*len ntypes*/, double initStepSize /*scalar*/, const SRRD * rdist = nullptr):
            TypedMoveInterface(nvecs*veclen, 0, ntypes, typeEnds, initStepSize), _nvecs(nvecs), _veclen(veclen),
            _rdidx(std::uniform_int_distribution<int>(0, _nvecs-1)),
            _rdmov( (rdist!=nullptr)? *rdist : createSymNormalRRD<SRRD>() /*fall-back*/ )
        {
            if (_nvecs < 1) { throw std::invalid_argument("[SRRDVecMove] Number of vectors must be at least 1."); }
            if (_veclen < 1) { throw std::invalid_argument("[SRRDVecMove] Vector length must be at least 1."); }
            if (_ntypes > 1) {
                for (int i=0; i<_ntypes; ++i) { // we rely on this later
                    if (_typeEnds[i]%_veclen != 0) {
                        throw std::invalid_argument("[SRRDVecMove] All type end indices must be multiples of vector length.");
                    }
                }
            }
        }

        // Full constructor, array step init
        SRRDVecMove(int nvecs, int veclen, int ntypes, const int typeEnds[], const double initStepSizes[] /*len ntypes*/, const SRRD * rdist = nullptr):
            SRRDVecMove(nvecs, veclen, ntypes, typeEnds, 0.) // reuse above constructor
        {
            std::copy(initStepSizes, initStepSizes+_ntypes, _stepSizes, rdist); // put the proper values in
        }

        // ntype=1 constructor
        SRRDVecMove(int nvecs, int veclen, double initStepSize, const SRRD * rdist = nullptr):
            SRRDVecMove(nvecs, veclen, 1, nullptr, initStepSize, rdist) // it is safe to use the constructor like this
        {}

        // Methods required for auto-calibration
        double getChangeRate() const final { return 1./_nvecs; } // equivalent to _veclen/_ndim
        void getUsedStepSizes(const WalkerState &wlkstate, int &nusedSizes, int usedSizeIdx[]) const final
        { // we know that we changed only a single vector
            nusedSizes = 1;
            for (int i=0; i<_ntypes; ++i) {
                if (wlkstate.changedIdx[0] < _typeEnds[i]) {
                    usedSizeIdx[0] = i;
                    return;
                }
            }
            // this method is not performance-critical, so let's check for this
            throw std::runtime_error("[SRRDVecMove::getUsedStepSizes] End of method reached, without result.");
        }

        void protoFunction(const double/*in*/[], double/*protovalues*/[]) final {} // not needed

        void onAcceptance(const SamplingFunctionContainer&/*pdfcont*/, double/*protoold*/[]) final {} // not needed

        double trialMove(WalkerState &wlkstate, const double/*protoold*/[], double/*protonew*/[]) final
        {
            // determine vector to change and its type
            const int vidx = _rdidx(*_rgen);
            const int xidx = vidx*_veclen; // first x index to change
            int tidx = 0; // type index
            while (tidx<_ntypes) { // This executes fast, when ntypes is small
                if (xidx < _typeEnds[tidx]) {
                    break;
                }
                ++tidx;
            }

            // do step
            for (int i=0; i<_veclen; ++i) {
                wlkstate.xnew[xidx + i] += _stepSizes[tidx] * _rdmov(*_rgen);
                wlkstate.changedIdx[i] = xidx + i;
            }
            wlkstate.nchanged = _veclen; // how many indices we changed

            return 1.; // uniform -> no move acceptance factor
        }
    };

    // Instantiations for applicable standard-library distributions
    using UniformVecMove = SRRDVecMove<std::uniform_real_distribution<double>>;
    using GaussianVecMove = SRRDVecMove<std::normal_distribution<double>>;
    using StudentVecMove = SRRDVecMove<std::student_t_distribution<double>>;
    using CauchyVecMove = SRRDVecMove<std::cauchy_distribution<double>>;

} // namespace mci

#endif
