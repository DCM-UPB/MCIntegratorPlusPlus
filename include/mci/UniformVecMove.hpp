#ifndef MCI_UNIFORMVECMOVE_HPP
#define MCI_UNIFORMVECMOVE_HPP

#include "mci/TrialMoveInterface.hpp"

#include <algorithm>
#include <stdexcept>

namespace mci
{
    // Uniform single-vector move
    // Use this if your input x consists of vectors of length _veclen and you want to move a
    // random one of them on each step. Just use veclen=1, if you want single-index moves.
    // If you have different types of particles/vectors, and want to calibrate and use different
    // step size per particle type, you can the full constructor and pass ntypes>1. Vectors of a
    // common type are expected to be found in a single block within the position arrays. You need
    // to pass the typeEnds array containing the end-indices of every type block. We follow the
    // convention that "end" means one behind the last element.
    class UniformVecMove: public TrialMoveInterface
    {
    private:
        const int _nvecs; // how many vector/particles are considered (calculated as ndim/veclen)
        const int _veclen; // how many indices does one particle/vector have? (i.e. space dimension)
        const int _ntypes; // how many different types of particles do you have?
        int * const _typeEnds; // end-indices of every type in x (i.e. the last index of type i is _typeEnds[i]-1 )
        double * const _stepSizes; // holds the step sizes, one per type
        std::uniform_int_distribution<int> _rdidx; // uniform integer distribution to choose vector index
        std::uniform_real_distribution<double> _rdmov; // uniform double distribution to move vector

    protected:
        TrialMoveInterface * _clone() const override {
            return new UniformVecMove(_nvecs, _veclen, _ntypes, _typeEnds, _stepSizes);
        }

    public:
        // Full constructor, scalar step init
        UniformVecMove(int nvecs, int veclen, int ntypes, const int typeEnds[] /*len ntypes*/, double initStepSize /*scalar init*/):
            TrialMoveInterface(nvecs*veclen, 0), _nvecs(nvecs), _veclen(veclen), _ntypes(ntypes),
            _typeEnds(new int[_ntypes]), _stepSizes(new double[_ntypes]),
            _rdidx(std::uniform_int_distribution<int>(0, _nvecs-1)),
            _rdmov(std::uniform_real_distribution<double>(-1.,1.))
        {
            if (_nvecs < 1) { throw std::invalid_argument("[UniformVecMove] Number of vectors must be at least 1."); }
            if (_veclen < 1) { throw std::invalid_argument("[UniformVecMove] Vector length must be at least 1."); }
            if (_ntypes < 1) { throw std::invalid_argument("[UniformVecMove] Number of types must be at least 1."); }
            if (_ntypes > 1) {
                for (int i=0; i<_ntypes; ++i) { // we rely on this later
                    if (typeEnds[i]%_veclen != 0) {
                        throw std::invalid_argument("[UniformVecMove] All type end indices must be multiples of vector length.");
                    }
                }
                std::copy(typeEnds, typeEnds+_ntypes, _typeEnds);
            } else {
                _typeEnds[0] = _ndim;
            }
            std::fill(_stepSizes, _stepSizes+_ntypes, initStepSize);
        }

        // Full constructor, array step init
        UniformVecMove(int nvecs, int veclen, int ntypes, const int typeEnds[], const double initStepSizes[] /*len ntypes*/):
            UniformVecMove(nvecs, veclen, ntypes, typeEnds, 0.) // reuse above constructor
        {
            std::copy(initStepSizes, initStepSizes+_ntypes, _stepSizes); // put the proper values in
        }

        // ntype=1 constructor
        UniformVecMove(int nvecs, int veclen, double initStepSize):
            UniformVecMove(nvecs, veclen, 1, nullptr, initStepSize) // it is safe to use the constructor like this
        {}

        ~UniformVecMove() override {
            delete [] _stepSizes;
            delete [] _typeEnds;
        }

        // Methods required for auto-callibration:
        int getNStepSizes() const override { return _ntypes; }
        void setStepSize(int i, double val) override { _stepSizes[i] = val; }
        double getStepSize(int i) const override { return _stepSizes[i]; }
        double getChangeRate() const override { return 1./_nvecs; } // equivalent to _veclen/_ndim
        void getUsedStepSizes(const WalkerState &wlkstate, int &nusedSizes, int usedSizeIdx[]) const override
        { // we know that we changed only a single vector
            nusedSizes = 1;
            for (int i=0; i<_ntypes; ++i) {
                if (wlkstate.changedIdx[0] < _typeEnds[i]) {
                    usedSizeIdx[0] = i;
                    return;
                }
            }
            // this method is not performance-critical, so let's check for this
            throw std::runtime_error("[UniformVecMove::getUsedStepSizes] End of method reached, without result.");
        }

        void protoFunction(const double/*in*/[], double/*protovalues*/[]) override {} // not needed

        void onAcceptance(const SamplingFunctionContainer&/*pdfcont*/, double/*protoold*/[]) override {} // not needed

        double trialMove(WalkerState &wlkstate, const double/*protoold*/[], double/*protonew*/[]) override
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

} // namespace mci

#endif
