#ifndef MCI_UNIFORMALLMOVE_HPP
#define MCI_UNIFORMALLMOVE_HPP

#include "mci/TrialMoveInterface.hpp"

#include <algorithm>

namespace mci
{
    // Uniform all-particle move. Probably the simplest kind of move,
    // which requires the least effort in your sampling functions/observables,
    // since you do not need to provide single-particle update methods for efficiency.
    // However, if the number of particles is large, all-particle moves suffer from
    // small step sizes and therefore lead to suboptimal auto-correlation times.
    class UniformAllMove: public TrialMoveInterface
    {
    private:
        const int _ntypes; // how many different types of particles do you have?
        int * const _typeEnds; // end-indices of every type in x (i.e. the last index of type i is _typeEnds[i]-1 )
        double * const _stepSizes; // holds the step sizes, one per type
        std::uniform_real_distribution<double> _rd; // uniform random distribution for move

    protected:
        TrialMoveInterface * _clone() const override {
            return new UniformAllMove(_ndim, _ntypes, _typeEnds, _stepSizes);
        }

    public:
        // Full constructor, scalar step init
        UniformAllMove(int ndim, int ntypes, const int typeEnds[] /*len ntypes*/, double initStepSize /*scalar init*/):
            TrialMoveInterface(ndim, 0), _ntypes(ntypes),
            _typeEnds(new int[_ntypes]), _stepSizes(new double[_ntypes]),
            _rd(std::uniform_real_distribution<double>(-1.,1.))
        {
            if (_ntypes < 1) { throw std::invalid_argument("[UniformAllMove] Number of types must be at least 1."); }
            if (_ntypes > 1) {
                std::copy(typeEnds, typeEnds+_ntypes, _typeEnds);
            } else {
                _typeEnds[0] = _ndim;
            }
            std::fill(_stepSizes, _stepSizes+_ntypes, initStepSize);
        }

        // Full constructor, array step init
        UniformAllMove(int ndim, int ntypes, const int typeEnds[], const double initStepSizes[] /*len ntypes*/):
            UniformAllMove(ndim, ntypes, typeEnds, 0.) // reuse above constructor
        {
            std::copy(initStepSizes, initStepSizes+_ntypes, _stepSizes); // put the proper values in
        }

        // ntype=1 constructor
        UniformAllMove(int ndim, double initStepSize):
            UniformAllMove(ndim, 1, nullptr, initStepSize) // it is safe to use the constructor like this
        {}

        ~UniformAllMove() override {
            delete [] _stepSizes;
            delete [] _typeEnds;
        }

        // Methods required for auto-calibration:
        int getNStepSizes() const override { return _ntypes; }
        void setStepSize(int i, double val) override { _stepSizes[i] = val; }
        double getStepSize(int i) const override { return _stepSizes[i]; }
        double getChangeRate() const override { return 1.; } // chance for a single index to change is 1 (because they all change)
        void getUsedStepSizes(int /*nchangedX*/, const int /*changedIdx*/[], int &nusedSizes, int usedSizeIdx[]) const override
        { // we always use all step sizes
            nusedSizes = _ntypes;
            std::iota(usedSizeIdx, usedSizeIdx+_ntypes, 0); // fill 0..._ntypes-1
        }

        void protoFunction(const double/*in*/[], double/*protovalues*/[]) override {} // not needed

        void onAcceptance(const SamplingFunctionContainer&/*pdfcont*/, double/*protoold*/[]) override {} // not needed

        double trialMove(double xnew[], int &nchanged, int/*changedIdx*/[], const double/*protoold*/[], double/*protonew*/[]) override
        {
            // do step
            int xidx = 0;
            for (int tidx=0; tidx<_ntypes; ++tidx) {
                while (xidx < _typeEnds[tidx]) {
                    xnew[xidx] += _stepSizes[tidx] * _rd( *(_rgen) );
                    ++xidx;
                }
            }
            nchanged = _ndim; // if we changed all, we don't need to fill changedIdx

            return 1.; // uniform -> no move acceptance factor
        }
    };

} // namespace mci

#endif

