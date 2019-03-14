#ifndef MCI_UNIFORMALLMOVE_HPP
#define MCI_UNIFORMALLMOVE_HPP

#include "mci/TrialMoveInterface.hpp"

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
        double _stepSize; // same step size for all indices
        std::uniform_real_distribution<double> _rd; // uniform random distribution

    protected:
        TrialMoveInterface * _clone() const override {
            return new UniformAllMove(this->getNDim(), 0.5*_stepSize);
        }

    public:
        UniformAllMove(int ndim, double initStepSize):
            TrialMoveInterface(ndim, 0), _stepSize(initStepSize),
            _rd(std::uniform_real_distribution<double>(-1.,1.)) // note that we chose the dist symmetric around 0
        {}

        // Methods required for auto-calibration:
        int getNStepSizes() const override { return 1; }
        void setStepSize(int, double val) override { _stepSize = val; }
        double getStepSize(int) const override { return _stepSize; }
        double getChangeRate() const override { return 1.; } // chance for a single index to change is 1 (because they all change)
        void getUsedStepSizes(int, const int[], int &nusedSizes, int usedSizeIdx[]) const override
        { // we have only one step size
            nusedSizes = 1;
            usedSizeIdx[0] = 0;
        }

        void protoFunction(const double[], double[]) override {} // not needed

        void onAcceptance(const SamplingFunctionContainer &, double[]) override {} // not needed

        double trialMove(double xnew[], int &nchanged, int[], const double[], double[]) override
        {
            // do step
            for (int i=0; i<this->getNDim(); ++i) {
                xnew[i] += _stepSize * _rd( *(this->getRGen()) );
            }
            nchanged = this->getNDim(); // if we changed all, we don't need to fill changedIdx

            return 1.; // uniform -> no move acceptance factor
        }
    };

} // namespace mci

#endif
