#ifndef MCI_UNIFORMNSIZEMOVE_HPP
#define MCI_UNIFORMNSIZEMOVE_HPP

#include "mci/TrialMoveInterface.hpp"

#include <algorithm>

namespace mci
{

    class UniformNSizeMove: public TrialMoveInterface
    {
    private:
        std::uniform_real_distribution<double> _rd; // uniform random distribution
        double * const _stepSizes; // holds the step sizes, one per dimension

    public:
        UniformTrialMove(int ndim, double initStepSize = 0.1):
            ProtoFunctionInterface(ndim, 1), _stepSizes(new double[ndim]),
            _rd(std::uniform_real_distribution<double>(-1.,1.))
        {
            std::fill(_stepSizes, _stepSizes+ndim, initStepSize);
        }

        ~UniformNSizeMove() { delete [] _stepSizes; }

        // Methods required for auto-callibration:
        int getNStepSizes() const override { return this->getNDim(); }
        void setStepSize(int i, double val) override { _stepSizes[i] = val; };

        void onAcceptance(const SamplingFunctionContainer &pdfcont, double protoold[]) override {} // not needed

        double trialMove(double xnew[], int &nchanged, int changedIdx[], const double protoold[], double protonew[]) override
        {
            for (int i=0; i<this->getNDim(); ++i) {
                xnew[i] += _stepSizes[i] * _rd( *(this->getRGen()) );
            }
            nchanged = this->getNDim(); // if we say this we don't need to set changedIdx

            return 1.; // uniform -> no move acceptance factor
        }

    }; // namespace mci

#endif
