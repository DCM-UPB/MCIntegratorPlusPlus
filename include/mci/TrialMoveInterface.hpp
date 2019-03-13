#ifndef MCI_TRIALMOVEINTERFACE_HPP
#define MCI_TRIALMOVEINTERFACE_HPP

#include "mci/ProtoFunctionInterface.hpp"
#include "mci/Clonable.hpp"
#include "mci/SamplingFunctionContainer.hpp"

namespace mci
{

    class TrialMoveInterface: public ProtoFunctionInterface, public Clonable<TrialMoveInterface>
    {
    protected:
        // we share the random generator with MCI, get's passed on bindRGen()
        std::mt19937_64 * _rgen; // ptr to MCI's rgen
        std::mt19937_64 * getRGen() { return _rgen; }

    public:
        TrialMoveInterface(int ndim, int nproto): ProtoFunctionInterface(ndim, nproto), _rgen(nullptr) {}

        // store ptr to MCI's rgen
        void bindRGen(const std::mt19937_64 &rgen) {
            _rgen = &rgen;
        }

        // Called after step accepted (after newToOld)
        void callOnAcceptance(const SamplingFunctionContainer &pdfcont) { this->onAcceptance(pdfcont, _protoold); }

        // compute move, for details see below
        double computeTrialMove(double xnew[], int &nchanged, int changedIdx[]) { return this->trialMove(xnew, nchanged, changedIdx, _protoold, _protonew); }

        // do we have step sizes to calibrate?
        bool hasStepSizes() const { return (this->getNStepSizes() > 0); }

        // --- METHODS THAT MUST BE IMPLEMENTED

        // Methods required for auto-callibration:
        virtual int getNStepSizes() const = 0; // how many step-size like parameters do you have? (typically 0, 1 or ndim)
        virtual void setStepSize(int i, double val) = 0; // set step size with index i to val

        // Remember to implement the protoFunction from ProtoFunctionInterface, if you can optimize
        // your trial moves by storing two sets of temporary values (new/old). For details, see
        // ProtoFunctionInterface.hpp. If you don't need that, in your constructor call the Trial-
        // MoveInterface constructor with nproto set to 0 and and implement protoFunction as:
        //     double protoFunction(const double[], double[]) const override {}

        // If your trial move requires knowledge of the previous sampling function value,
        // you may extract it on acceptance and store it into a designated proto value element.
        // Otherwise, leave the method empty.
        // Note: This method will also be called before the first step, to initialize.
        virtual void onAcceptance(const SamplingFunctionContainer &pdfcont, double protoold[]) = 0; // e.g. protoold[0] = pdfcont.getOldSamplingFunction();

        // Propose a trial move
        // For performance reasons, this method does everything related to a new trial move, at once.
        // That means: Update xnew (and protonew), count changed indices in nchanged, store the indices in
        // ascending order in changedIdx (is allocated to length ndim) and finally return the acceptance
        // factor of your trial move.
        // Note: The passed xnew array is in/out and originally contains the previous positions.
        virtual double trialMove(double xnew[], int &nchanged, int changedIdx[], const double protoold[], double protonew[]) = 0;
    };

}; // namespace mci

#endif
