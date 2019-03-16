#ifndef MCI_TRIALMOVEINTERFACE_HPP
#define MCI_TRIALMOVEINTERFACE_HPP

#include "mci/Clonable.hpp"
#include "mci/ProtoFunctionInterface.hpp"
#include "mci/SamplingFunctionContainer.hpp"
#include "mci/WalkerState.hpp"

#include <random>

namespace mci
{
    // Base interface for trial moves
    //
    // For details see the methods&comments below.
    //
    class TrialMoveInterface: public ProtoFunctionInterface, public Clonable<TrialMoveInterface>
    {
    protected:
        // we share the random generator with MCI, get's passed on bindRGen()
        std::mt19937_64 * _rgen; // ptr to MCI's rgen
        std::mt19937_64 * getRGen() const { return _rgen; }

    public:
        TrialMoveInterface(int ndim, int nproto): ProtoFunctionInterface(ndim, nproto), _rgen(nullptr) {}

        // store ptr to MCI's rgen
        void bindRGen(std::mt19937_64 &rgen) {
            _rgen = &rgen;
        }

        // Called after step accepted (after newToOld)
        void callOnAcceptance(const SamplingFunctionContainer &pdfcont) { this->onAcceptance(pdfcont, _protoold); }

        // compute move, for details see below
        double computeTrialMove(WalkerState &wlk) { return this->trialMove(wlk, _protoold, _protonew); }

        // do we have step sizes to calibrate?
        bool hasStepSizes() const { return (this->getNStepSizes() > 0); }

        // --- METHODS THAT MUST BE IMPLEMENTED

        // Methods only used for auto-calibration of step sizes:
        virtual int getNStepSizes() const = 0; // how many step-size like parameters do you have? (usually 0, 1 or ndim)
        virtual double getStepSize(int i) const = 0; // get step size with index i
        virtual void setStepSize(int i, double val) = 0; // set step size with index i to val
        virtual double getChangeRate() const = 0; // average probability of a single x index to change on move (e.g. 1./ndim on single-index moves)
        // tell which step sizes were used in a step described by wlk
        virtual void getUsedStepSizes(const WalkerState &wlk, int &nusedSizes, int usedSizeIdx[]) const = 0;

        void scaleStepSize(int i, double fac) { this->setStepSize(i, this->getStepSize(i)*fac); } // scale step size with index i by fac

        // Methods used during sampling:

        // Proto-value function
        // Remember to implement the protoFunction from ProtoFunctionInterface, if you can optimize
        // calculation of your trial move acceptance factor, by storing two sets of temporaries (new/old).
        // For details, see ProtoFunctionInterface.hpp. If you don't need that, in your constructor call
        // the TrialMoveInterface constructor with nproto set to 0 and and implement protoFunction as:
        //     double protoFunction(const double[], double[]) const final {}

        // Callback on acceptance
        // If your trial move requires knowledge of the previous sampling function value,
        // you may extract it on acceptance and store it into a designated proto value element.
        // Otherwise, leave the method empty.
        // Note: This method will also be called before the first step, to initialize.
        virtual void onAcceptance(const SamplingFunctionContainer &pdfcont, double protoold[]) = 0; // e.g. protoold[0] = pdfcont.getOldSamplingFunction();

        // Propose a trial move
        // For performance reasons, this method does everything related to a new trial move, at once.
        // We pass the walker state to update (i.e. xnew, nchanged and changedIdx) and your protovalues.
        // Update xnew (and your protonew), count changed indices in nchanged, store the indices in
        // ascending order in changedIdx (is allocated to length ndim) and finally return the acceptance
        // factor of your trial move.
        virtual double trialMove(WalkerState &wlk, const double protoold[], double protonew[]) = 0;
    };


    // Useful helper for generic derived classes
    // A helper template used to default-initialize real random distributions from the standard library,
    // to be symmetric around 0 and, if possible, standard deviation 1 (except: interval dists are [-1,1]).
    // The template is made to compile only for the standard-library distributions that allow to have
    // these properties. If such distributions are used for MC moves, the detailed balance condition
    // is fullfilled automatically.
    template < class SRRD >/* should be applicable stdlib random dist <double> */
    auto createSymNormalRRD() = delete; // cannot be used without specialization

    // Specialization for uniform
    template<>
    inline auto createSymNormalRRD< std::uniform_real_distribution<double> >() {
        return std::uniform_real_distribution<double>(-1., 1.); // the odd case with sigma!=1
    }

    // Specialization for gaussian
    template<>
    inline auto createSymNormalRRD< std::normal_distribution<double> >() {
        return std::normal_distribution<double>(); // defaults to mean 0 std 1
    }

    // Specialization for student-t
    template<>
    inline auto createSymNormalRRD< std::student_t_distribution<double> >() {
        return std::student_t_distribution<double>(); // default is symmetric around 0, but stddev undefined
    }

    // Specialization for cauchy
    template<>
    inline auto createSymNormalRRD< std::cauchy_distribution<double> >() {
        return std::cauchy_distribution<double>(); // default is symmetric around 0, but mean&stddev undefined
    }

}; // namespace mci

#endif
