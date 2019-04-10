#ifndef MCI_TRIALMOVEINTERFACE_HPP
#define MCI_TRIALMOVEINTERFACE_HPP

#include "mci/Clonable.hpp"
#include "mci/ProtoFunctionInterface.hpp"
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
    void bindRGen(std::mt19937_64 &rgen)
    {
        _rgen = &rgen;
    }

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
    virtual int getStepSizeIndex(int xidx) const = 0; // return the step size index that corresponds to given position index xidx

    void scaleStepSize(int i, double fac) { this->setStepSize(i, this->getStepSize(i)*fac); } // scale step size with index i by fac
    void scaleStepSizes(double fac)
    { // scale all step sizes by fac
        for (int i = 0; i < this->getNStepSizes(); ++i) { this->scaleStepSize(i, fac); }
    }

    // Methods used during sampling:

    // Proto-value function
    // Remember to implement the protoFunction from ProtoFunctionInterface, if you can optimize
    // calculation of your trial move acceptance factor, by storing two sets of temporaries (new/old).
    // For details, see ProtoFunctionInterface.hpp. If you don't need that, in your constructor call
    // the TrialMoveInterface constructor with nproto set to 0 and and implement protoFunction as:
    //     void protoFunction(const double[], double[]) final {}


    // Propose a trial move
    // For performance reasons, this method does everything related to a new trial move, at once.
    // We pass the walker state to update (i.e. xnew, nchanged and changedIdx) and your protovalues.
    // Update xnew (and your protonew), count changed indices in nchanged, store the indices in
    // ascending order in changedIdx (is allocated to length ndim) and finally return the acceptance
    // factor of your trial move.
    virtual double trialMove(WalkerState &wlk, const double protoold[], double protonew[]) = 0;
};


// --- Useful stuff for derived classes
// NOTE: Below we consider std::mt19937_64 as the only possible random generator.


// Template to turn non-symmetric distributions generating x>0 into
// distributions which are symmetric around 0. This is done by first
// generating a x>0 from the given distribution, and then decide on
// the sign by a draw from bernoulli-distribution.
template <class PRRD> /* should be positive-real-valued stdlib random dist <double>*/
struct SymmetrizedPRRD
{
private:
    std::bernoulli_distribution bd; // bernoulli dist, yielding true and false with equal probability (default)

public:
    PRRD prrd; // externally passed positive-real-valued random distribution

    SymmetrizedPRRD() = default;
    explicit SymmetrizedPRRD(PRRD myPRRD) { prrd = myPRRD; }

    double operator()(std::mt19937_64 &rgen)
    { // this overload allows it to be used like other random distributions
        const double val = prrd(rgen);
        return (bd(rgen)) ? val : -val; // return + or - val, with same probability
    }
};


// Default initialization template to obtain symmetric distributions
// A helper template used to default-initialize real random distributions,
// e.g. (symmetric/uniform) distributions from the standard library or
// others when wrapped with SymmetrizedPRRD(). At the moment the helper
// is only really needed to make an exception for the uniform distribution,
// which is not symmetric around 0 by default.
// Remember that distributions used for simple MC moves should fulfill this
// condition in order to maintain detailed balance (automatically).
//
// NOTE: All applicable standard distributions have explicit specializations here!
//       You may still use different types if you know they are applicable,
//       as long as their default constructed version is symmetric around 0.
//
template <class SRRD /* applicable random dist */>
inline auto createSymNormalRRD() // default specialization
{
    return SRRD(); // return default-constructed version
}

// Specialization for uniform
template <>
inline auto createSymNormalRRD<std::uniform_real_distribution<double> >()
{
    return std::uniform_real_distribution<double>(-1., 1.); // the odd case with sigma!=1
}

// Specialization for gaussian
template <>
inline auto createSymNormalRRD<std::normal_distribution<double> >()
{
    return std::normal_distribution<double>(); // defaults to mean 0 std 1
}

// Specialization for student-t
template <>
inline auto createSymNormalRRD<std::student_t_distribution<double> >()
{
    return std::student_t_distribution<double>(); // default is symmetric around 0, but stddev undefined
}

// Specialization for cauchy
template <>
inline auto createSymNormalRRD<std::cauchy_distribution<double> >()
{
    return std::cauchy_distribution<double>(); // default is symmetric around 0, but mean&stddev undefined
}

// Specialization for symmetrized exponential distribution
template <>
inline auto createSymNormalRRD<SymmetrizedPRRD<std::exponential_distribution<double> > >()
{
    return SymmetrizedPRRD<std::exponential_distribution<double> >(); // use default (lambda=1.)
}

// Specialization for symmetrized gamma distribution
template <>
inline auto createSymNormalRRD<SymmetrizedPRRD<std::gamma_distribution<double> > >()
{
    return SymmetrizedPRRD<std::gamma_distribution<double> >(); // use default (alpha=1., beta=1.)
}

// Specialization for symmetrized weibull distribution
template <>
inline auto createSymNormalRRD<SymmetrizedPRRD<std::weibull_distribution<double> > >()
{
    return SymmetrizedPRRD<std::weibull_distribution<double> >(); // use default (alpha=1., beta=1.)
}

// Specialization for symmetrized lognormal distribution
template <>
inline auto createSymNormalRRD<SymmetrizedPRRD<std::lognormal_distribution<double> > >()
{
    return SymmetrizedPRRD<std::lognormal_distribution<double> >(); // use default (mu=0., s=1.)
}

// Specialization for symmetrized chi squared distribution
template <>
inline auto createSymNormalRRD<SymmetrizedPRRD<std::chi_squared_distribution<double> > >()
{
    return SymmetrizedPRRD<std::chi_squared_distribution<double> >(); // use default n=1
}

// Specialization for symmetrized fisher-f distribution
template <>
inline auto createSymNormalRRD<SymmetrizedPRRD<std::fisher_f_distribution<double> > >()
{
    return SymmetrizedPRRD<std::fisher_f_distribution<double> >(); // use default m=1, n=1
}
}; // namespace mci

#endif
