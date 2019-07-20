#ifndef MCI_TESTMCIFUNCTIONS_HPP
#define MCI_TESTMCIFUNCTIONS_HPP

#include "mci/ObservableFunctionInterface.hpp"
#include "mci/SamplingFunctionInterface.hpp"
#include "mci/WalkerState.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>


class TestWalk1s
{ // helps to generate random walk corresponding
    // to N particles in one-dimensional 1s orbital
protected:
    int _acc, _rej; // counters to calculate acceptance ratio

    bool _isStepAccepted(const double oldWFVal, double newWFVal) const
    {   // standard VMC acceptance criterion
        if (oldWFVal == 0) {
            return newWFVal != 0;
        }
        const double threshold = (newWFVal*newWFVal)/(oldWFVal*oldWFVal);
        if (threshold >= 1.) {
            return true;
        }
        return (rand()*(1.0/RAND_MAX) <= threshold);
    }

    double _calcWFVal(const double position[]) const
    {   // product of 1s orbitals in 1D
        double wfval = 0.;
        for (int i = 0; i < _ndim; ++i) {
            wfval += fabs(position[i]);
        }
        return exp(-wfval);
    }

    bool _generatePosition(const double oldPosition[], double newPosition[], int * nchanged, int changedIdx[]) const
    {
        double oldWFVal = _calcWFVal(oldPosition);

        if (nchanged != nullptr) { *nchanged = 0; }
        for (int i = 0; i < _ndim; ++i) {
            bool doMove = true;
            if (_changeProb < 1.) {
                doMove = (rand()*(1.0/RAND_MAX) <= _changeProb);
            }
            if (doMove) {
                newPosition[i] = oldPosition[i] + 2.*_stepSize*(rand()*(1.0/RAND_MAX) - 0.5);
                if (changedIdx != nullptr) { changedIdx[*nchanged] = i; }
                if (nchanged != nullptr) { *nchanged += 1; }
            }
            else {
                newPosition[i] = oldPosition[i];
            }
        }
        const double newWFVal = _calcWFVal(newPosition);

        return _isStepAccepted(oldWFVal, newWFVal);
    }

public:
    int _NMC;
    int _ndim;
    double _stepSize;
    double _changeProb; // 0..1, probability for single index to change on move

    TestWalk1s(int NMC, int ndim, double stepSize = 0.1, double changeProb = 1.):
            _acc(0), _rej(0), _NMC(NMC), _ndim(ndim), _stepSize(stepSize), _changeProb(changeProb) {}

    void generateWalk(double datax[] /*NMC*ndim length*/,
                      bool datacc[] = nullptr, /*if passed (NMC length), remember which steps were accepted*/
                      int nchanged[] = nullptr, /*if passed (NMC length), remember how many indices were changed at each step*/
                      int changedIdx[] = nullptr /*if passed (NMC*ndim length), remember which indices were changed*/
                     )
    {
        _acc = 0;
        _rej = 0;
        for (int j = 0; j < _ndim; ++j) { datax[j] = rand()*(1.0/RAND_MAX) - 0.5; } // set initial pos in [-0.5, 0.5)
        if (datacc != nullptr) { datacc[0] = true; }
        if (nchanged != nullptr) { nchanged[0] = _ndim; }
        if (changedIdx != nullptr) { std::iota(changedIdx, changedIdx + _ndim, 0); } // fill 0..ndim-1 on first step

        for (int i = 1; i < _NMC; ++i) {
            const bool accepted = _generatePosition(datax + (i - 1)*_ndim, datax + i*_ndim,
                                                    nchanged != nullptr ? nchanged + i : nullptr,
                                                    changedIdx != nullptr ? changedIdx + i*_ndim : nullptr);
            if (accepted) {
                ++_acc;
            }
            else {
                ++_rej;
                std::copy(datax + (i - 1)*_ndim, datax + i*_ndim, datax + i*_ndim); // copy old position
            }
            if (datacc != nullptr) { datacc[i] = accepted; }
        }
    }

    double getAcceptanceRate() { return static_cast<double>(_acc)/(static_cast<double>(_acc) + _rej); }
};


// --- SAMPLING FUNCTIONS

class ThreeDimGaussianPDF final: public mci::SamplingFunctionInterface
{
protected:
    mci::SamplingFunctionInterface * _clone() const final
    {
        return new ThreeDimGaussianPDF();
    }

public:
    ThreeDimGaussianPDF(): mci::SamplingFunctionInterface(3, 1) {}

    void protoFunction(const double in[], double protovalues[]) final
    {
        protovalues[0] = in[0]*in[0] + in[1]*in[1] + in[2]*in[2];
    }

    double samplingFunction(const double protov[]) const final
    {
        return exp(-protov[0]);
    }

    double acceptanceFunction(const double protoold[], const double protonew[]) const final
    {
        return exp(-protonew[0] + protoold[0]);
    }
};


class Gauss final: public mci::SamplingFunctionInterface
{
protected:
    mci::SamplingFunctionInterface * _clone() const final
    {
        return new Gauss(_ndim);
    }

public:
    explicit Gauss(const int ndim): mci::SamplingFunctionInterface(ndim, ndim) {}

    void protoFunction(const double in[], double out[]) final
    {
        for (int i = 0; i < _ndim; ++i) {
            out[i] = in[i]*in[i];
        }
    }

    double samplingFunction(const double protov[]) const final
    {
        return exp(-std::accumulate(protov, protov + _nproto, 0.));
    }

    double acceptanceFunction(const double protoold[], const double protonew[]) const final
    {
        double expf = std::accumulate(protoold, protoold + _nproto, 0.);
        expf -= std::accumulate(protonew, protonew + _nproto, 0.);
        return exp(expf);
    }

    double updatedAcceptance(const mci::WalkerState &wlk, const double pvold[], double pvnew[]) final
    {
        double expf = 0.;
        for (int i = 0; i < wlk.nchanged; ++i) {
            pvnew[wlk.changedIdx[i]] = wlk.xnew[wlk.changedIdx[i]]*wlk.xnew[wlk.changedIdx[i]];
            expf += pvnew[wlk.changedIdx[i]] - pvold[wlk.changedIdx[i]];
        }
        return exp(-expf);
    }
};

class Exp1DPDF final: public mci::SamplingFunctionInterface
{
protected:
    mci::SamplingFunctionInterface * _clone() const final
    {
        return new Exp1DPDF();
    }

public:
    Exp1DPDF(): mci::SamplingFunctionInterface(1, 1) {}

    void protoFunction(const double in[], double protovalues[]) final
    {
        protovalues[0] = fabs(in[0]);
    }

    double samplingFunction(const double protov[]) const final
    {
        return exp(-protov[0]);
    }

    double acceptanceFunction(const double protoold[], const double protonew[]) const final
    {
        return exp(-protonew[0] + protoold[0]);
    }
};

class ExpNDPDF final: public mci::SamplingFunctionInterface
{
protected:
    mci::SamplingFunctionInterface * _clone() const final
    {
        return new ExpNDPDF(_ndim);
    }

public:
    explicit ExpNDPDF(const int ndim): mci::SamplingFunctionInterface(ndim, ndim) {}

    void protoFunction(const double in[], double protovalues[]) final
    {
        for (int i = 0; i < _ndim; ++i) { protovalues[i] = fabs(in[i]); }
    }

    double samplingFunction(const double protov[]) const final
    {
        return exp(-std::accumulate(protov, protov + _nproto, 0.));
    }

    double acceptanceFunction(const double protoold[], const double protonew[]) const final
    {
        double expf = std::accumulate(protoold, protoold + _nproto, 0.);
        expf -= std::accumulate(protonew, protonew + _nproto, 0.);
        return exp(expf);
    }

    double updatedAcceptance(const mci::WalkerState &wlk, const double pvold[], double pvnew[]) final
    {
        double expf = 0.;
        for (int i = 0; i < wlk.nchanged; ++i) {
            pvnew[wlk.changedIdx[i]] = fabs(wlk.xnew[wlk.changedIdx[i]]);
            expf += pvnew[wlk.changedIdx[i]] - pvold[wlk.changedIdx[i]];
        }
        return exp(-expf);
    }
};

// --- OBSERVABLE FUNCTIONS

class XSquared final: public mci::ObservableFunctionInterface
{
protected:
    mci::ObservableFunctionInterface * _clone() const final
    {
        return new XSquared();
    }

public:
    XSquared(): mci::ObservableFunctionInterface(3, 1, false) {}

    void observableFunction(const double in[], double out[]) final
    {
        out[0] = (in[0]*in[0] + in[1]*in[1] + in[2]*in[2])/3.;
    }
};

// like XSquared, but multiplied by normalized Gaussian
class GaussXSquared final: public mci::ObservableFunctionInterface
{
protected:
    mci::ObservableFunctionInterface * _clone() const final
    {
        return new GaussXSquared();
    }
    const double _normf = 1./sqrt(M_PI*M_PI*M_PI)/3.; // /3. comes from the xsquared
public:
    GaussXSquared(): mci::ObservableFunctionInterface(3, 1, false) {}

    void observableFunction(const double in[], double out[]) final
    {
        const double x2 = in[0]*in[0] + in[1]*in[1] + in[2]*in[2];
        out[0] = exp(-x2)*x2*_normf;
    }
};


class XYZSquared final: public mci::ObservableFunctionInterface
{
protected:
    mci::ObservableFunctionInterface * _clone() const final
    {
        return new XYZSquared();
    }

public:
    XYZSquared(): mci::ObservableFunctionInterface(3, 3, false) {}

    void observableFunction(const double in[], double out[]) final
    {
        out[0] = in[0]*in[0];
        out[1] = in[1]*in[1];
        out[2] = in[2]*in[2];
    }
};


class X1D final: public mci::ObservableFunctionInterface
{
protected:
    mci::ObservableFunctionInterface * _clone() const final
    {
        return new X1D();
    }

public:
    X1D(): mci::ObservableFunctionInterface(1, 1, false) {}

    void observableFunction(const double in[], double out[]) final
    {
        out[0] = in[0];
    }
};


class XND final: public mci::ObservableFunctionInterface
{
protected:
    mci::ObservableFunctionInterface * _clone() const final
    {
        return new XND(_ndim);
    }

public:
    explicit XND(int nd): mci::ObservableFunctionInterface(nd, nd, false) {}

    void observableFunction(const double in[], double out[]) final
    {
        std::copy(in, in + _ndim, out);
    }
};

class UpdateableXND final: public mci::ObservableFunctionInterface
{
protected:
    mci::ObservableFunctionInterface * _clone() const final
    {
        return new UpdateableXND(_ndim);
    }

public:
    explicit UpdateableXND(int nd): mci::ObservableFunctionInterface(nd, nd, true) {}

    void observableFunction(const double in[], double out[]) final
    {
        std::copy(in, in + _ndim, out);
    }
    void updatedObservable(const double in[], const int/*nchanged*/, const bool flags[], double out[]) final
    {
        for (int i = 0; i < _ndim; ++i) { // this is likely slower in any case, but used for testing
            if (flags[i]) { out[i] = in[i]; }
        }
    }
};


class Constval final: public mci::ObservableFunctionInterface
{
protected:
    mci::ObservableFunctionInterface * _clone() const final
    {
        return new Constval(_ndim);
    }

public:
    explicit Constval(const int ndim): mci::ObservableFunctionInterface(ndim, 1, false) {}

    void observableFunction(const double /*in*/[], double out[]) final
    {
        out[0] = 1.3;
    }
};


class Polynom final: public mci::ObservableFunctionInterface
{
protected:
    mci::ObservableFunctionInterface * _clone() const final
    {
        return new Polynom(_ndim);
    }

public:
    explicit Polynom(const int ndim): mci::ObservableFunctionInterface(ndim, 1, false) {}

    void observableFunction(const double in[], double out[]) final
    {
        out[0] = 0.;
        for (int i = 0; i < _ndim; ++i) {
            out[0] += in[i];
        }
    }
};


class X2Sum final: public mci::ObservableFunctionInterface
{
protected:
    mci::ObservableFunctionInterface * _clone() const final
    {
        return new X2Sum(_ndim);
    }

public:
    explicit X2Sum(const int ndim): mci::ObservableFunctionInterface(ndim, 1, false) {}

    void observableFunction(const double in[], double out[]) final
    {
        out[0] = 0.;
        for (int i = 0; i < _ndim; ++i) {
            out[0] += in[i]*in[i];
        }
    }
};


class X2 final: public mci::ObservableFunctionInterface
{
protected:
    mci::ObservableFunctionInterface * _clone() const final
    {
        return new X2(_ndim);
    }

public:
    explicit X2(const int ndim): mci::ObservableFunctionInterface(ndim, ndim, true) {}

    void observableFunction(const double in[], double out[]) final
    {
        for (int i = 0; i < _ndim; ++i) {
            out[i] = in[i]*in[i];
        }
    }

    void updatedObservable(const double in[], const int/*nchanged*/, const bool flags[], double out[]) final
    {
        for (int i = 0; i < _ndim; ++i) {
            if (flags[i]) { // this may actually be faster for small nchanged and large _ndim
                out[i] = in[i]*in[i];
            }
        }
    }
};

#endif
