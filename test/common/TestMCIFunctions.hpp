#ifndef MCI_TESTMCIFUNCTIONS_HPP
#define MCI_TESTMCIFUNCTIONS_HPP

#include "mci/ObservableFunctionInterface.hpp"
#include "mci/SamplingFunctionInterface.hpp"
#include "mci/UpdateableObservableInterface.hpp"
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

    bool _isStepAccepted(const double oldWFVal, double newWFVal)
    {   // standard VMC acceptance criterion
        if (oldWFVal == 0) {
            return true;
        }
        if (newWFVal == 0) {
            return false;
        }
        const double threshold = (newWFVal*newWFVal)/(oldWFVal*oldWFVal);
        if (threshold >= 1.) {
            return true;
        }
        return ( rand()*(1.0 / RAND_MAX) <= threshold );
    }

    double _calcWFVal(const double position[])
    {   // product of 1s orbitals in 1D
        double wfval = 0.;
        for (int i=0; i<_ndim; ++i) {
            wfval += fabs(position[i]);
        }
        return exp(-wfval);
    }

    bool _generatePosition(const double oldPosition[], double newPosition[])
    {
        double oldWFVal = _calcWFVal(oldPosition);

        for (int i=0; i<_ndim; ++i) {
            newPosition[i] = oldPosition[i] + 2.*_stepSize*(rand()*(1.0 / RAND_MAX) - 0.5);
        }
        const double newWFVal = _calcWFVal(newPosition);

        return _isStepAccepted(oldWFVal, newWFVal);
    }

public:
    int _NMC;
    int _ndim;
    double _stepSize;

    TestWalk1s(int NMC, int ndim, double stepSize = 0.1):
        _acc(0), _rej(0), _NMC(NMC), _ndim(ndim), _stepSize(stepSize) {}

    void generateWalk(double datax[] /*NMC*ndim shape*/, bool datacc[] = nullptr /*if passed (NMC length), remember which steps where new ones*/)
    {
        _acc = 0;
        _rej = 0;
        for (int j=0; j<_ndim; ++j) { datax[j] = rand()*(1.0 / RAND_MAX) - 0.5; } // set initial pos in [-0.5, 0.5)
        if (datacc != nullptr) { datacc[0] = true; }

        for (int i=1; i<_NMC; ++i) {
            const bool accepted = _generatePosition(datax+(i-1)*_ndim, datax+i*_ndim);
            if (accepted) {
                ++_acc;
            } else {
                ++_rej;
                std::copy(datax+(i-1)*_ndim, datax+i*_ndim, datax+i*_ndim); // copy old position
            }
            if (datacc != nullptr) { datacc[i] = accepted; }
        }
    }

    double getAcceptanceRate() { return static_cast<double>(_acc)/( static_cast<double>(_acc) + _rej ); }
};


// --- SAMPLING FUNCTIONS

class ThreeDimGaussianPDF: public mci::SamplingFunctionInterface
{
protected:
    mci::SamplingFunctionInterface * _clone() const final {
        return new ThreeDimGaussianPDF();
    }

public:
    ThreeDimGaussianPDF(): mci::SamplingFunctionInterface(3, 1){}
    ~ThreeDimGaussianPDF() final= default;

    void protoFunction(const double in[], double protovalues[]) final{
        protovalues[0] = (in[0]*in[0]) + (in[1]*in[1]) + (in[2]*in[2]);
    }

    double samplingFunction(const double protov[]) const final
    {
        return exp(-protov[0]);
    }

    double acceptanceFunction(const double protoold[], const double protonew[]) const final{
        return exp(-protonew[0]+protoold[0]);
    }
};


class Gauss: public mci::SamplingFunctionInterface
{
protected:
    mci::SamplingFunctionInterface * _clone() const final {
        return new Gauss(_ndim);
    }

public:
    explicit Gauss(const int ndim): mci::SamplingFunctionInterface(ndim,ndim)
    {}

    void protoFunction(const double in[], double out[]) final
    {
        for (int i=0; i<_ndim; ++i) {
            out[i] = in[i]*in[i];
        }
    }

    double samplingFunction(const double protov[]) const final
    {
        return exp(-std::accumulate(protov, protov+_nproto, 0.));
    }

    double acceptanceFunction(const double protoold[], const double protonew[]) const final
    {
        double expf = std::accumulate(protoold, protoold+_nproto, 0.);
        expf -= std::accumulate(protonew, protonew+_nproto, 0.);
        return exp(expf);
    }

    double updatedAcceptance(const mci::WalkerState &wlk, const double pvold[], double pvnew[]) final
    {
        double expf = 0.;
        for (int i=0; i<wlk.nchanged; ++i) {
            pvnew[wlk.changedIdx[i]] = wlk.xnew[wlk.changedIdx[i]] * wlk.xnew[wlk.changedIdx[i]];
            expf += pvnew[wlk.changedIdx[i]] - pvold[wlk.changedIdx[i]];
        }
        return exp(-expf);
    }
};

class Exp1DPDF: public mci::SamplingFunctionInterface
{
protected:
    mci::SamplingFunctionInterface * _clone() const final {
        return new Exp1DPDF();
    }

public:
    Exp1DPDF(): mci::SamplingFunctionInterface(1, 1){}
    ~Exp1DPDF() final= default;

    void protoFunction(const double in[], double protovalues[]) final{
        protovalues[0] = fabs(in[0]);
    }

    double samplingFunction(const double protov[]) const final
    {
        return exp(-protov[0]);
    }

    double acceptanceFunction(const double protoold[], const double protonew[]) const final{
        return exp(-protonew[0]+protoold[0]);
    }
};

class ExpNDPDF: public mci::SamplingFunctionInterface
{
protected:
    mci::SamplingFunctionInterface * _clone() const final {
        return new ExpNDPDF(_ndim);
    }

public:
    explicit ExpNDPDF(const int ndim): mci::SamplingFunctionInterface(ndim,ndim)
    {}

    void protoFunction(const double in[], double protovalues[]) final
    {
        for (int i=0; i<_ndim; ++i) { protovalues[i] = fabs(in[i]); }
    }

    double samplingFunction(const double protov[]) const final
    {
        return exp(-std::accumulate(protov, protov+_nproto, 0.));
    }

    double acceptanceFunction(const double protoold[], const double protonew[]) const final
    {
        double expf = std::accumulate(protoold, protoold+_nproto, 0.);
        expf -= std::accumulate(protonew, protonew+_nproto, 0.);
        return exp(expf);
    }

    double updatedAcceptance(const mci::WalkerState &wlk, const double pvold[], double pvnew[]) final
    {
        double expf = 0.;
        for (int i=0; i<wlk.nchanged; ++i) {
            pvnew[wlk.changedIdx[i]] = wlk.xnew[wlk.changedIdx[i]];
            expf += pvnew[wlk.changedIdx[i]] - pvold[wlk.changedIdx[i]];
        }
        return exp(-expf);
    }
};

// --- OBSERVABLE FUNCTIONS

class XSquared: public mci::ObservableFunctionInterface
{
protected:
    mci::ObservableFunctionInterface * _clone() const final {
        return new XSquared();
    }

public:
    XSquared(): mci::ObservableFunctionInterface(3, 1){}
    ~XSquared() final= default;

    void observableFunction(const double in[], double out[]) final{
        out[0] = in[0] * in[0];
    }
};


class XYZSquared: public mci::ObservableFunctionInterface
{
protected:
    mci::ObservableFunctionInterface * _clone() const final {
        return new XYZSquared();
    }

public:
    XYZSquared(): mci::ObservableFunctionInterface(3, 3){}
    ~XYZSquared() final= default;

    void observableFunction(const double in[], double out[]) final {
        out[0] = in[0] * in[0];
        out[1] = in[1] * in[1];
        out[2] = in[2] * in[2];
    }
};


class X1D: public mci::ObservableFunctionInterface
{
protected:
    mci::ObservableFunctionInterface * _clone() const final {
        return new X1D();
    }

public:
    X1D(): mci::ObservableFunctionInterface(1, 1){}
    ~X1D() final= default;

    void observableFunction(const double in[], double out[]) final{
        out[0] = in[0];
    }
};

class XND: public mci::UpdateableObservableInterface
{
protected:
    mci::ObservableFunctionInterface * _clone() const final {
        return new XND(_ndim);
    }

public:
    explicit XND(int nd): mci::UpdateableObservableInterface(nd, nd){}
    ~XND() final= default;

    void observableFunction(const double in[], double out[]) final {
        std::copy(in, in+_ndim, out);
    }
    void updatedObservable(const double in[], const int/*nchanged*/, const bool flags[], double out[]) final {
        for (int i=0; i<_ndim; ++i) { // this is likely slower in any case, but used for testing
            if (flags[i]) { out[i] = in[i]; }
        }
    }
};


class Constval: public mci::ObservableFunctionInterface
{
protected:
    mci::ObservableFunctionInterface * _clone() const final {
        return new Constval(_ndim);
    }

public:
    explicit Constval(const int ndim): mci::ObservableFunctionInterface(ndim, 1) {}

    void observableFunction(const double /*in*/[], double out[]) final
    {
        out[0] = 1.3;
    }
};


class Polynom: public mci::ObservableFunctionInterface
{
protected:
    mci::ObservableFunctionInterface * _clone() const final {
        return new Polynom(_ndim);
    }

public:
    explicit Polynom(const int ndim): mci::ObservableFunctionInterface(ndim, 1) {}

    void observableFunction(const double in[], double out[]) final
    {
        out[0]=0.;
        for (int i=0; i<_ndim; ++i) {
            out[0] += in[i];
        }
    }
};


class X2Sum: public mci::ObservableFunctionInterface
{
protected:
    mci::ObservableFunctionInterface * _clone() const final {
        return new X2Sum(_ndim);
    }

public:
    explicit X2Sum(const int ndim): mci::ObservableFunctionInterface(ndim,1) {}

    void observableFunction(const double in[], double out[]) final
    {
        out[0] = 0.;
        for (int i=0; i<_ndim; ++i) {
            out[0] += in[i]*in[i];
        }
    }
};


class X2: public mci::UpdateableObservableInterface
{
protected:
    mci::ObservableFunctionInterface * _clone() const final {
        return new X2Sum(_ndim);
    }

public:
    explicit X2(const int ndim): mci::UpdateableObservableInterface(ndim,ndim) {}

    void observableFunction(const double in[], double out[]) final
    {
        for (int i=0; i<_ndim; ++i) {
            out[i] = in[i]*in[i];
        }
    }

    void updatedObservable(const double in[], const int/*nchanged*/, const bool flags[], double out[]) final
    {
        for (int i=0; i<_ndim; ++i) {
            if (flags[i]) { // this may actually be faster for small nchanged and large _ndim
                out[i] = in[i]*in[i];
            }
        }
    }
};

#endif
