#include "mci/MCIObservableFunctionInterface.hpp"
#include "mci/MCISamplingFunctionInterface.hpp"

#include <cmath>
#include <algorithm>
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

    double _calcWFVal(const double * position)
    {   // product of 1s orbitals in 1D
        double wfval = 0.;
        for (int i=0; i<_ndim; ++i) {
            wfval += fabs(position[i]);
        }
        return exp(-wfval);
    }

    bool _generatePosition(const double * oldPosition, double * newPosition)
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

    void generateWalk(double * datax /* NMC*ndim shape */, bool * datacc = nullptr /* if passed (NMC length), remember which steps where new ones */)
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


class ThreeDimGaussianPDF: public MCISamplingFunctionInterface{
public:
    ThreeDimGaussianPDF(): MCISamplingFunctionInterface(3, 1){}
    ~ThreeDimGaussianPDF() override= default;

    void samplingFunction(const double *in, double * protovalues) override{
        protovalues[0] = (in[0]*in[0]) + (in[1]*in[1]) + (in[2]*in[2]);
    }


    double getAcceptance(const double * protoold, const double * protonew) override{
        return exp(-protonew[0]+protoold[0]);
    }
};

class XSquared: public MCIObservableFunctionInterface{
public:
    XSquared(): MCIObservableFunctionInterface(3, 1){}
    ~XSquared() override= default;

    void observableFunction(const double * in, double * out) override{
        out[0] = in[0] * in[0];
    }
};


class XYZSquared: public MCIObservableFunctionInterface{
public:
    XYZSquared(): MCIObservableFunctionInterface(3, 3){}
    ~XYZSquared() override= default;

    void observableFunction(const double * in, double * out) override{
        out[0] = in[0] * in[0];
        out[1] = in[1] * in[1];
        out[2] = in[2] * in[2];
    }
};

class XND: public MCIObservableFunctionInterface{
public:
    XND(int nd): MCIObservableFunctionInterface(nd, nd){}
    ~XND() override= default;

    void observableFunction(const double * in, double * out) override{
        std::copy(in, in+_ndim, out);
    }
};
