#ifndef MCI_EXAMPLEFUNCTIONS_HPP
#define MCI_EXAMPLEFUNCTIONS_HPP

#include "mci/SamplingFunctionInterface.hpp"
#include "mci/ObservableFunctionInterface.hpp"

#include <cmath>


// Observable functions
class Parabola final: public mci::ObservableFunctionInterface
{
protected:
    // Observables need to be copyable, so we need to provide this
    // protected method. In return, there will be a public method
    // clone() returning a std::unique_ptr<ObservableFunctionInterface> .
    mci::ObservableFunctionInterface * _clone() const final
    {
        return new Parabola();
    }
public:
    Parabola(): mci::ObservableFunctionInterface(1 /*1D input*/, 1 /*1D output*/, false/*no partial updating*/) {}

    // here we calculate the observable function
    void observableFunction(const double in[], double out[]) final
    {
        out[0] = 4.*in[0] - in[0]*in[0];
    }
};


class NormalizedParabola final: public mci::ObservableFunctionInterface
{
protected:
    // same as above
    mci::ObservableFunctionInterface * _clone() const final
    {
        return new NormalizedParabola();
    }
public:
    explicit NormalizedParabola(): mci::ObservableFunctionInterface(1, 1, false/*no partial updating*/) {}

    void observableFunction(const double in[], double out[]) final
    {
        out[0] = (4. - in[0])*5.;
        if (std::signbit(in[0])) { out[0] = -out[0]; }
    }
};


// Sampling function
class NormalizedLine final: public mci::SamplingFunctionInterface
{
protected:
    mci::SamplingFunctionInterface * _clone() const final
    {
        return new NormalizedLine();
    }
public:
    explicit NormalizedLine(): mci::SamplingFunctionInterface(1, 1) {}

    void protoFunction(const double in[], double protovalue[]) final
    {
        protovalue[0] = 0.2*fabs(in[0]);
    }

    double samplingFunction(const double protovalue[]) const final
    {
        return protovalue[0];
    }

    double acceptanceFunction(const double protoold[], const double protonew[]) const final
    { // don't forget the const!
        return protonew[0]/protoold[0];
    }
};

#endif
