#include "mci/MCIObservableFunctionInterface.hpp"
#include "mci/MCISamplingFunctionInterface.hpp"

#include <cmath>

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
