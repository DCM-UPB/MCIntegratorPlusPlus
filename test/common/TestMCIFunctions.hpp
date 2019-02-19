#ifndef TEST_MCI_FUNCTIONS
#define TEST_MCI_FUNCTIONS

class ThreeDimGaussianPDF: public MCISamplingFunctionInterface{
public:
    ThreeDimGaussianPDF(): MCISamplingFunctionInterface(3, 1){}
    ~ThreeDimGaussianPDF(){}

    void samplingFunction(const double *in, double * protovalues){
        protovalues[0] = (in[0]*in[0]) + (in[1]*in[1]) + (in[2]*in[2]);
    }


    double getAcceptance(const double * protoold, const double * protonew){
        return exp(-protonew[0]+protoold[0]);
    }
};

class XSquared: public MCIObservableFunctionInterface{
public:
    XSquared(): MCIObservableFunctionInterface(3, 1){}
    ~XSquared(){}

    void observableFunction(const double * in, double * out){
        out[0] = in[0] * in[0];
    }
};


class XYZSquared: public MCIObservableFunctionInterface{
public:
    XYZSquared(): MCIObservableFunctionInterface(3, 3){}
    ~XYZSquared(){}

    void observableFunction(const double * in, double * out){
        out[0] = in[0] * in[0];
        out[1] = in[1] * in[1];
        out[2] = in[2] * in[2];
    }
};

#endif
