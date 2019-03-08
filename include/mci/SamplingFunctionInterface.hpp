#ifndef MCI_SAMPLINGFUNCTIONINTERFACE_HPP
#define MCI_SAMPLINGFUNCTIONINTERFACE_HPP

namespace mci
{
    class SamplingFunctionInterface
    {
    protected:
        const int _ndim; //dimension of the input array (walker position)
        int _nproto; //number of proto sampling functions given as output
        double * _protonew; //array containing the new proto sampling functions
        double * _protoold; //array containing the old proto sampling functions

        // Setters
        void setNProto(int nproto); // you may freely choose the amount of values you need

    public:
        SamplingFunctionInterface(int ndim, int nproto);

        virtual ~SamplingFunctionInterface();

        // Getters
        int getNDim() const { return _ndim;}
        int getNProto() const { return _nproto;}

        // Main operational methods
        void newToOld(); // swap old and new protovalues
        void computeNewSamplingFunction(const double in[]) { samplingFunction(in, _protonew); }
        double getAcceptance() const { return getAcceptance(_protoold, _protonew); }


        // --- METHODS THAT MUST BE IMPLEMENTED
        // Function that MCI uses for the proto-sampling function. Computes _protonew
        virtual void samplingFunction(const double in[], double protovalues[]) = 0;
        //                                      ^walker position  ^resulting proto-values

        // Acceptance function, that uses the new and old values of the proto sampling function
        virtual double getAcceptance(const double protoold[], const double protonew[]) const = 0;
    };
}  // namespace mci

#endif
