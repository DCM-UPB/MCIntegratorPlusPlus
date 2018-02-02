#ifndef MCI_SAMPLING_FUNCTION_INTERFACE
#define MCI_SAMPLING_FUNCTION_INTERFACE


class MCISamplingFunctionInterface
{
   protected:
      int _ndim; //dimension of the input array (walker position)
      int _nproto; //number of proto sampling functions given as output
      double * _protonew; //array containing the new proto sampling functions
      double * _protoold; //array containing the old proto sampling functions

   public:
      MCISamplingFunctionInterface(const int &ndim, const int &nproto)
      {
         _ndim=ndim; _nproto=nproto;
         _protonew = new double[nproto]; _protoold = new double[nproto];
         for (int i=0; i<nproto; ++i){ _protonew[i]=0.; }
         for (int i=0; i<nproto; ++i){ _protoold[i]=0.; }
      }
      virtual ~MCISamplingFunctionInterface()
      {
         delete[] _protonew; delete[] _protoold;
      }

      // Getters
      int getNDim(){ return _ndim;}
      int getNProto(){ return _nproto;}
      double getProtoNew(const int &i){ return _protonew[i]; }
      double getProtoOld(const int &i){ return _protoold[i]; }

      // Utilities
      void newToOld()
      {
         for(int i=0; i<_nproto; ++i){ double * foo; foo=_protonew; _protonew=_protoold; _protoold=foo; }
      }

      void computeNewSamplingFunction(const double * in)
      {
         samplingFunction(in, _protonew);
      }

      // --- METHODS THAT MUST BE IMPLEMENTED
      // Function that MCI uses for the proto-sampling function. Computes _protonew
      virtual void samplingFunction(const double *in, double * protovalues) = 0;
      //                                      ^walker position  ^resulting proto-values

      // Acceptance function, that uses the new and old values of the proto sampling function
      virtual double getAcceptance() = 0;
};


#endif
