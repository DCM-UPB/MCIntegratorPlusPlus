#include "Estimators.hpp"

#include <iostream>
#include <math.h>
#include <cmath>

namespace mci
{
   void UncorrelatedEstimator(const long &n, const double * x, double * average, double * error)
   {
      using namespace std;
      
      const double SMALLEST_ERROR=1.e-300;

      if ( n < 2)
      {
         cout << "MCI error UncorrelatedEstimator() : n must be larger than 1";
         exit(EXIT_FAILURE);
      }

      *average=0.;
      *error=0.;
      for (long i=0; i<n; ++i)
      {
         *average+=*(x+i);
         *error+=(*(x+i))*(*(x+i));
      }

      double norm = 1./((double) n);
      *average*=norm;
      *error*=norm;
      if ( ((*(error))-(*(average))*(*(average))) < SMALLEST_ERROR )
      {
         *(error)=0.;
      } else 
      {
         *(error)=sqrt( ((*(error))-(*(average))*(*(average)))/((double) (n-1)) );
      }
   }

   
   void BlockEstimator(const long &n, const double * x, const int &nblocks, double * average, double * error)
   {
      using namespace std;

      if ( n < 2)
      {
         cout << "MCI error BlockEstimator() : n must be larger than 1";
         exit(EXIT_FAILURE);
      }

      long nperblock=(n)/(nblocks);
      double norm = 1./((double) (nperblock));

      double *av = new double[nblocks];
      for (int i1=0; i1<(nblocks); ++i1)
      {
           *(av+i1)=0;
           for (long i2=(i1)*(nperblock); i2<(i1+1)*(nperblock); ++i2)
           {
               *(av+i1)+=*(x+i2);
           }
           *(av+i1)*=norm;
      }

      long n_ue=nblocks;
      UncorrelatedEstimator(n_ue, av, average, error);

      delete [] av;
   }


   void CorrelatedEstimator(const long &n, const double * x, double * average, double * error)
   {
      const int MIN_BLOCKS=6, MAX_BLOCKS=50;
      const int MAX_PLATEAU_AVERAGE=4;

      using namespace std;

      if ( n < 2)
      {
         cout << "MCI error CorrelatedEstimator() : n must be larger than 1";
         exit(EXIT_FAILURE);
      }

      double *av = new double[MAX_BLOCKS-MIN_BLOCKS+1];
      double *err = new double[MAX_BLOCKS-MIN_BLOCKS+1];
      int nblocks;
      for (int i1=0; i1<MAX_BLOCKS-MIN_BLOCKS+1; ++i1)
      {
         nblocks=i1+MIN_BLOCKS;
         BlockEstimator(n, x, nblocks, av+i1, err+i1);
         //cout << "Nblocks = " << nblocks << "   average = " << *(av+i1) << "   error = " << *(err+i1) << endl;
      }

      double delta = 0.;
      double *accdelta = new double[MAX_BLOCKS-MIN_BLOCKS+1-2*MAX_PLATEAU_AVERAGE];
      for (int i2=0; i2<MAX_BLOCKS-MIN_BLOCKS+1-2*MAX_PLATEAU_AVERAGE; ++i2)
      {
         *(accdelta+i2) = 0.;
      }
      for (int i1=1; i1<=MAX_PLATEAU_AVERAGE; ++i1)
      {
         for (int i2=MAX_PLATEAU_AVERAGE; i2<MAX_BLOCKS-MIN_BLOCKS+1-MAX_PLATEAU_AVERAGE; ++i2)
         {
            switch(i1)
            {
               case 1:
                  delta = ( -0.5*(*(err+i2-1)) + 0.5*(*(err+i2+1)) );
                  break;
               case 2:
                  delta = ( (1./12.)*(*(err+i2-2)) -(2./3.)*(*(err+i2-1)) 
                            +(2./3.)*(*(err+i2+1)) -(1./12.)*(*(err+i2+2)) );
                  break;
               case 3:
                  delta = ( -(1./60.)*(*(err+i2-3)) +(3./20.)*(*(err+i2-2)) -0.75*(*(err+i2-1)) 
                            +0.75*(*(err+i2+1)) -(3./20.)*(*(err+i2+2)) +(1./60.)*(*(err+i2+3))  );
                  break;
               case 4:
                  delta = ( (1./280.)*(*(err+i2-4)) -(4./105.)*(*(err+i2-3)) +0.2*(*(err+i2-2)) -0.8*(*(err+i2-1))
                            +0.8*(*(err+i2+1)) -0.2*(*(err+i2+2)) +(4./105.)*(*(err+i2+3)) -(1./280.)*(*(err+i2+4)) );
                  break;
            }
            *(accdelta+i2-MAX_PLATEAU_AVERAGE)+=delta;
         }
      }

      /*for (int i2=MAX_PLATEAU_AVERAGE; i2<MAX_BLOCKS-MIN_BLOCKS+1-MAX_PLATEAU_AVERAGE; ++i2)
      {
         cout << "i = " << i2+MIN_BLOCKS << "   accdelta = " << *(accdelta+i2-MAX_PLATEAU_AVERAGE) << endl;
      }*/

      int i_min = 0;
      for (int i2=1; i2<MAX_BLOCKS-MIN_BLOCKS+1-2*MAX_PLATEAU_AVERAGE; ++i2)
      {
         if ( abs(*(accdelta+i2)) < abs(*(accdelta+i_min)) ) i_min=i2;
      }

      i_min+=MAX_PLATEAU_AVERAGE;
      //cout << "The plateau has been detected at nblocks = " << i_min+MIN_BLOCKS << endl;
      *average = 0.2*( *(av+i_min-2) + *(av+i_min-1) + *(av+i_min) + *(av+i_min+1) + *(av+i_min+2) );
      *error = 0.2*( *(err+i_min-2) + *(err+i_min-1) + *(err+i_min) + *(err+i_min+1) + *(err+i_min+2) );

      delete [] av;
      delete [] err;
      delete [] accdelta;
   }


   void MultiDimUncorrelatedEstimator(const long &n, const int &ndim, const double * const * x, double * average, double * error)
   {
      using namespace std;
      
      const double SMALLEST_ERROR=1.e-300;

      if ( n < 2)
      {
         cout << "MCI error UncorrelatedEstimator() : n must be larger than 1";
         exit(EXIT_FAILURE);
      }

      for (int j=0; j<ndim; ++j)
      {
         *(average+j)=0.;
         *(error+j)=0.;
      }

      for (long i=0; i<n; ++i)
      {
         for (int j=0; j<ndim; ++j)
         {
            *(average+j) += *(*(x+i)+j) ;
            *(error+j) +=( (*(*(x+i)+j)) * (*(*(x+i)+j)) );
         }
      }

      double norm = 1./((double) n);
      double norm2 = (double) (n-1);
      for (int j=0; j<ndim; ++j)
      {
         *(average+j)*=norm;
         *(error+j)*=norm;
         if ( ((*(error+j))-(*(average+j))*(*(average+j))) < SMALLEST_ERROR )
         {
            *(error+j)=0.;
         } else 
         {
            *(error+j)=sqrt( ((*(error+j))-(*(average+j))*(*(average+j)))/(norm2) );
         }
      }
   }


   void MultiDimBlockEstimator(const long &n, const int &ndim, const double * const * x, const int &nblocks, double * average, double * error)
   {
      using namespace std;

      if ( n < 2)
      {
         cout << "MCI error BlockEstimator() : n must be larger than 1";
         exit(EXIT_FAILURE);
      }

      long nperblock=(n)/(nblocks);
      double norm = 1./((double) nperblock);

      double ** av = new double*[nblocks];
      for (int j=0; j<nblocks; ++j)
      { 
         *(av+j) = new double[ndim]; 
      }

      for (int i1=0; i1<(nblocks); ++i1)
      {
         for (int j=0; j<ndim; ++j)
         {
            *(*(av+i1)+j)=0;
         }

         for (long i2=(i1)*(nperblock); i2<(i1+1)*(nperblock); ++i2)
         {
            for (int j=0; j<ndim; ++j)
            {
               (*(*(av+i1)+j)) += (*(*(x+i2)+j));
            }
         }

         for (int j=0; j<ndim; ++j)
         {
            *(*(av+i1)+j) *= norm;
         }
      }

      long n_ue=(nblocks);
      MultiDimUncorrelatedEstimator(n_ue, ndim, av, average, error);

      for (int j=0; j<nblocks; ++j){ delete [] *(av+j); }
      delete [] av;

   }


   void MultiDimCorrelatedEstimator(const long &n, const int &ndim, const double * const * x, double * average, double * error)
   {
      const int MIN_BLOCKS=6, MAX_BLOCKS=50;
      const int MAX_PLATEAU_AVERAGE=4;

      using namespace std;

      if ( n < 2)
      {
         cout << "MCI error CorrelatedEstimator() : n must be larger than 1";
         exit(EXIT_FAILURE);
      }

      double ** av = new double*[MAX_BLOCKS-MIN_BLOCKS+1];
      for (int j=0; j<MAX_BLOCKS-MIN_BLOCKS+1; ++j){ *(av+j) = new double[ndim]; }
      
      double ** err = new double*[MAX_BLOCKS-MIN_BLOCKS+1];
      for (int j=0; j<MAX_BLOCKS-MIN_BLOCKS+1; ++j){ *(err+j) = new double[ndim]; }

      int nblocks;
      for (int i1=0; i1<MAX_BLOCKS-MIN_BLOCKS+1; ++i1)
      {
         nblocks=i1+MIN_BLOCKS;
         MultiDimBlockEstimator(n, ndim, x, nblocks, *(av+i1), *(err+i1));
      }

      double *delta=new double[ndim];
      for (int j=0; j<ndim; ++j) {*(delta+j)=0.;}

      double **accdelta = new double*[MAX_BLOCKS-MIN_BLOCKS+1-2*MAX_PLATEAU_AVERAGE];
      for (int j=0; j<MAX_BLOCKS-MIN_BLOCKS+1-2*MAX_PLATEAU_AVERAGE; ++j){ *(accdelta+j)=new double[ndim]; }
      for (int i2=0; i2<MAX_BLOCKS-MIN_BLOCKS+1-2*MAX_PLATEAU_AVERAGE; ++i2)
      {
         for (int j=0; j<ndim; ++j){*(*(accdelta+i2)+j)=0.;}
      }

      for (int i1=1; i1<=MAX_PLATEAU_AVERAGE; ++i1)
      {
         for (int i2=MAX_PLATEAU_AVERAGE; i2<MAX_BLOCKS-MIN_BLOCKS+1-MAX_PLATEAU_AVERAGE; ++i2)
         {
            for (int j=0; j<ndim; ++j)
            {
               switch(i1)
               {
                  case 1:
                     *(delta+j) = ( -0.5*(*(*(err+i2-1)+j)) + 0.5*(*(*(err+i2+1)+j)) );
                     break;
                  case 2:
                     *(delta+j) = ( (1./12.)*(*(*(err+i2-2)+j)) -(2./3.)*(*(*(err+i2-1)+j)) 
                               +(2./3.)*(*(*(err+i2+1)+j)) -(1./12.)*(*(*(err+i2+2)+j)) );
                     break;
                  case 3:
                     *(delta+j) = ( -(1./60.)*(*(*(err+i2-3)+j)) +(3./20.)*(*(*(err+i2-2)+j)) -0.75*(*(*(err+i2-1)+j)) 
                               +0.75*(*(*(err+i2+1)+j)) -(3./20.)*(*(*(err+i2+2)+j)) +(1./60.)*(*(*(err+i2+3)+j))  );
                     break;
                  case 4:
                     *(delta+j) = ( (1./280.)*(*(*(err+i2-4)+j)) -(4./105.)*(*(*(err+i2-3)+j)) +0.2*(*(*(err+i2-2)+j)) -0.8*(*(*(err+i2-1)+j))
                               +0.8*(*(*(err+i2+1)+j)) -0.2*(*(*(err+i2+2)+j)) +(4./105.)*(*(*(err+i2+3)+j)) -(1./280.)*(*(*(err+i2+4)+j)) );
                     break;
               }
               *(*(accdelta+i2-MAX_PLATEAU_AVERAGE)+j) += *(delta+j);
            }
         }
      }

      int * i_min = new int[ndim];
      for (int j=0; j<ndim; ++j){ *(i_min+j)=0; }
      for (int i2=1; i2<MAX_BLOCKS-MIN_BLOCKS+1-2*MAX_PLATEAU_AVERAGE; ++i2)
      {
         for (int j=0; j<ndim; ++j)
         {
            if ( abs(*(*(accdelta+i2)+j)) < abs( *(*(accdelta+(*(i_min+j)))+j) ) ) *(i_min+j)=i2;
         }
      }
      for (int j=0; j<ndim; ++j){ *(i_min+j) += MAX_PLATEAU_AVERAGE; }

      for (int j=0; j<ndim; ++j)
      {
         *(average+j) = 0.2*( *(*(av+(*(i_min+j))-2)+j) + *(*(av+(*(i_min+j))-1)+j) + *(*(av+(*(i_min+j)))+j) + *(*(av+(*(i_min+j))+1)+j) + *(*(av+(*(i_min+j))+2)+j) );
         *(error+j) = 0.2*( *(*(err+(*(i_min+j))-2)+j) + *(*(err+(*(i_min+j))-1)+j) + *(*(err+(*(i_min+j)))+j) + *(*(err+(*(i_min+j))+1)+j) + *(*(err+(*(i_min+j))+2)+j) );
      }

      delete[] i_min;

      for (int j=0; j<MAX_BLOCKS-MIN_BLOCKS+1-2*MAX_PLATEAU_AVERAGE; ++j){ delete [] *(accdelta+j); };
      delete [] accdelta;

      delete [] delta;

      for (int j=0; j<MAX_BLOCKS-MIN_BLOCKS+1 ; ++j){ delete [] *(av+j); }
      delete [] av;

      for (int j=0; j<MAX_BLOCKS-MIN_BLOCKS+1 ; ++j){ delete[] *(err+j); }
      delete[] err;
   }

}
