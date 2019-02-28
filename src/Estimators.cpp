#include "mci/Estimators.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

double calcErrDelta(const int mode, const double * err /* length 9 */)
{   // for Francesco's plateau finding algorithm
    const int im = 4; // index of middle element
    switch(mode)
        {
        case 1:
            return ( -0.5*err[im-1] + 0.5*err[im+1] );
        case 2:
            return ( (1./12.)*err[im-2] -(2./3.)*err[im-1]
                     +(2./3.)*err[im+1] -(1./12.)*err[im+2] );
        case 3:
            return ( -(1./60.)*err[im-3] +(3./20.)*err[im-2] -0.75*err[im-1]
                     +0.75*err[im+1] -(3./20.)*err[im+2] +(1./60.)*err[im+3]  );
        case 4:
            return ( (1./280.)*err[im-4] -(4./105.)*err[im-3] +0.2*err[im-2] -0.8*err[im-1]
                     +0.8*err[im+1] -0.2*err[im+2] +(4./105.)*err[im+3] -(1./280.)*err[im+4] );
        }
    return 0.;
}


namespace mci
{
    void UncorrelatedEstimator(const int &n, const double * x, double & average, double & error)
    {
        using namespace std;

        const double SMALLEST_ERROR=1.e-300;

        if ( n < 2) {
            cout << "MCI error UncorrelatedEstimator() : n must be larger than 1";
            exit(EXIT_FAILURE);
        }

        average = accumulate(x, x+n, 0.);
        error = inner_product(x, x+n, x, 0.);

        const double norm = 1./n;
        average*=norm;
        error=error*norm - average*average;
        if ( error > SMALLEST_ERROR ){
            error=sqrt( error/(n-1.) );
        } else {
            error=0.;
        }
    }


    void BlockEstimator(const int &n, const double * x, const int &nblocks, double & average, double & error)
    {
        using namespace std;

        if ( n < nblocks) {
            cout << "MCI error BlockEstimator() : n must be >= nblocks";
            exit(EXIT_FAILURE);
        }

        const int nperblock=n/nblocks; // if there is a rest, it is ignored
        const double norm = 1./nperblock;

        double av[nblocks]; // nblocks wont be huge number, so we can just stack-allocate
        for (int i1=0; i1<nblocks; ++i1) {
            av[i1] = accumulate(x+i1*nperblock, x+(i1+1)*nperblock, 0.);
            av[i1] *= norm;
        }

        UncorrelatedEstimator(nblocks, av, average, error);
    }


    void CorrelatedEstimator(const int &n, const double * x, double & average, double & error)
    {
        const int MIN_BLOCKS=6, MAX_BLOCKS=50;
        const int MAX_PLATEAU_AVERAGE=4;

        using namespace std;

        if ( n < MAX_BLOCKS) {
            cout << "MCI error CorrelatedEstimator() : n must be >= " << MAX_BLOCKS;
            exit(EXIT_FAILURE);
        }

        const int nav = MAX_BLOCKS-MIN_BLOCKS+1;
        double av[nav];
        double err[nav];
        for (int i1=0; i1<nav; ++i1) {
            const int nblocks=i1+MIN_BLOCKS;
            BlockEstimator(n, x, nblocks, av[i1], err[i1]);
            //cout << "Nblocks = " << nblocks << "   average = " << av[i1] << "   error = " << err[i1] << endl;
        }

        const int naccd = nav-2*MAX_PLATEAU_AVERAGE;
        double accdelta[naccd];
        fill(accdelta, accdelta+naccd, 0.);
        for (int i2=MAX_PLATEAU_AVERAGE; i2<naccd+MAX_PLATEAU_AVERAGE; ++i2) {
            for (int i1=1; i1<=MAX_PLATEAU_AVERAGE; ++i1) {
                accdelta[i2-MAX_PLATEAU_AVERAGE] += calcErrDelta(i1, &err[i2-4]);
            }
        }

        /*for (int i2=MAX_PLATEAU_AVERAGE; i2<MAX_BLOCKS-MIN_BLOCKS+1-MAX_PLATEAU_AVERAGE; ++i2) {
          cout << "i = " << i2+MIN_BLOCKS << "   accdelta = " << accdelta[i2-MAX_PLATEAU_AVERAGE] << endl;
          }*/

        int i_min = 0;
        for (int i2=1; i2<naccd; ++i2) {
            if ( fabs(accdelta[i2]) < fabs(accdelta[i_min]) ) { i_min=i2; }
        }

        i_min+=MAX_PLATEAU_AVERAGE;
        //cout << "The plateau has been detected at nblocks = " << i_min+MIN_BLOCKS << endl;
        average = 0.2*( av[i_min-2] + av[i_min-1] + av[i_min] + av[i_min+1] + av[i_min+2] );
        error = 0.2*( err[i_min-2] + err[i_min-1] + err[i_min] + err[i_min+1] + err[i_min+2] );
    }


    void MultiDimUncorrelatedEstimator(const int &n, const int &ndim, const double * x, double * average, double * error)
    {   // we create an explicit multidimensional implementation, for better efficiency in memory access
        using namespace std;

        const double SMALLEST_ERROR=1.e-300;

        if ( n < 2) {
            cout << "MCI error MultiDimUncorrelatedEstimator() : n must be larger than 1";
            exit(EXIT_FAILURE);
        }

        fill(average, average+ndim, 0.);
        fill(error, error+ndim, 0.);

        for (int i=0; i<n; ++i) {
            for (int j=0; j<ndim; ++j) {
                average[j] += x[i*ndim+j];
                error[j] += x[i*ndim+j]*x[i*ndim+j];
            }
        }

        const double norm = 1./n;
        const double norm2 = 1./(n-1.);
        for (int j=0; j<ndim; ++j) {
            average[j]*=norm;
            error[j]=error[j]*norm - average[j]*average[j];
            if ( error[j] > SMALLEST_ERROR ) {
                error[j]=sqrt( error[j]*norm2 );
            } else {
                error[j]=0.;
            }
        }
    }


    void MultiDimBlockEstimator(const int &n, const int &ndim, const double * x, const int &nblocks, double * average, double * error)
    {   // we create an explicit multidimensional implementation, for better efficiency in memory access
        using namespace std;

        if ( n < nblocks) {
            cout << "MCI error MultiDimBlockEstimator() : n must be >= nblocks";
            exit(EXIT_FAILURE);
        }

        const int nperblock=n/nblocks; // if there is a rest, it is ignored
        const double norm = 1./nperblock;
        const int ndata = nblocks*ndim;

        auto * av = new double[ndata]; // let's play it safe and heap-allocate
        fill(av, av+ndata, 0.);

        for (int i1=0; i1<nblocks; ++i1) {
            for (int i2=i1*nperblock; i2<(i1+1)*nperblock; ++i2) {
                for (int j=0; j<ndim; ++j) {
                    av[i1*ndim+j] += x[i2*ndim+j];
                }
            }
        }
        for (int i=0; i<ndata; ++i) {
            av[i] *= norm;
        }

        MultiDimUncorrelatedEstimator(nblocks, ndim, av, average, error);

        delete [] av;
    }


    void MultiDimCorrelatedEstimator(const int &n, const int &ndim, const double * x, double * average, double * error)
    {   // we create an explicit multidimensional implementation, for better efficiency in memory access
        const int MIN_BLOCKS=6, MAX_BLOCKS=50;
        const int MAX_PLATEAU_AVERAGE=4;

        using namespace std;

        if ( n < MAX_BLOCKS) {
            cout << "MCI error MultiDimCorrelatedEstimator() : n must be >= " << MAX_BLOCKS;
            exit(EXIT_FAILURE);
        }

        const int nav = MAX_BLOCKS-MIN_BLOCKS+1;
        const int nav_total = nav*ndim;
        auto * av = new double[nav_total];
        auto * err = new double[nav_total];

        for (int i1=0; i1<nav; ++i1) {
            const int nblocks=i1+MIN_BLOCKS;
            MultiDimBlockEstimator(n, ndim, x, nblocks, &av[i1*ndim], &err[i1*ndim]);
        }

        double delta[ndim];
        fill(delta, delta+ndim, 0.);

        const int naccd = nav-2*MAX_PLATEAU_AVERAGE;
        const int naccd_total = naccd*ndim;
        auto * accdelta = new double[naccd_total];
        fill(accdelta, accdelta+naccd_total, 0.);

        double errh[9]; // unfortunately we need to copy some values for passing
        for (int i2=MAX_PLATEAU_AVERAGE; i2<naccd+MAX_PLATEAU_AVERAGE; ++i2) {
            for (int j=0; j<ndim; ++j) {
                for (int ishift=-4; ishift<5; ++ishift) { errh[ishift+4] = err[(i2+ishift)*ndim+j]; }
                for (int i1=1; i1<=MAX_PLATEAU_AVERAGE; ++i1) {
                    accdelta[(i2-MAX_PLATEAU_AVERAGE)*ndim + j] += calcErrDelta(i1, errh);
                }
            }
        }

        int i_min[ndim];
        fill(i_min, i_min+ndim, 0);
        for (int i2=1; i2<MAX_BLOCKS-MIN_BLOCKS+1-2*MAX_PLATEAU_AVERAGE; ++i2) {
            for (int j=0; j<ndim; ++j) {
                if ( fabs( accdelta[i2*ndim+j] ) < fabs( accdelta[i_min[j]*ndim+j] ) ) {
                    i_min[j]=i2;
                }
            }
        }
        for (int j=0; j<ndim; ++j){
            const int index = i_min[j] + MAX_PLATEAU_AVERAGE;
            average[j] = 0.2*( av[(index-2)*ndim+j] + av[(index-1)*ndim+j] + av[index*ndim+j] + av[(index+1)*ndim+j] + av[(index+2)*ndim+j] );
            error[j] = 0.2*( err[(index-2)*ndim+j] + err[(index-1)*ndim+j] + err[index*ndim+j] + err[(index+1)*ndim+j] + err[(index+2)*ndim+j] );
        }

        delete [] accdelta;
        delete [] err;
        delete [] av;
    }

}  // namespace mci
