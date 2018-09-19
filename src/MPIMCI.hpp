#ifndef MPIMCI
#define MPIMCI

#include <stdexcept>
#include <mpi.h>
#include "MCIntegrator.hpp"

namespace MPIMCI
{
    int init() // return mpi rank of process
    {
        int isinit;
        MPI_Initialized(&isinit);
        if (isinit==1) throw std::runtime_error("MPI already initialized!");
        MPI_Init(NULL, NULL);
        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        return myrank;
    }

    void integrate(MCI * const mci, const long &Nmc, double * average, double * error, bool findMRT2step=true, bool initialdecorrelation=true, bool use_mpi=true) // by setting use_mpi to false you can use this without requiring MPI
    {
        if (use_mpi) {
            // make sure the user has MPI in the correct state
            int isinit, isfinal;
            MPI_Initialized(&isinit);
            if (isinit==0) throw std::runtime_error("MPI not initialized!");
            MPI_Finalized(&isfinal);
            if (isfinal==1) throw std::runtime_error("MPI already finalized!");

            int nranks;
            MPI_Comm_size(MPI_COMM_WORLD, &nranks);

            // the results are stored in myAverage/Error and then reduced into the same average/error for all processes
            double myAverage[mci->getNObsDim()];
            double myError[mci->getNObsDim()];

            mci->integrate(Nmc, myAverage, myError, findMRT2step, initialdecorrelation);

            for (int i=0; i<mci->getNObsDim(); ++i) {
                MPI_Allreduce(&myAverage[i], &average[i], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                myError[i] *= myError[i];
                MPI_Allreduce(&myError[i], &error[i], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                average[i] /= nranks;
                error[i] = sqrt(error[i]) / nranks;
            }
        }
        else {
            mci->integrate(Nmc, average, error, findMRT2step, initialdecorrelation); // regular single thread call
        }
    }

    void finalize()
    {
        // make sure the user has MPI in the correct state
        int isinit, isfinal;
        MPI_Initialized(&isinit);
        if (isinit==0) throw std::runtime_error("MPI not initialized!");
        MPI_Finalized(&isfinal);
        if (isfinal==1) throw std::runtime_error("MPI already finalized!");

        MPI_Finalize();
    }
};

#endif
