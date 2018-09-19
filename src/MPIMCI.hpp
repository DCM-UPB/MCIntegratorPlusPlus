#ifndef MPIMCI
#define MPIMCI


#include <mpi.h>
#include "MCIntegrator.hpp"

namespace MPIMCI
{
    int init() // return mpi rank of process
    {
        if (MPI::Is_initialized()) throw std::runtime_error("MPI already initialized!");
        MPI::Init();
        MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
        return MPI::COMM_WORLD.Get_rank();
    }

    void integrate(MCI * const mci, const long &Nmc, double * average, double * error, bool findMRT2step=true, bool initialdecorrelation=true, bool use_mpi=true) // by setting use_mpi to false you can use this without requiring MPI
    {
        if (use_mpi) {
            // make sure the user has MPI in the correct state
            if (!MPI::Is_initialized()) throw std::runtime_error("MPI not initialized!");
            if (MPI::Is_finalized()) throw std::runtime_error("MPI already finalized!");

            const int nranks = MPI::COMM_WORLD.Get_size();

            // the results are stored in myAverage/Error and then reduced into the same average/error for all processes
            double myAverage[mci->getNObsDim()];
            double myError[mci->getNObsDim()];

            mci->integrate(Nmc, myAverage, myError, findMRT2step, initialdecorrelation);

            for (int i=0; i<mci->getNObsDim(); ++i) {
                MPI::COMM_WORLD.Allreduce(&myAverage[i], &average[i], 1, MPI::DOUBLE, MPI::SUM);

                myError[i] *= myError[i];
                MPI::COMM_WORLD.Allreduce(&myError[i], &error[i], 1, MPI::DOUBLE, MPI::SUM);

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
        if (!MPI::Is_initialized()) throw std::runtime_error("MPI not initialized!");
        if (MPI::Is_finalized()) throw std::runtime_error("MPI already finalized!");

        MPI::Finalize();
    }
};

#endif
