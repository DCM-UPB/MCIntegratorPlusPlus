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

    void integrate(MCI * const mci, const long &Nmc, double * average, double * error)
    {
        // make sure the user has MPI in the correct state
        if (!MPI::Is_initialized()) throw std::runtime_error("MPI not initialized!");
        if (MPI::Is_finalized()) throw std::runtime_error("MPI already finalized!");

        const int myrank = MPI::COMM_WORLD.Get_rank();
        const int nranks = MPI::COMM_WORLD.Get_size();

        // the results are stored in myAverage/Error and then reduced into average/error for root process
        double myAverage[mci->getNObsDim()];
        double myError[mci->getNObsDim()];

        mci->integrate(Nmc, myAverage, myError);

        for (int i=0; i<mci->getNObsDim(); ++i) {
            MPI::COMM_WORLD.Reduce(&myAverage[i], &average[i], 1, MPI::DOUBLE, MPI::SUM, 0);

            myError[i] *= myError[i];
            MPI::COMM_WORLD.Reduce(&myError[i], &error[i], 1, MPI::DOUBLE, MPI::SUM, 0);

            if (myrank == 0) {
                average[i] /= nranks;
                error[i] = sqrt(error[i]) / nranks;
            }
        }
    }

    void finalize(MCI * mci)
    {
        // make sure the user has MPI in the correct state
        if (!MPI::Is_initialized()) throw std::runtime_error("MPI not initialized!");
        if (MPI::Is_finalized()) throw std::runtime_error("MPI already finalized!");

        delete mci;
        MPI::Finalize();
    }
};
