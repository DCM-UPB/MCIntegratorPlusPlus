#if USE_MPI==1

#include "mci/MPIMCI.hpp"

#include "mci/MCIntegrator.hpp"

#include <mpi.h>
#include <stdexcept>
#include <string>

namespace MPIMCI
{
    int myrank()
    {
        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        return myrank;
    }


    int size()
    {
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        return size;
    }


    int init()
    {
        int isinit;
        MPI_Initialized(&isinit);
        if (isinit==1) { throw std::runtime_error("MPI already initialized!");}
        MPI_Init(nullptr, nullptr);

        return myrank();
    }


    void setSeed(MCI * mci, const std::string &filename, const int offset)
    {
        int myrank = MPIMCI::myrank();
        int nranks = MPIMCI::size();

        uint_fast64_t seeds[nranks];
        uint_fast64_t myseed;
        if (myrank == 0) {
            std::ifstream seedfile;
            seedfile.open(filename);

            if (!seedfile.good()) {throw std::runtime_error("Random seed file could not be found.");}

            for (int i=0; i<offset; ++i) {
                if (seedfile.eof()) {throw std::runtime_error("Chosen seed offset is already larger than the number of seeds in seed file.");}
                uint_fast64_t skip;
                seedfile >> skip;
            }
            for (int i=0; i<nranks; ++i) {
                if (seedfile.eof()) {throw std::runtime_error("Seed file doesn't provide enough seeds for the chosen number of ranks and offset.");}
                seedfile >> seeds[i];
            }
            seedfile.close();
        }

        MPI_Scatter(seeds, 1, MPI_UNSIGNED_LONG, &myseed, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        mci->setSeed(myseed);
    }


    void integrate(MCI * mci, int Nmc, double * average, double * error, const bool doFindMRT2Step, const bool doDecorrelation)
    {
        // make sure the user has MPI in the correct state
        int isinit, isfinal;
        MPI_Initialized(&isinit);
        if (isinit==0) {throw std::runtime_error("MPI not initialized!");}
        MPI_Finalized(&isfinal);
        if (isfinal==1) {throw std::runtime_error("MPI already finalized!");}

        int nranks = size();

        // the results are stored in myAverage/Error and then reduced into the same average/error for all processes
        double myAverage[mci->getNObsDim()];
        double myError[mci->getNObsDim()];

        mci->integrate(Nmc, myAverage, myError, doFindMRT2Step, doDecorrelation);

        for (int i=0; i<mci->getNObsDim(); ++i) {
            MPI_Allreduce(&myAverage[i], &average[i], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            myError[i] *= myError[i];
            MPI_Allreduce(&myError[i], &error[i], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            average[i] /= nranks;
            error[i] = sqrt(error[i]) / nranks;
        }
    }


    void finalize()
    {
        // make sure the user has MPI in the correct state
        int isinit, isfinal;
        MPI_Initialized(&isinit);
        if (isinit==0) {throw std::runtime_error("MPI not initialized!");}
        MPI_Finalized(&isfinal);
        if (isfinal==1) {throw std::runtime_error("MPI already finalized!");}

        MPI_Finalize();
    }
} // namespace MPIMCI

#endif
