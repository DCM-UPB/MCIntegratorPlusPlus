#if USE_MPI==1

#ifndef MPIMCI
#define MPIMCI

#include <stdexcept>
#include <string>
#include <fstream>
#include <mpi.h>
#include "mci/MCIntegrator.hpp"

namespace MPIMCI
{
    // return my rank
    int rank()
    {
        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        return myrank;
    }

    // return size
    int size()
    {
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        return size;
    }

    // init MPI and return rank of process
    int init()
    {
        int isinit;
        MPI_Initialized(&isinit);
        if (isinit==1) throw std::runtime_error("MPI already initialized!");
        MPI_Init(NULL, NULL);
        return rank();
    }

    // set different random seeds per thread from a file
    void setSeed(MCI * const mci, const std::string &filename, const int &offset = 0) // with offset you can control how many seeds to skip initially
    {
        int myrank = rank();
        int nranks = size();

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


    // integrate in parallel and accumulate results

    void integrate(MCI * const mci, const long &Nmc, double * average, double * error, int NfindMRT2stepIterations, int NdecorrelationSteps, size_t nblocks=0, bool use_mpi=true) // by setting use_mpi to false you can use this without requiring MPI
    {
        if (use_mpi) {
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

            mci->integrate(Nmc, myAverage, myError, NfindMRT2stepIterations, NdecorrelationSteps, nblocks);

            for (int i=0; i<mci->getNObsDim(); ++i) {
                MPI_Allreduce(&myAverage[i], &average[i], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                myError[i] *= myError[i];
                MPI_Allreduce(&myError[i], &error[i], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                average[i] /= nranks;
                error[i] = sqrt(error[i]) / nranks;
            }
        }
        else {
            mci->integrate(Nmc, average, error, NfindMRT2stepIterations, NdecorrelationSteps, nblocks); // regular single thread call
        }
    }

    void integrate(MCI * const mci, const long &Nmc, double * average, double * error, bool findMRT2step=true, bool initialdecorrelation=true, size_t nblocks=0, bool use_mpi=true) // auto-mode-wrapper
    {
        int stepsMRT2 = findMRT2step ? -1 : 0;
        int stepsDecorr = initialdecorrelation ? -1 : 0;
        integrate(mci, Nmc, average, error, stepsMRT2, stepsDecorr, nblocks, use_mpi);
    }


    // finalize MPI
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
};

#endif

#endif
