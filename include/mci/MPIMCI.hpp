#ifndef MCI_MPIMCI_HPP
#define MCI_MPIMCI_HPP

#if USE_MPI==1

#include "mci/MCIntegrator.hpp"

#include <cstdint>
#include <string>

namespace MPIMCI
{
    // return my rank
    int myrank(); // we don't name it rank() here because it would collide with "using namespace std"

    // return size
    int size();

    // init MPI and return rank of process
    int init();

    // set different random seeds per thread from a file
    void setSeed(mci::MCI &mci, const std::string &filename, int offset = 0); // with offset you can control how many seeds to skip initially

    // integrate in parallel and accumulate results
    void integrate(mci::MCI &mci, int64_t Nmc, double average[], double error[], bool doFindMRT2Step = true, bool doDecorrelation = true);

    // finalize MPI
    void finalize();

} // namespace MPIMCI

#endif

#endif
