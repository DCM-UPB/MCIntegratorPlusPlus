#ifndef MCI_MJBLOCKER_HPP
#define MCI_MJBLOCKER_HPP

#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace mci
{

class MJBlocker
    // Custom implementation of the automatic blocking technique published by Marius Jonsson (PhysRevE 98, 043304 (2018)).
    // Implemented by Jan Kessler, by adapting the original ( https://github.com/computative/block/tree/master/c%2B%2B )
    // to better suit our specific needs. Notably, this means that the method was modified to allow
    // processing multi-dimensional data in one pass, with data passed as a flat C-style array. The data array
    // is expected to consist of ndata blocks of ndim doubles, just like what is produced by MCI observables.
    //
    // NOTE 1: Just like the original, this algorithm requires the number of samples, ndata, to be a power of 2.
    // NOTE 2: All required intermediate arrays are allocated to const pointers on object creation and the memory
    //         is not freed until object deletion. This means that a sequence of data arrays with the same layout
    //         can be processed without any reallocation in between.
    //
{
public:
    // --- Public Consts
    const int64_t ndata; // the number of samples, must be power of 2!
    const int ndim; // number of dimensions per sample
    const int npow; // number of powers of 2 to go through

private:
    // --- Internals

    // arrays storing variations of input x (length ndata*ndim)
    double * const _x; // first copy of input x
    double * const _X; // second copy of input x

    // smaller helper arrays to be pre-allocated (length npow*ndim)
    double * const _var; // for results of _gamma0() (variance)
    double * const _gamma; // for results of gamma1()
    double * const _M; // test statistics Mi
    double * const _Msum; // cumulative sum of test statistics Mi


    // Init
    void _computeMean(double mean[]) const; // compute multi-dimensional mean of input x
    void _initX(const double mean[]); // set internal _X = x - mean

    // Operations involving the copied data
    void _gamma0(double var[]/*part of _var*/, int64_t nred); // compute variance (gamma_h(0))
    void _gamma1(double gamma[]/*part of _gamma*/, int64_t nred); // compute gamma_h(1)
    int64_t _transform(const double mean[], int64_t nred); // perform blocking transform on _x and _X, return nred/2

    // Post-Processing
    void _generateM(); // performs cumulative sum of _M (and store in _M)

public:
    // --- User

    MJBlocker(int64_t n_data, int n_dim); // the constructor checks if ndata is power of 2 (else throw)
    ~MJBlocker();

    // Estimates average and variance of data array x, containing samples with dimension ndim
    void estimate(const double x[], double avg[], double err[]); // arrays have flat layout (ndata*ndim)
};
} // namespace mci

#endif
