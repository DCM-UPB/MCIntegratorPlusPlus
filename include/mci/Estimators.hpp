#ifndef MCI_ESTIMATORS_HPP
#define MCI_ESTIMATORS_HPP

#include <cstdint>

namespace mci
{
    // Compute average and standard deviation (error) of a set of data x[N], assuming that they are not correlated
    void OneDimUncorrelatedEstimator(int64_t n, const double x[], double & average, double & error);

    // Compute average and error, using the blocking technique (only used by FCBlockerEstimator, not by MCI)
    void OneDimBlockEstimator(int64_t n, const double x[], int64_t nblocks, double & average, double & error);

    // Compute average and error for correlated data, using auto blocking technique (by Francesco Calcavecchia)
    void OneDimFCBlockerEstimator(int64_t n, const double x[], double & average, double & error);


    // Estimators for multidimensional observable data
    // Useful when you have a n-sample data set of arrays with ndim elements
    //
    // Data are supposed to be organized in this way:
    // double data[n*ndim], where data[0]...data[ndim-1] is the first array, data[ndim]...data[2*ndim-1] the second, and so on.
    //
    // The output will be given by arrays: average[ndim] and error[ndim]

    // Compute average and standard deviation (error) of a set of data x[N], assuming that they are not correlated
    void MultiDimUncorrelatedEstimator(int64_t n, int ndim, const double x[], double average[], double error[]);

    // Compute average and error, using fixed blocking technique (only used by FCBlockerEstimator, not by MCI)
    void MultiDimBlockEstimator(int64_t n, int ndim, const double x[], int64_t nblocks, double average[], double error[]);

    // Compute average and error for correlated data, using auto blocking technique (by Francesco Calcavecchia)
    void MultiDimFCBlockerEstimator(int64_t n, int ndim, const double x[], double average[], double error[]);


    // any-dim wrappers for above functions and other estimators
    void UncorrelatedEstimator(int64_t n, int ndim, const double x[], double average[], double error[]);
    void CorrelatedEstimator(int64_t n, int ndim, const double x[], double average[], double error[]); // if n is power of 2, use MJBlocker, else FCBlocker
    void FCBlockerEstimator(int64_t n, int ndim, const double x[], double average[], double error[]);

    // Calls our implementation Marius Jonsson's auto-blocking algorithm
    void MJBlockerEstimator(int64_t n, int ndim, const double x[], double average[], double error[]);

    // no-op estimator (used when data contains the averages already and error is irrelevant)
    void NoopEstimator(int64_t/*n*/, int ndim, const double x[], double average[], double error[]);

} // namespace mci

#endif
