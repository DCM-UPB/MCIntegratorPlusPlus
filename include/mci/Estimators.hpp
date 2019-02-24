#ifndef MCI_ESTIMATORS_HPP
#define MCI_ESTIMATORS_HPP

namespace mci
{
    // Compute average and standard deviation (error) of a set of data x[N], assuming that they are not correlated
    void UncorrelatedEstimator(const long &n, const double * x, double * average, double * error);

    // Compute average and error, using the blocking technique
    void BlockEstimator(const long &n, const double * x, const int &nblocks, double * average, double * error);

    // Compute average and error for correlated data, using the blocking technique
    void CorrelatedEstimator(const long &n, const double * x, double * average, double * error);


    // Multidimensional data set
    // Useful when you have an array of data sets
    // Data are supposed to be organized in this way:
    // data[n][ndim]
    // i.e. data[i][0] contains the element i for the first data set of which one wants to compute average and error
    //      data[i][1] for the second data set, etc.
    // The output will be given by arrays: average[ndim] and error[ndim]

    // Compute average and standard deviation (error) of a set of data x[N], assuming that they are not correlated
    void MultiDimUncorrelatedEstimator(const long &n, const int &ndim, const double * const * x, double * average, double *error);

    // Compute average and error, using the blocking technique
    void MultiDimBlockEstimator(const long &n, const int &ndim, const double * const * x, const int &nblocks, double * average, double * error);

    // Compute average and error for correlated data, using the blocking technique
    void MultiDimCorrelatedEstimator(const long &n, const int &ndim, const double * const * x, double * average, double * error);

} // namespace mci

#endif
