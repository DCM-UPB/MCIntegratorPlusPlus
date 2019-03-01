#ifndef MCI_ESTIMATORS_HPP
#define MCI_ESTIMATORS_HPP

namespace mci
{
    // Compute average and standard deviation (error) of a set of data x[N], assuming that they are not correlated
    void OneDimUncorrelatedEstimator(const int &n, const double * x, double & average, double & error);

    // Compute average and error, using the blocking technique
    void OneDimBlockEstimator(const int &n, const double * x, const int &nblocks, double & average, double & error);

    // Compute average and error for correlated data, using the blocking technique
    void OneDimCorrelatedEstimator(const int &n, const double * x, double & average, double & error);


    // Estimators for multidimensional observable data
    // Useful when you have a n-sample data set of arrays with ndim elements
    //
    // Data are supposed to be organized in this way:
    // double data[n*ndim], where data[0]...data[ndim-1] is the first array, data[ndim]...data[2*ndim-1] the second, and so on.
    //
    // The output will be given by arrays: average[ndim] and error[ndim]

    // Compute average and standard deviation (error) of a set of data x[N], assuming that they are not correlated
    void MultiDimUncorrelatedEstimator(const int &n, const int &ndim, const double * x, double * average, double * error);

    // Compute average and error, using the blocking technique
    void MultiDimBlockEstimator(const int &n, const int &ndim, const double * x, const int &nblocks, double * average, double * error);

    // Compute average and error for correlated data, using the blocking technique
    void MultiDimCorrelatedEstimator(const int &n, const int &ndim, const double * x, double * average, double * error);


    // wrappers for any dim
    void UncorrelatedEstimator(const int &n, const int &ndim, const double * x, double * average, double * error);
    void BlockEstimator(const int &n, const int &ndim, const double * x, const int &nblocks, double * average, double * error);
    void CorrelatedEstimator(const int &n, const int &ndim, const double * x, double * average, double * error);

} // namespace mci

#endif
