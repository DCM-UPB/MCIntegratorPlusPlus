#ifndef MCI_ESTIMATORS_HPP
#define MCI_ESTIMATORS_HPP

namespace mci
{
    // Compute average and standard deviation (error) of a set of data x[N], assuming that they are not correlated
    void OneDimUncorrelatedEstimator(int n, const double * x, double & average, double & error);

    // Compute average and error, using the blocking technique
    void OneDimBlockEstimator(int n, const double * x, int nblocks, double & average, double & error);

    // Compute average and error for correlated data, using the blocking technique
    void OneDimCorrelatedEstimator(int n, const double * x, double & average, double & error);


    // Estimators for multidimensional observable data
    // Useful when you have a n-sample data set of arrays with ndim elements
    //
    // Data are supposed to be organized in this way:
    // double data[n*ndim], where data[0]...data[ndim-1] is the first array, data[ndim]...data[2*ndim-1] the second, and so on.
    //
    // The output will be given by arrays: average[ndim] and error[ndim]

    // Compute average and standard deviation (error) of a set of data x[N], assuming that they are not correlated
    void MultiDimUncorrelatedEstimator(int n, int ndim, const double * x, double * average, double * error);

    // Compute average and error, using the blocking technique
    void MultiDimBlockEstimator(int n, int ndim, const double * x, int nblocks, double * average, double * error);

    // Compute average and error for correlated data, using the blocking technique
    void MultiDimCorrelatedEstimator(int n, int ndim, const double * x, double * average, double * error);


    // wrappers for any dim
    void UncorrelatedEstimator(int n, int ndim, const double * x, double * average, double * error);
    void BlockEstimator(int n, int ndim, const double * x, int nblocks, double * average, double * error);
    void CorrelatedEstimator(int n, int ndim, const double * x, double * average, double * error);

} // namespace mci

#endif
