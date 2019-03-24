#include "mci/MJBlocker.hpp"

#include <algorithm>

// --- A static array of magic numbers
static constexpr std::array<double, 64> quantile{3.841459, 5.991465, 7.814728, 9.487729, 11.070498, 12.591587, 14.067140, 15.507313,
                                                 16.918978, 18.307038, 19.675138, 21.026070, 22.362032, 23.684791, 24.995790, 26.296228,
                                                 27.587112, 28.869299, 30.143527, 31.410433, 32.670573, 33.924438, 35.172462, 36.415029,
                                                 37.652484, 38.885139, 40.113272, 41.337138, 42.556968, 43.772972, 44.985343, 46.194260,
                                                 47.399884, 48.602367, 49.801850, 50.998460, 52.192320, 53.383541, 54.572228, 55.758479,
                                                 56.942387, 58.124038, 59.303512, 60.480887, 61.656233, 62.829620, 64.001112, 65.170769,
                                                 66.338649, 67.504807, 68.669294, 69.832160, 70.993453, 72.153216, 73.311493, 74.468324,
                                                 75.623748, 76.777803, 77.930524, 79.081944, 80.232098, 81.381015, 82.528727, 83.675261};

namespace mci
{
// --- Constructor/Destructor

MJBlocker::MJBlocker(const int64_t n_data, const int n_dim):
        ndata(n_data), ndim(n_dim), npow(static_cast<int>(log2(n_data))), // below we check that n is a power of two, so d is a (small) integer
        _x(new double[ndata*ndim]), _X(new double[ndata*ndim]),
        _var(new double[npow*ndim]), _gamma(new double[npow*ndim]),
        _M(new double[npow*ndim]), _Msum(new double[npow*ndim])
{
    if (ndata <= 0) {
        throw std::invalid_argument("[MJBlocker] ndata must be a natural number.");
    }
    if ((static_cast<uint64_t>(ndata) & (static_cast<uint64_t>(ndata) - 1)) != 0) { // ndata must be power of two
        throw std::invalid_argument("[MJBlocker] ndata must be a power of two.");
    }
}

MJBlocker::~MJBlocker()
{
    delete[] _Msum;
    delete[] _M;
    delete[] _gamma;
    delete[] _var;
    delete[] _X;
    delete[] _x;
}

// --- Initializers

// estimates mean of x
void MJBlocker::_computeMean(double mean[]) const
{
    std::fill(mean, mean + ndim, 0.);
    for (int64_t i = 0; i < ndata; ++i) {
        for (int j = 0; j < ndim; ++j) {
            mean[j] += _x[i*ndim + j];
        }
    }
    for (int j = 0; j < ndim; ++j) { mean[j] /= ndata; }
}

// stores x minus mean
void MJBlocker::_initX(const double mean[])
{
    for (int64_t i = 0; i < ndata; ++i) {
        for (int j = 0; j < ndim; ++j) {
            _X[i*ndim + j] = _x[i*ndim + j] - mean[j];
        }
    }
}

// --- Operations involving the copies of input data

// estimates gamma_h(0) for all h
void MJBlocker::_gamma0(double var[], const int64_t nred)
{
    std::fill(var, var + ndim, 0.);
    for (int64_t i = 0; i < nred; ++i) {
        for (int j = 0; j < ndim; ++j) {
            var[j] += _X[i*ndim + j]*_X[i*ndim + j];
        }
    }
    for (int j = 0; j < ndim; ++j) { var[j] /= nred; }
}

// estimates gamma_h(1) for all h
void MJBlocker::_gamma1(double gamma[], const int64_t nred)
{
    std::fill(gamma, gamma + ndim, 0.);
    for (int64_t i = 0; i < nred - 1; ++i) {
        for (int j = 0; j < ndim; ++j) {
            gamma[j] += _X[i*ndim + j]*_X[(i + 1)*ndim + j];
        }
    }
    for (int j = 0; j < ndim; ++j) { gamma[j] /= nred; }
}

// performs blocking transformation
int64_t MJBlocker::_transform(const double mean[], const int64_t nred)
{
    for (int64_t i = 0; i < nred/2; ++i) {
        for (int j = 0; j < ndim; ++j) {
            _x[i*ndim + j] = 0.5*(_x[(2*i)*ndim + j] + _x[(2*i + 1)*ndim + j]);
            _X[i*ndim + j] = _x[i*ndim + j] - mean[j];
        }
    }
    return nred/2;
}

// generate test statistics Mk (in reverse order) and perform cumulative sum (stored in _Msum)
void MJBlocker::_generateM()
{
    for (int i = 0; i < npow; ++i) {
        for (int k = 0; k < ndim; ++k) {
            _M[(npow - i - 1)*ndim + k] = pow(_gamma[i*ndim + k]/_var[i*ndim + k], 2)*pow(2, npow - i);
        }
    }

    // compute cumulative sum
    std::fill(_Msum, _Msum + npow*ndim, 0.);
    for (int i = 0; i < npow; ++i) {
        for (int j = 0; j <= i; ++j) {
            for (int k = 0; k < ndim; ++k) {
                _Msum[i*ndim + k] += _M[j*ndim + k];
            }
        }
    }
}


// the algorithm which computes the variance of the sample mean.
void MJBlocker::estimate(const double x[], double avg[], double err[])
{
    std::copy(x, x + ndata*ndim, _x); // first copy our input data ( we don't want to modify the input )
    this->_computeMean(avg); // store average of x in avg
    this->_initX(avg); // store x minus avg in _X

    // compute covariance and variance and apply blocking transform
    int64_t nred = ndata; // will go through powers of 2
    for (int k = 0; k < npow; ++k) {
        this->_gamma0(_var + k*ndim, nred);
        this->_gamma1(_gamma + k*ndim, nred);
        nred = this->_transform(avg, nred);
    }

    this->_generateM(); // compute cumulative test statistics sum
    for (int l = 0; l < ndim; ++l) { // find the first k such that Mk < quantile[k]
        int k;
        for (k = npow - 1; k >= 0; --k) {
            if (_Msum[k*ndim + l] < quantile[k]) {
                break;
            }
        }
        k = npow - (k + 1);

        // and finally compute the errors
        err[l] = sqrt(_var[k*ndim + l]/pow(2, npow - k));
    }
}
} // namespace mci
