#ifndef MCI_FACTORIES_HPP
#define MCI_FACTORIES_HPP

#include "mci/AccumulatorInterface.hpp"
#include "mci/BlockAccumulator.hpp"
#include "mci/FullAccumulator.hpp"
#include "mci/SimpleAccumulator.hpp"

#include "mci/Estimators.hpp"

#include <memory>
#include <algorithm>
#include <functional>
#include <stdexcept>

// All factory functions go here

namespace mci
{

    // --- Create Accumulators

    std::unique_ptr<AccumulatorInterface> createAccumulator(const ObservableFunctionInterface &obs, int blocksize = 1, int nskip = 1)
    {
        // sanity
        blocksize = std::max(0, blocksize);
        nskip = std::max(1, nskip);

        if (blocksize == 0) {
            return std::unique_ptr<AccumulatorInterface>( new SimpleAccumulator(obs.clone(), nskip) );
        }
        else {
            if (blocksize == 1) {
                return std::unique_ptr<AccumulatorInterface>( new FullAccumulator(obs.clone(), nskip) );
            } else {
                return std::unique_ptr<AccumulatorInterface>( new BlockAccumulator(obs.clone(), nskip, blocksize) );
            }
        }
    }


    // --- Create Estimator Functions

    enum class Estimator {
                          // Enumeration of built-in any-dim estimators with the same general interface,
                          // which are either not-blocking (to be used when data already consists of
                          // block-averages or is uncorrelated) or auto-blocking estimators.
                          Noop,
                          Uncorrelated,
                          Correlated
    };
    std::function<void(int/*nstore*/, int/*nobs*/, const double[]/*data*/, double[]/*avg*/, double[]/*error*/)> createEstimator(Estimator estim /*from Estimators enumeration*/)
    {
        switch (estim) {
        case Estimator::Noop:
            return NoopEstimator;

        case Estimator::Uncorrelated:
            return UncorrelatedEstimator;

        case Estimator::Correlated:
            return CorrelatedEstimator;

        default:
            throw std::domain_error("[createEstimator] Unhandled estimator enumerator.");
        }
    }
    std::function<void(int, int, const double[], double[], double[])> createEstimator(const bool flag_correlated, const bool flag_error = true)
    {
        if (flag_correlated) {
            if (!flag_error) {
                throw std::invalid_argument("[createEstimator] Error calculation was set off, but correlated error estimation was set on.");
            }
            return createEstimator(Estimator::Correlated);
        }
        else {
            return createEstimator(flag_error? Estimator::Uncorrelated : Estimator::Noop);
        }
    }

} // namespace mci

#endif
