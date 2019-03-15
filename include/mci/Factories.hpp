#ifndef MCI_FACTORIES_HPP
#define MCI_FACTORIES_HPP

#include "mci/AccumulatorInterface.hpp"
#include "mci/BlockAccumulator.hpp"
#include "mci/FullAccumulator.hpp"
#include "mci/SimpleAccumulator.hpp"

#include "mci/Estimators.hpp"

#include <algorithm>
#include <functional>
#include <memory>
#include <stdexcept>

// All factory functions go here

namespace mci
{

    // --- Create Accumulators

    inline std::unique_ptr<AccumulatorInterface> createAccumulator(const ObservableFunctionInterface &obs, int blocksize = 1, int nskip = 1)
    {
        // sanity
        blocksize = std::max(0, blocksize);
        nskip = std::max(1, nskip);

        if (blocksize == 0) {
            return std::unique_ptr<AccumulatorInterface>( new SimpleAccumulator(obs.clone(), nskip) );
        }
        if (blocksize == 1) {
            return std::unique_ptr<AccumulatorInterface>( new FullAccumulator(obs.clone(), nskip) );
        }

        return std::unique_ptr<AccumulatorInterface>( new BlockAccumulator(obs.clone(), nskip, blocksize) );
    }


    // --- Create Estimator Functions

    enum class EstimatorType {
                              // Enumeration of built-in any-dim estimators with the same general interface,
                              // which are either not-blocking (to be used when data already consists of
                              // block-averages or is uncorrelated) or auto-blocking estimators.
                              Noop,
                              Uncorrelated,
                              Correlated
    };
    inline EstimatorType selectEstimatorType(const bool flag_correlated, const bool flag_error = true)
    {
        if (flag_correlated) {
            if (!flag_error) {
                throw std::invalid_argument("[selectEstimatorType] Error calculation is set off, but correlated error estimation is set on.");
            }
            return EstimatorType::Correlated;
        }

        return flag_error? EstimatorType::Uncorrelated : EstimatorType::Noop;
    }

    inline std::function<void(int/*nstore*/, int/*nobs*/, const double[]/*data*/, double[]/*avg*/, double[]/*error*/)> createEstimator(EstimatorType estimType /*from Estimators enumeration*/)
    {
        switch (estimType) {
        case EstimatorType::Noop:
            return NoopEstimator;

        case EstimatorType::Uncorrelated:
            return UncorrelatedEstimator;

        case EstimatorType::Correlated:
            return CorrelatedEstimator;

        default:
            throw std::domain_error("[createEstimator] Unhandled estimator enumerator.");
        }
    }
    inline std::function<void(int, int, const double[], double[], double[])> createEstimator(const bool flag_correlated, const bool flag_error = true)
    {
        return createEstimator(selectEstimatorType(flag_correlated, flag_error));
    }

} // namespace mci

#endif
