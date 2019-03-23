#ifndef MCI_FACTORIES_HPP
#define MCI_FACTORIES_HPP

#include "mci/AccumulatorInterface.hpp"
#include "mci/BlockAccumulator.hpp"
#include "mci/FullAccumulator.hpp"
#include "mci/SimpleAccumulator.hpp"

#include "mci/Estimators.hpp"

#include "mci/MultiStepMove.hpp"
#include "mci/SRRDAllMove.hpp"
#include "mci/SRRDVecMove.hpp"
#include "mci/TrialMoveInterface.hpp"

#include <algorithm>
#include <cstdint>
#include <functional>
#include <memory>
#include <stdexcept>

// All factory functions go here

namespace mci
{

    // --- Create Accumulators

    inline std::unique_ptr<AccumulatorInterface> createAccumulator(std::unique_ptr<ObservableFunctionInterface> obs, int blocksize = 1, int nskip = 1)
    {
        // sanity
        blocksize = std::max(0, blocksize);
        nskip = std::max(1, nskip);

        if (blocksize == 0) {
            return std::unique_ptr<AccumulatorInterface>( new SimpleAccumulator(std::move(obs), nskip) );
        }
        if (blocksize == 1) {
            return std::unique_ptr<AccumulatorInterface>( new FullAccumulator(std::move(obs), nskip) );
        }

        return std::unique_ptr<AccumulatorInterface>( new BlockAccumulator(std::move(obs), nskip, blocksize) );
    }



    // --- Create Estimator Functions

    enum class EstimatorType {
                              // Enumeration of built-in any-dim estimators with the same general interface,
                              // which are either not-blocking (to be used when data already consists of
                              // block-averages or is uncorrelated) or auto-blocking estimators.
                              Noop,
                              Uncorrelated,
                              Correlated, /* uses MJBlocker, if ndata power of 2, and else FCBlocker */
                              FCBlocker, /* Francesco's auto blocker implementation */
                              MJBlocker /* Our implementation of Marius Jonsson's auto blocking */
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

    inline std::function<void(int64_t/*nstore*/, int/*nobs*/, const double[]/*data*/, double[]/*avg*/, double[]/*error*/)> createEstimator(EstimatorType estimType /*from Estimators enumeration*/)
    {
        switch (estimType) {
        case EstimatorType::Noop:
            return NoopEstimator;

        case EstimatorType::Uncorrelated:
            return UncorrelatedEstimator;

        case EstimatorType::Correlated:
            return CorrelatedEstimator;

        case EstimatorType::FCBlocker:
            return FCBlockerEstimator;

        case EstimatorType::MJBlocker:
            return MJBlockerEstimator;

        default:
            throw std::domain_error("[createEstimator] Unhandled estimator enumerator.");
        }
    }

    inline std::function<void(int64_t, int, const double[], double[], double[])> createEstimator(const bool flag_correlated, const bool flag_error = true)
    {
        return createEstimator(selectEstimatorType(flag_correlated, flag_error));
    }



    // --- Create Trial Moves

    // Enumeration of move type classes
    enum class MoveType {
                         All,
                         Vec,
                         MultiStep
    };
    static constexpr std::initializer_list<MoveType> list_all_MoveType = {MoveType::All,
                                                                          MoveType::Vec,
                                                                          MoveType::MultiStep};

    // Enumeration of usable symmetric real valued random distribution
    enum class SRRDType {
                         // inherently symmetric distributions
                         Uniform,
                         Gaussian,
                         Student,
                         Cauchy,
                         // symmetrized distributions
                         Exponential,
                         Gamma,
                         Weibull,
                         Lognormal,
                         Chisq,
                         Fisher
    };
    static constexpr std::initializer_list<SRRDType> list_all_SRRDType = {SRRDType::Uniform,
                                                                          SRRDType::Gaussian,
                                                                          SRRDType::Student,
                                                                          SRRDType::Cauchy,
                                                                          SRRDType::Exponential,
                                                                          SRRDType::Gamma,
                                                                          SRRDType::Weibull,
                                                                          SRRDType::Lognormal,
                                                                          SRRDType::Chisq,
                                                                          SRRDType::Fisher};

    static constexpr double DEFAULT_MRT2STEP = 0.05; // step size default to fall-back to

    // common sanity check
    inline void checkTrialMoveSanity(int ndim, int ntypes = 1, const int typeEnds[] = nullptr)
    {
        if (ndim<1) { throw std::invalid_argument("[checkTrialMoveSanity] ndim must be at least 1."); }
        if (ntypes > 1 && typeEnds == nullptr) { throw std::invalid_argument("[checkTrialMoveSanity] ntypes>1 requires passed typeEnds array."); }
    }

    // create a all-index move with chosen SRRDType
    inline std::unique_ptr<TrialMoveInterface> createSRRDAllMove(
                                                                 SRRDType srrd /*from enum*/, int ndim,
                                                                 int ntypes = 1, const int typeEnds[] = nullptr /*needs to be passed if ntypes > 1*/
                                                                 )
    {
        // some sanity
        ntypes = std::max(1, ntypes);
        checkTrialMoveSanity(ndim, ntypes, typeEnds);

        // create chosen move
        switch (srrd) {
        case (SRRDType::Uniform):
            return std::unique_ptr<TrialMoveInterface>( new UniformAllMove(ndim, ntypes, typeEnds, DEFAULT_MRT2STEP) );
        case (SRRDType::Gaussian):
            return std::unique_ptr<TrialMoveInterface>( new GaussianAllMove(ndim, ntypes, typeEnds, DEFAULT_MRT2STEP) );
        case (SRRDType::Student):
            return std::unique_ptr<TrialMoveInterface>( new StudentAllMove(ndim, ntypes, typeEnds, DEFAULT_MRT2STEP) );
        case (SRRDType::Cauchy):
            return std::unique_ptr<TrialMoveInterface>( new CauchyAllMove(ndim, ntypes, typeEnds, DEFAULT_MRT2STEP) );
        case (SRRDType::Exponential):
            return std::unique_ptr<TrialMoveInterface>( new ExponentialAllMove(ndim, ntypes, typeEnds, DEFAULT_MRT2STEP) );
        case (SRRDType::Gamma):
            return std::unique_ptr<TrialMoveInterface>( new GammaAllMove(ndim, ntypes, typeEnds, DEFAULT_MRT2STEP) );
        case (SRRDType::Weibull):
            return std::unique_ptr<TrialMoveInterface>( new WeibullAllMove(ndim, ntypes, typeEnds, DEFAULT_MRT2STEP) );
        case (SRRDType::Lognormal):
            return std::unique_ptr<TrialMoveInterface>( new LognormalAllMove(ndim, ntypes, typeEnds, DEFAULT_MRT2STEP) );
        case (SRRDType::Chisq):
            return std::unique_ptr<TrialMoveInterface>( new ChisqAllMove(ndim, ntypes, typeEnds, DEFAULT_MRT2STEP) );
        case (SRRDType::Fisher):
            return std::unique_ptr<TrialMoveInterface>( new FisherAllMove(ndim, ntypes, typeEnds, DEFAULT_MRT2STEP) );

        default:
            throw std::domain_error("[createSRRDAllMove] Unhandled SRRDType enumerator.");
        }
    }

    // create a single-vector move with chosen SRRDType
    inline std::unique_ptr<TrialMoveInterface> createSRRDVecMove(
                                                                 SRRDType srrd /*from enum*/, int nvecs, /*number of vectors*/
                                                                 int veclen = 1 /*vec dim*/, int ntypes = 1, const int typeEnds[] = nullptr
                                                                 )
    {
        // some sanity
        veclen = std::max(1, veclen);
        ntypes = std::max(1, ntypes);
        checkTrialMoveSanity(nvecs*veclen, ntypes, typeEnds);

        // create chosen move
        switch (srrd) {
        case (SRRDType::Uniform):
            return std::unique_ptr<TrialMoveInterface>( new UniformVecMove(nvecs, veclen, ntypes, typeEnds, DEFAULT_MRT2STEP) );
        case (SRRDType::Gaussian):
            return std::unique_ptr<TrialMoveInterface>( new GaussianVecMove(nvecs, veclen, ntypes, typeEnds, DEFAULT_MRT2STEP) );
        case (SRRDType::Student):
            return std::unique_ptr<TrialMoveInterface>( new StudentVecMove(nvecs, veclen, ntypes, typeEnds, DEFAULT_MRT2STEP) );
        case (SRRDType::Cauchy):
            return std::unique_ptr<TrialMoveInterface>( new CauchyVecMove(nvecs, veclen, ntypes, typeEnds, DEFAULT_MRT2STEP) );
        case (SRRDType::Exponential):
            return std::unique_ptr<TrialMoveInterface>( new ExponentialVecMove(nvecs, veclen, ntypes, typeEnds, DEFAULT_MRT2STEP) );
        case (SRRDType::Gamma):
            return std::unique_ptr<TrialMoveInterface>( new GammaVecMove(nvecs, veclen, ntypes, typeEnds, DEFAULT_MRT2STEP) );
        case (SRRDType::Weibull):
            return std::unique_ptr<TrialMoveInterface>( new WeibullVecMove(nvecs, veclen, ntypes, typeEnds, DEFAULT_MRT2STEP) );
        case (SRRDType::Lognormal):
            return std::unique_ptr<TrialMoveInterface>( new LognormalVecMove(nvecs, veclen, ntypes, typeEnds, DEFAULT_MRT2STEP) );
        case (SRRDType::Chisq):
            return std::unique_ptr<TrialMoveInterface>( new ChisqVecMove(nvecs, veclen, ntypes, typeEnds, DEFAULT_MRT2STEP) );
        case (SRRDType::Fisher):
            return std::unique_ptr<TrialMoveInterface>( new FisherVecMove(nvecs, veclen, ntypes, typeEnds, DEFAULT_MRT2STEP) );

        default:
            throw std::domain_error("[createSRRDVecMove] Unhandled SRRDType enumerator.");
        }
    }

    // create a default version of selected move type
    inline std::unique_ptr<TrialMoveInterface> createMoveDefault(MoveType mtype, int ndim)
    {
        switch (mtype) {
        case (MoveType::All):
            return createSRRDAllMove(SRRDType::Uniform, ndim);
        case (MoveType::Vec):
            return createSRRDVecMove(SRRDType::Uniform, ndim); // default to single index moves
        case (MoveType::MultiStep):
            return std::unique_ptr<TrialMoveInterface>( new MultiStepMove(ndim) ); // contains no pdf per default, so that should be set before use

        default:
            throw std::domain_error("[createMoveDefault] Unhandled MoveType enumerator.");
        }
    }

} // namespace mci

#endif
