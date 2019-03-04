#ifndef MCI_MCIOBSERVABLECONTAINER_HPP
#define MCI_MCIOBSERVABLECONTAINER_HPP

#include "mci/MCIAccumulatorInterface.hpp"

#include <functional>
#include <stdexcept>

struct MCIObservableContainer
{
    // Accumulator
    MCIAccumulatorInterface * const accu;

    // Estimator function used to obtain result of MC integration
    const std::function< void (double * /*avg*/, double * /*error*/) > estim; // calculating avg/error of stored data, with all parameters bound already (e.g. nstored, nobs)

    // flags
    const bool flag_equil; // should this observable be equilibrated when using automatic equilibration?


    MCIObservableContainer(MCIAccumulatorInterface * accumulator,
                           const std::function< void (int /*nstored*/, int /*nobs*/, const double * /*data*/, double * /*avg*/, double * /*error*/) > &estimator,
                           bool needsEquil):
        accu(accumulator),
        estim( [estimator, accumulator](double * average, double * error) { // lambda functional
                   if(!accumulator->isFinalized()) {
                       throw std::runtime_error("[MCIObservableContainer.estim] Estimator was called, but accumulator is not finalized.");
                   }
                   estimator(accumulator->getNStore(), accumulator->getNObs(), accumulator->getData(), average, error);
               } ),
        flag_equil(needsEquil)
    {}
};


#endif
