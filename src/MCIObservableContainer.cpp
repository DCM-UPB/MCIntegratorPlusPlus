#include "mci/MCIObservableContainer.hpp"

#include <stdexcept>

void MCIObservableContainer::addObservable(MCIAccumulatorInterface * accumulator,
                                      const std::function< void (int, int, const double *, double *, double *) > &estimator,
                                      bool needsEquil)
{
    _accus.emplace_back(accumulator);
    _nobsdim+=accumulator->getNObs();

    _estims.emplace_back( [accumulator, estimator](double * average, double * error) { // lambda functional
                              if(!accumulator->isFinalized()) {
                                  throw std::runtime_error("[MCIObservableContainer.estim] Estimator was called, but accumulator is not finalized.");
                              }
                              estimator(accumulator->getNStore(), accumulator->getNObs(), accumulator->getData(), average, error);
                          } );
    _flags_equil.emplace_back(needsEquil);
}


void MCIObservableContainer::allocateObservables(const int Nmc)
{   // allocate observable accumulators for Nmc steps
    for (auto & accu : _accus) {
        accu->allocate(Nmc);
    }
}


void MCIObservableContainer::accumulateObservables(const double * x, const bool flagacc)
{   // let the accumulators process the step
    for (auto & accu : _accus) {
        accu->accumulate(x, flagacc);
    }
}


void MCIObservableContainer::printObservableValues(std::ofstream &file)
{
    for (auto & accu : _accus) {
        MCIObservableFunctionInterface * const obs = accu->getObservableFunction(); // acquire ptr to obsfun
        for (int j=0; j<obs->getNObs(); ++j) {
            file << " " << obs->getValue(j);
        }
    }
    file << " ";
}


void MCIObservableContainer::finalizeObservables()
{   // apply normalization, if necessary
    for (auto & accu : _accus) {
        accu->finalize();
    }
}


void MCIObservableContainer::estimateObservables(double * average, double * error)
{
    int iobs = 0;
    int offset = 0;
    for (auto & estim : _estims) { // use all estimators
        estim(average+offset, error+offset);
        offset += _accus[iobs]->getNObs();
        ++iobs;
    }
}


void MCIObservableContainer::resetObservables()
{   // reset without deallocation
    for (auto & accu : _accus) {
        accu->reset();
    }
}

void MCIObservableContainer::deallocateObservables()
{   // reset & free memory
    for (auto & accu : _accus) {
        accu->deallocate();
    }
}

void MCIObservableContainer::clearObservables()
{
    for (auto & accu : _accus) {
        delete accu;
    }
    _accus.clear();
    _nobsdim=0;

    _estims.clear();
    _flags_equil.clear();
}
