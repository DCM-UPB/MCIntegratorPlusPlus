#include "mci/ObservableContainer.hpp"

#include <stdexcept>

namespace mci
{

    void ObservableContainer::addObservable(AccumulatorInterface * accumulator,
                                            const std::function< void (int, int, const double [], double [], double []) > &estimator,
                                            bool needsEquil)
    {
        _accus.emplace_back(accumulator);
        _nobsdim+=accumulator->getNObs();

        _estims.emplace_back( [accumulator, estimator](double average[], double error[]) { // lambda functional
                                  if(!accumulator->isFinalized()) {
                                      throw std::runtime_error("[ObservableContainer.estim] Estimator was called, but accumulator is not finalized.");
                                  }
                                  estimator(accumulator->getNStore(), accumulator->getNObs(), accumulator->getData(), average, error);
                              } );
        _flags_equil.emplace_back(needsEquil ? 1 : 0);
    }


    void ObservableContainer::allocate(const int Nmc)
    {
        for (auto & accu : _accus) {
            accu->allocate(Nmc);
        }
    }


    void ObservableContainer::accumulate(const double x[], const bool flagacc, const bool flags_xchanged[])
    {
        for (auto & accu : _accus) {
            accu->accumulate(x, flagacc, flags_xchanged);
        }
    }


    void ObservableContainer::printObsValues(std::ofstream &file)
    {
        for (auto & accu : _accus) {
            ObservableFunctionInterface * const obs = accu->getObservableFunction(); // acquire ptr to obsfun
            for (int j=0; j<obs->getNObs(); ++j) {
                file << " " << obs->getValue(j);
            }
        }
        file << " ";
    }


    void ObservableContainer::finalize()
    {
        for (auto & accu : _accus) {
            accu->finalize();
        }
    }


    void ObservableContainer::estimate(double average[], double error[])
    {
        int iobs = 0;
        int offset = 0;
        for (auto & estim : _estims) { // use all estimators
            estim(average+offset, error+offset);
            offset += _accus[iobs]->getNObs();
            ++iobs;
        }
    }


    void ObservableContainer::reset()
    {
        for (auto & accu : _accus) {
            accu->reset();
        }
    }

    void ObservableContainer::deallocate()
    {
        for (auto & accu : _accus) {
            accu->deallocate();
        }
    }

    void ObservableContainer::clear()
    {
        for (auto & accu : _accus) {
            delete accu;
        }
        _accus.clear();
        _nobsdim=0;

        _estims.clear();
        _flags_equil.clear();
    }

}
