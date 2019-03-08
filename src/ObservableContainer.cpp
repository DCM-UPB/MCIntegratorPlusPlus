#include "mci/ObservableContainer.hpp"

#include <stdexcept>

namespace mci
{

    void ObservableContainer::addObservable(std::unique_ptr<AccumulatorInterface> accumulator,
                                            const std::function< void (int, int, const double [], double [], double []) > &estimator,
                                            bool needsEquil)
    {
        _estims.emplace_back( [accu=accumulator.get(), estimator](double average[], double error[]) { // lambda functional
                                  if(!accu->isFinalized()) {
                                      throw std::runtime_error("[ObservableContainer.estim] Estimator was called, but accumulator is not finalized.");
                                  }
                                  estimator(accu->getNStore(), accu->getNObs(), accu->getData(), average, error);
                              } );
        _nobsdim+=accumulator->getNObs();
        _accus.emplace_back(std::move(accumulator)); // now accumulator is owned by _accus vector

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


    void ObservableContainer::printObsValues(std::ofstream &file) const
    {
        for (auto & accu : _accus) {
            for (int j=0; j<accu->getNObs(); ++j) {
                file << " " << accu->getObsValue(j);
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


    void ObservableContainer::estimate(double average[], double error[]) const
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
        //for (auto & accu : _accus) {
        //    delete accu;
        //}
        _accus.clear();
        _nobsdim=0;

        _estims.clear();
        _flags_equil.clear();
    }

}  // namespace mci
