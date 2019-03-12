#include "mci/ObservableContainer.hpp"

#include <stdexcept>

namespace mci
{

    void ObservableContainer::addObservable(std::unique_ptr<AccumulatorInterface> accumulator,
                                            const std::function< void (int, int, const double [], double [], double []) > &estimator,
                                            bool needsEquil)
    {
        _nobsdim+=accumulator->getNObs();
        ObservableContainerElement newElement;
        newElement.estim = [accu=accumulator.get() /*OK*/, estimator](double average[], double error[]) // lambda functional
                           {
                               if(!accu->isFinalized()) {
                                   throw std::runtime_error("[ObservableContainer.estim] Estimator was called, but accumulator is not finalized.");
                               }
                               estimator(accu->getNStore(), accu->getNObs(), accu->getData(), average, error);
                           };
        newElement.accu = std::move(accumulator); // ownership goes to new element
        newElement.flag_equil = needsEquil;
        _cont.emplace_back(std::move(newElement)); // and then into container
    }


    void ObservableContainer::allocate(const int Nmc)
    {
        for (auto & el : _cont) {
            el.accu->allocate(Nmc);
        }
    }


    void ObservableContainer::accumulate(const double x[], const int nchanged, const int changedIdx[])
    {
        for (auto & el : _cont) {
            el.accu->accumulate(x, nchanged, changedIdx);
        }
    }


    void ObservableContainer::printObsValues(std::ofstream &file) const
    {
        for (auto & el : _cont) {
            for (int j=0; j<el.accu->getNObs(); ++j) {
                file << " " << el.accu->getObsValue(j);
            }
        }
        file << " ";
    }


    void ObservableContainer::finalize()
    {
        for (auto & el : _cont) {
            el.accu->finalize();
        }
    }


    void ObservableContainer::estimate(double average[], double error[]) const
    {
        int offset = 0;
        for (auto & el : _cont) { // go through estimators and write to avg/error blocks
            el.estim(average+offset, error+offset);
            offset += el.accu->getNObs();
        }
    }


    void ObservableContainer::reset()
    {
        for (auto & el : _cont) {
            el.accu->reset();
        }
    }

    void ObservableContainer::deallocate()
    {
        for (auto & el : _cont) {
            el.accu->deallocate();
        }
    }

    void ObservableContainer::clear()
    {
        _cont.clear();
        _nobsdim=0;
    }

}  // namespace mci
