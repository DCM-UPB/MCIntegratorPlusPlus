#include "mci/ObservableContainer.hpp"

namespace mci
{

void ObservableContainer::addObservable(std::unique_ptr<AccumulatorInterface> accumulator,
                                        const std::function<void(int64_t, int, const double [], double [], double [])> &estimator,
                                        bool needsEquil)
{
    _nobsdim += accumulator->getNObs();
    ObservableContainerElement newElement;
    newElement.estim = [accu = accumulator.get() /*OK*/, estimator](double average[], double error[]) // lambda functional
    {
        if (!accu->isFinalized()) {
            throw std::runtime_error("[ObservableContainer.estim] Estimator was called, but accumulator is not finalized.");
        }
        estimator(accu->getNStore(), accu->getNObs(), accu->getData(), average, error);
    };
    newElement.accu = std::move(accumulator); // ownership goes to new element
    newElement.flag_equil = needsEquil;
    _cont.emplace_back(std::move(newElement)); // and then into container
}


void ObservableContainer::allocate(const int64_t Nmc)
{
    for (auto &el : _cont) {
        el.accu->allocate(Nmc);
    }
}


void ObservableContainer::accumulate(const WalkerState &wlk, const SamplingFunctionContainer &pdfcont)
{
    for (auto &el : _cont) {
        el.accu->accumulate(wlk, pdfcont);
    }
}


void ObservableContainer::printObsValues(std::ofstream &file) const
{
    for (auto &el : _cont) {
        for (int j = 0; j < el.accu->getNObs(); ++j) {
            file << " " << el.accu->getObsValue(j);
        }
    }
    file << " ";
}


void ObservableContainer::finalize()
{
    for (auto &el : _cont) {
        el.accu->finalize();
    }
}


void ObservableContainer::estimate(double average[], double error[]) const
{
    int offset = 0;
    for (auto &el : _cont) { // go through estimators and write to avg/error blocks
        el.estim(average + offset, error + offset);
        offset += el.accu->getNObs();
    }
}


void ObservableContainer::reset()
{
    for (auto &el : _cont) {
        el.accu->reset();
    }
}

void ObservableContainer::deallocate()
{
    for (auto &el : _cont) {
        el.accu->deallocate();
    }
}

std::unique_ptr<AccumulatorInterface> ObservableContainer::pop_back()
{
    auto accu = std::move(_cont.back().accu); // move last accu out
    _nobsdim -= accu->getNObs(); // adjust nobsdim
    _cont.pop_back(); // resize vector
    return accu;
}

void ObservableContainer::clear()
{
    _cont.clear();
    _nobsdim = 0;
}
}  // namespace mci
