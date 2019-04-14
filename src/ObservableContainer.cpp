#include "mci/ObservableContainer.hpp"

namespace mci
{

void ObservableContainer::addObservable(std::unique_ptr<ObservableFunctionInterface> obs,
                                        const int blocksize, const int nskip, const bool needsEquil, const EstimatorType estimType)
{
    ObservableContainerElement newElement;
    newElement.obs = std::move(obs); // ownership by element
    _nobsdim += newElement.obs->getNObs();
    newElement.accu = createAccumulator(*newElement.obs, blocksize, nskip); // use create from Factories.hpp

    // lambda functional (again use create from Factories.hpp)
    newElement.estim = [accu = newElement.accu.get() /*OK*/, estimator=createEstimator(estimType)](double average[], double error[])
    {
        if (!accu->isFinalized()) {
            throw std::runtime_error("[ObservableContainer.estim] Estimator was called, but accumulator is not finalized.");
        }
        estimator(accu->getNStore(), accu->getNObs(), accu->getData(), average, error);
    };

    newElement.flag_equil = needsEquil;
    _cont.push_back(std::move(newElement)); // and then into container
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

std::unique_ptr<ObservableFunctionInterface> ObservableContainer::pop_back()
{
    auto obs = std::move(_cont.back().obs); // move last obs out
    _nobsdim -= obs->getNObs(); // adjust nobsdim
    _cont.pop_back(); // resize vector
    return obs;
}

void ObservableContainer::clear()
{
    _cont.clear();
    _nobsdim = 0;
}
}  // namespace mci
