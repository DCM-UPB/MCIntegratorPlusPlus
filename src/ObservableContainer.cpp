#include "mci/ObservableContainer.hpp"

namespace mci
{

int gcd_helper(int a, int b) { // simple recursive greatest common divisor computation
    if (a > b) {
        return gcd_helper(a - b, b);
    }
    if (b > a) {
        return gcd_helper(a, b - a);
    }
    return a; // a == b
}

void ObservableContainer::_setDependsOnPDF()
{
    int gcd = 0;
    for (auto &el : _cont) {
        if (el.depobs != nullptr) {
            if (el.depobs->dependsOnPDF()) {
                if (gcd == 0) {
                    gcd = el.accu->getNSkip();
                }
                else if (gcd == 1) {
                    break; // no need to continue
                }
                else {
                    gcd = gcd_helper(gcd, el.accu->getNSkip());
                }
            }
        }
    }
    _nskip_PDF = gcd;
}

void ObservableContainer::addObservable(std::unique_ptr<ObservableFunctionInterface> obs,
                                        const int blocksize, const int nskip, const bool needsEquil, const EstimatorType estimType)
{
    ObservableContainerElement newElement;
    // obs+accu
    newElement.obs = std::move(obs); // ownership by element
    _nobsdim += newElement.obs->getNObs();
    newElement.depobs = dynamic_cast<DependentObservableInterface *>(newElement.obs.get()); // might be nullptr
    newElement.accu = createAccumulator(*newElement.obs, blocksize, nskip); // use create from Factories.hpp

    // estimator lambda functional (again use create from Factories.hpp)
    newElement.estim = [accu = newElement.accu.get() /*OK*/, estimator = createEstimator(estimType)](double average[], double error[])
    {
        if (!accu->isFinalized()) {
            throw std::runtime_error("[ObservableContainer.estim] Estimator was called, but accumulator is not finalized.");
        }
        estimator(accu->getNStore(), accu->getNObs(), accu->getData(), average, error);
    };

    newElement.flag_equil = needsEquil;
    _cont.push_back(std::move(newElement)); // and then into container
    this->_setDependsOnPDF(); // keep it simple and call this to update the depend flag
}


void ObservableContainer::allocate(const int64_t Nmc, const SamplingFunctionContainer &pdfcont)
{
    std::vector<AccumulatorInterface *> accuvec; // vectors of accu pointers for obs to register
    accuvec.reserve(_cont.size());
    for (auto &el : _cont) {
        el.accu->allocate(Nmc);
        accuvec.push_back(el.accu.get());
    }
    // let dependent obs register
    for (int i = 0; i < this->getNObs(); ++i) {
        if (_cont[i].depobs != nullptr) { _cont[i].depobs->registerDeps(pdfcont, accuvec, i); }
    }
}


void ObservableContainer::accumulate(const WalkerState &wlk)
{
    for (auto &el : _cont) {
        el.accu->accumulate(wlk);
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
    // let dependent obs deregister
    for (int i = 0; i < this->getNObs(); ++i) {
        if (_cont[i].depobs != nullptr) { _cont[i].depobs->deregisterDeps(); }
    }
}

std::unique_ptr<ObservableFunctionInterface> ObservableContainer::pop_back()
{
    auto obs = std::move(_cont.back().obs); // move last obs out
    _cont.pop_back(); // resize vector
    _nobsdim -= obs->getNObs(); // adjust nobsdim
    this->_setDependsOnPDF(); // adjust depend flag
    return obs;
}

void ObservableContainer::clear()
{
    _cont.clear();
    _nobsdim = 0;
    _nskip_PDF = 0;
}
}  // namespace mci
