#include "mci/SamplingFunctionContainer.hpp"

namespace mci
{

    void SamplingFunctionContainer::addSamplingFunction(std::unique_ptr<SamplingFunctionInterface> sf /* we acquire ownership */)
    {
        _pdfs.emplace_back(std::move(sf)); // now sf is owned by _pdfs vector
    }

    void SamplingFunctionContainer::updateSamplingFunctions()
    {
        for (auto & sf : _pdfs) {
            sf->newToOld();
        }
    }

    void SamplingFunctionContainer::computeOldSamplingFunctions(const double xold[])
    {
        for (auto & sf : _pdfs) {
            sf->computeNewSamplingFunction(xold);
            sf->newToOld();
        }
    }

    void SamplingFunctionContainer::computeNewSamplingFunctions(const double xnew[])
    {
        for (auto & sf : _pdfs) {
            sf->computeNewSamplingFunction(xnew);
        }
    }

    void SamplingFunctionContainer::computeNewSamplingFunctions(const double xold[], const double xnew[], const int nchanged, const int changedIdx[])
    {
        for (auto & sf : _pdfs) {
            sf->computeNewSamplingFunction(xold, xnew, nchanged, changedIdx);
        }
    }


    double SamplingFunctionContainer::computeAcceptance() const
    {
        double acceptance=1.;
        for (auto & sf : _pdfs) {
            acceptance*=sf->getAcceptance();
        }
        return acceptance;
    }

    void SamplingFunctionContainer::clear()
    {
        _pdfs.clear();
    }

}  // namespace mci
