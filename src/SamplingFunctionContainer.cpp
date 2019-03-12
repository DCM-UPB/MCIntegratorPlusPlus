#include "mci/SamplingFunctionContainer.hpp"

namespace mci
{

    void SamplingFunctionContainer::addSamplingFunction(std::unique_ptr<SamplingFunctionInterface> sf /* we acquire ownership */)
    {
        _pdfs.emplace_back(std::move(sf)); // now sf is owned by _pdfs vector
    }

    void SamplingFunctionContainer::newToOld()
    {
        for (auto & sf : _pdfs) {
            sf->newToOld();
        }
    }

    void SamplingFunctionContainer::computeOldSamplingFunctions(const double xold[])
    {
        for (auto & sf : _pdfs) {
            sf->computeOldSamplingFunction(xold);
        }
    }

    double SamplingFunctionContainer::computeAcceptance(const double xnew[])
    {
        double acceptance = 1.;
        for (auto & sf : _pdfs) {
            acceptance *= sf->computeAcceptance(xnew);
        }
        return acceptance;
    }

    double SamplingFunctionContainer::computeAcceptance(const double xold[], const double xnew[], const int nchanged, const int changedIdx[])
    {
        double acceptance = 1.;
        for (auto & sf : _pdfs) {
            acceptance *= sf->computeAcceptance(xold, xnew, nchanged, changedIdx);
        }
        return acceptance;
    }


    void SamplingFunctionContainer::clear()
    {
        _pdfs.clear();
    }

}  // namespace mci
