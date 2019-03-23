#include "mci/SamplingFunctionContainer.hpp"

namespace mci
{

    void SamplingFunctionContainer::addSamplingFunction(std::unique_ptr<SamplingFunctionInterface> sf /* we acquire ownership */)
    {
        _pdfs.emplace_back(std::move(sf)); // now sf is owned by _pdfs vector
    }

    void SamplingFunctionContainer::newToOld(const WalkerState &wlk)
    {
        for (auto & sf : _pdfs) {
            sf->newToOld(wlk);
        }
    }

    void SamplingFunctionContainer::oldToNew()
    {
        for (auto & sf : _pdfs) {
            sf->oldToNew();
        }
    }

    void SamplingFunctionContainer::initializeProtoValues(const WalkerState &wlk)
    {
        for (auto & sf : _pdfs) {
            sf->initializeProtoValues(wlk);
        }
    }

    double SamplingFunctionContainer::getOldSamplingFunction() const
    {
        double sampf = 1.;
        for (auto & sf : _pdfs) {
            sampf *= sf->getOldSamplingFunction();
        }
        return sampf;
    }

    double SamplingFunctionContainer::computeAcceptance(const WalkerState &wlk)
    {
        double acceptance = 1.;
        for (auto & sf : _pdfs) {
            acceptance *= sf->computeAcceptance(wlk);
        }
        return acceptance;
    }

    std::unique_ptr<SamplingFunctionInterface> SamplingFunctionContainer::pop_back()
    {
        auto pdf = std::move(_pdfs.back()); // move last pdf out of vector
        _pdfs.pop_back();
        return pdf;
    }

}  // namespace mci
