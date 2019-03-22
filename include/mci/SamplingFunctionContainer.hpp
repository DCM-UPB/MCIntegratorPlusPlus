#ifndef MCI_SAMPLINGFUNCTIONCONTAINER_HPP
#define MCI_SAMPLINGFUNCTIONCONTAINER_HPP

#include "mci/SamplingFunctionInterface.hpp"
#include "mci/WalkerState.hpp"

#include <memory>
#include <vector>

namespace mci
{
    class SamplingFunctionContainer
    { // Container to store sampling functions for MCI
    private:
        // Sampling Functions
        std::vector< std::unique_ptr<SamplingFunctionInterface> > _pdfs;

    public:
        explicit SamplingFunctionContainer() = default;
        ~SamplingFunctionContainer() = default;

        // simple getters

        int size() const { return _pdfs.size(); }
        int getNPDF() const { return this->size(); }

        bool empty() const { return _pdfs.empty(); }
        bool hasPDF() const { return !this->empty(); }

        SamplingFunctionInterface & getSamplingFunction(int i) const { return *_pdfs[i]; }

        // operational methods

        void addSamplingFunction(std::unique_ptr<SamplingFunctionInterface> sf); // we acquire ownership

        void newToOld(const WalkerState &wlk); // copy new to old protovalues
        void oldToNew(); // copy old to new protovalues
        void initializeProtoValues(const WalkerState &wlk); // initialize the proto sampling values, given the WalkerState
        double getOldSamplingFunction() const; // returns the combined true sampling function value of the old step (potential use in trial moves)
        double computeAcceptance(const WalkerState &wlk); //compute then new sampling function and return acceptance of new coordinates

        //void printProtoValues(std::ofstream &file) const; // write last protovalues to filestream
        void pop_back() { _pdfs.pop_back(); } // remove last obs
        void clear() { _pdfs.clear(); } // clear everything
    };

} // namespace mci


#endif
