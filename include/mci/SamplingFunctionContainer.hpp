#ifndef MCI_SAMPLINGFUNCTIONCONTAINER_HPP
#define MCI_SAMPLINGFUNCTIONCONTAINER_HPP

#include "mci/SamplingFunctionInterface.hpp"

#include <vector>
#include <memory>

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
        int getNSampF() const { return _pdfs.size(); }
        int size() const { return this->getNSampF(); }

        bool isEmpty() const { return ( _pdfs.size()==0 ); }
        bool empty() const { return this->isEmpty(); }

        const SamplingFunctionInterface & getSamplingFunction(int i) const { return *_pdfs[i]; }

        // operational methods

        void addSamplingFunction(std::unique_ptr<SamplingFunctionInterface> sampf); // we acquire ownership

        void updateSamplingFunctions(); // swap old and new protovalues
        void computeOldSamplingFunctions(const double xold[]); //compute the old sampling function with the old coordinates
        void computeNewSamplingFunctions(const double xnew[]); //compute the new sampling function with new coordinates
        // compute new sampling from xold,xnew and _protoolds, with information about change (nchanged indices in changedIdx)
        void computeNewSamplingFunctions(const double xold[], const double xnew[], int nchanged, const int changedIdx[]);
        double computeAcceptance() const; // compute the full acceptance probability

        //void printProtoValues(std::ofstream &file) const; // write last protovalues to filestream
        void clear(); // clear everything
    };


} // namespace mci


#endif
