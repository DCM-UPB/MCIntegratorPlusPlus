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

        void newToOld(); // copy new to old protovalues
        void computeOldProtoValues(const double xold[]); //compute the old proto sampling function with the old coordinates
        double computeAcceptance(const double xnew[]); //compute then new sampling function and return acceptance of new coordinates
        // compute new sampling from xold,xnew and _protoolds, with nchanged indices stored in changedIdx. return acceptance
        double computeAcceptance(const double xold[], const double xnew[], int nchanged, const int changedIdx[]);

        //void printProtoValues(std::ofstream &file) const; // write last protovalues to filestream
        void clear(); // clear everything
    };

} // namespace mci


#endif
