#ifndef MCI_DEPENDENTOBSERVABLEINTERFACE_HPP
#define MCI_DEPENDENTOBSERVABLEINTERFACE_HPP

#include "mci/AccumulatorInterface.hpp"
#include "mci/SamplingFunctionContainer.hpp"

#include <vector>

namespace mci
{
// Extension interface for MC observables depending on sampling functions and/or other observables
//
// An implementation must provide the same methods as in ObservableFunctionInterface and
// also the methods registerDeps(..) and deregisterDeps(), to set and null pointers to other objects.
// These methods will be called at the beginning and end of a sampling run, respectively.
//
// RULES:
// 1) If you depend on the sampling function(s), flag_pdfdep must be set (to let preparation take place)
// 2) You may only depend on accumulators that reside at an earlier position (i.e. lower index) than the
//    accumulator of "this" observable, because computation of observables happens in that order.
// 3) The nskip values of the respective accumulators of two dependent observables must lead to synced computation.
//

class DependentObservableInterface
{
protected:
    const bool _flag_pdfdep; // obs depends on sampling function (i.e. triggers their prepareObservation())

    explicit DependentObservableInterface(bool dependsOnPDF): _flag_pdfdep(dependsOnPDF) {}

public:
    bool dependsOnPDF() const { return _flag_pdfdep; }

    // The latter two rules above may be checked in registerDeps() by using the helper function below:
    static bool isObsDepValid(const std::vector<AccumulatorInterface *> &accuvec, int selfIdx, int depIdx) // when index thisIdx wants to depend on depIdx
    {
        const bool isOrdered = (depIdx < selfIdx); // depIdx observable must be computed before thisIdx
        const int thisNskip = accuvec[selfIdx]->getNSkip();
        const int depNskip = accuvec[depIdx]->getNSkip();
        const bool isSynced = (thisNskip >= depNskip) ? (thisNskip%depNskip == 0) : false; // nskips must be synced

        return (isOrdered && isSynced);
    }

    // --- MUST BE IMPLEMENTED
    // When initializing a sampling run, registerDeps(..) method will be called once.
    // In this method you may scan the provided containers for objects you depend on,
    // e.g. the sampling functions, and save pointers to use in the observableFunction methods.
    // These are guaranteed to be valid until deregisterDeps() is called by MCI (end of sampling).
    // In this method any pointers set on registerDeps() should be invalidated.
    virtual void registerDeps(const SamplingFunctionContainer &pdfcont, const std::vector<AccumulatorInterface *> &accuvec, int selfIdx/*index of this*/) = 0;
    virtual void deregisterDeps() = 0; // set all pointers to nullptr here
};
}  // namespace mci


#endif
