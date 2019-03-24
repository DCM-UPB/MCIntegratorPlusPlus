#ifndef MCI_FULLACCUMULATOR_HPP
#define MCI_FULLACCUMULATOR_HPP

#include "mci/AccumulatorInterface.hpp"

namespace mci
{
// Class to handle accumulation of observables, when storing every single sample is desired
// Typically you want this for automatic blocking techniques after the sampling run
//
// NOTE: If the dimension of the observable is very large, consider using BlockAccumulator
// or SimpleAccumulator (if no error is required) instead, because the memory requirements
// of the FullAccumulator may become very large with a large number of MC steps.
//
class FullAccumulator: public AccumulatorInterface
{
protected:
    int64_t _nstore; // number of allocated storage elements with _nobs length each
    int64_t _storeidx; // storage index offset for next write

    // --- storage method to be implemented
    void _allocate() final;
    void _accumulate() final;
    void _finalize() final {} // nothing to do
    void _reset() final;
    void _deallocate() final;

public:
    FullAccumulator(std::unique_ptr<ObservableFunctionInterface> obs, int nskip):
            AccumulatorInterface(std::move(obs), nskip), _nstore(0), _storeidx(0) {}

    ~FullAccumulator() final { this->_deallocate(); }

    int64_t getNStore() const final { return _nstore; }
};
}  // namespace mci

#endif
