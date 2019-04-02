#ifndef MCI_SIMPLEACCUMULATOR_HPP
#define MCI_SIMPLEACCUMULATOR_HPP

#include "mci/AccumulatorInterface.hpp"

namespace mci
{
// Class to handle accumulation of observables that don't require an error (typically gradients)
// This mean that all we have to do is to sum up the individual samples
class SimpleAccumulator final: public AccumulatorInterface
{
protected:
    bool _flag_alloc; // to determine proper nstored/ndata

    // --- storage method to be implemented
    void _allocate() final;
    void _accumulate() final;
    void _finalize() final;
    void _reset() final;
    void _deallocate() final;

public:
    SimpleAccumulator(std::unique_ptr<ObservableFunctionInterface> obs, int nskip):
            AccumulatorInterface(std::move(obs), nskip), _flag_alloc(false) {}

    ~SimpleAccumulator() final { this->_deallocate(); }

    int64_t getNStore() const final
    {
        return _flag_alloc ? 1 : 0; // we don't store old values
    }
};
}  // namespace mci

#endif
