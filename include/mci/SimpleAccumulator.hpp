#ifndef MCI_SIMPLEACCUMULATOR_HPP
#define MCI_SIMPLEACCUMULATOR_HPP

#include "mci/AccumulatorInterface.hpp"

namespace mci
{
    // Class to handle accumulation of observables that don't require an error (typically gradients)
    // This mean that all we have to do is to sum up the individual samples
    class SimpleAccumulator: public AccumulatorInterface
    {
    protected:
        bool _flag_alloc; // to determine proper nstored/ndata

        // --- storage method to be implemented
        void _allocate() override;
        void _accumulate() override;
        void _finalize() override;
        void _reset() override;
        void _deallocate() override;

    public:
        SimpleAccumulator(const ObservableFunctionInterface &obs, int nskip):
            AccumulatorInterface(obs, nskip), _flag_alloc(false)
        {}

        ~SimpleAccumulator() override { this->_deallocate(); }

        int getNStore() const override {
            return _flag_alloc ? 1 : 0; // we don't store old values
        }
    };
}  // namespace mci

#endif
