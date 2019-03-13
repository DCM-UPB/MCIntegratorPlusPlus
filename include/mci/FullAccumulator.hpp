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
        int _nstore; // number of allocated storage elements with _nobs length each
        int _storeidx; // storage index offset for next write

        // --- storage method to be implemented
        void _allocate() override;
        void _accumulate() override;
        void _finalize() override {} // nothing to do
        void _reset() override;
        void _deallocate() override;

    public:
        FullAccumulator(std::unique_ptr<ObservableFunctionInterface> obs, int nskip):
            AccumulatorInterface(std::move(obs), nskip), _nstore(0), _storeidx(0)
        {}

        ~FullAccumulator() override { this->_deallocate(); }

        int getNStore() const override { return _nstore; }
    };
}  // namespace mci

#endif
