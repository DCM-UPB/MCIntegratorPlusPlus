#ifndef MCI_BLOCKACCUMULATOR_HPP
#define MCI_BLOCKACCUMULATOR_HPP

#include "mci/AccumulatorInterface.hpp"

#include <stdexcept>

namespace mci
{
    // Class to handle accumulation of observables, when averaging samples in blocks of fixed size
    // is desired. Typical use case is if you know how large the blocks have to be for uncorrelated samples
    // and want to avoid the memory&CPU overhead of using automatic blocking.
    // NOTE: The planned number of steps must be a multiple of the chosen blocksize.
    class BlockAccumulator: public AccumulatorInterface
    {
    protected:
        const int _blocksize; // how many samples to accumulate per block
        int _nblocks; // this will be set properly on allocation

        int _bidx; // counter to determine when block is finished
        int _storeidx; // storage index offset for next write

        // --- storage method to be implemented
        void _allocate() override;
        void _accumulate() override;
        void _finalize() override;
        void _reset() override;
        void _deallocate() override;

    public:
        BlockAccumulator(ObservableFunctionInterface * obs, int nskip, int blocksize):
            AccumulatorInterface(obs, nskip), _blocksize(blocksize), _nblocks(0), _bidx(0), _storeidx(0)
        {
            if (_blocksize < 1) { throw std::invalid_argument("[BlockAccumulator] Requested blocksize was < 1 ."); }
        }

        ~BlockAccumulator() override { this->_deallocate(); }

        int getBlockSize() { return _blocksize; }
        int getNStore() override { return _nblocks; }
    };
}

#endif
