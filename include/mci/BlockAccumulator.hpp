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
class BlockAccumulator final: public AccumulatorInterface
{
protected:
    const int _blocksize; // how many samples to accumulate per block
    int64_t _nblocks; // this will be set properly on allocation

    int _bidx; // counter to determine when block is finished
    int64_t _storeidx; // storage index offset for next write

    // --- storage method to be implemented
    void _allocate() final;
    void _accumulate() final;
    void _finalize() final;
    void _reset() final;
    void _deallocate() final;

public:
    BlockAccumulator(ObservableFunctionInterface &obs, int nskip, int blocksize):
            AccumulatorInterface(obs, nskip), _blocksize(blocksize), _nblocks(0), _bidx(0), _storeidx(0)
    {
        if (_blocksize < 1) { throw std::invalid_argument("[BlockAccumulator] Requested blocksize was < 1 ."); }
    }

    ~BlockAccumulator() final { this->_deallocate(); }

    int getBlockSize() const { return _blocksize; }
    int64_t getNStore() const final { return _nblocks; }
};
}  // namespace mci

#endif
