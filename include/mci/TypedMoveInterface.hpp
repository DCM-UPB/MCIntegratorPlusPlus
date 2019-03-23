#ifndef MCI_TYPEDMOVEINTERFACE_HPP
#define MCI_TYPEDMOVEINTERFACE_HPP

#include "mci/TrialMoveInterface.hpp"

#include <algorithm>
#include <stdexcept>

namespace mci
{
    // Specialized interface for moves supporting different stepSizes for specified blocks
    // within the position array. Use this with ntypes>1 if you have e.g. different types
    // of particles and want to calibrate and use different step size per particle typeypes>1.
    // You need to pass the typeEnds array containing the end-indices of every type block. We
    // follow the convention that "end" means one behind the last element.
    //
    // Example for types a,b,c :
    // x = (a0, a1, a2, b0, b1, c0) ->  ntypes=3, typeEnds = (4,6,7)
    //
    class TypedMoveInterface: public TrialMoveInterface
    {
    protected:
        const int _ntypes; // how many different types of particles do you have?
        int * const _typeEnds; // end-indices of every type in x (i.e. the last index of type i is _typeEnds[i]-1 )
        double * const _stepSizes; // holds the step sizes, one per type

    public:
        TypedMoveInterface(int ndim, int nproto, int ntypes, const int typeEnds[] /*len ntypes*/, double initStepSize /*scalar init*/):
            TrialMoveInterface(ndim, nproto), _ntypes(ntypes),
            _typeEnds(new int[_ntypes]), _stepSizes(new double[_ntypes])
        {
            if (_ntypes < 1) { throw std::invalid_argument("[TypedMoveInterface] Number of types must be at least 1."); }
            if (_ntypes > 1) {
                if (typeEnds==nullptr) { throw std::invalid_argument("[TypedMoveInterface] When ntypes>1, passed typeEnds must not be null."); }
                std::copy(typeEnds, typeEnds+_ntypes, _typeEnds);
            } else {
                _typeEnds[0] = _ndim;
            }
            std::fill(_stepSizes, _stepSizes+_ntypes, initStepSize);
        }

        ~TypedMoveInterface() override {
            delete [] _stepSizes;
            delete [] _typeEnds;
        }

        // Methods that we can implement final:
        int getNStepSizes() const final { return _ntypes; }
        void setStepSize(int i, double val) final { _stepSizes[i] = val; }
        double getStepSize(int i) const final { return _stepSizes[i]; }

        int getStepSizeIndex(int xidx) const final {
            for (int i=0; i<_ntypes; ++i) {
                if (xidx < _typeEnds[i]) {
                    return i;
                }
            }
            throw std::runtime_error("[TypedMoveInterface::getUsedStepSize] Passed xidx exceeds expected range.");
        }
    };

} // namespace mci


#endif
