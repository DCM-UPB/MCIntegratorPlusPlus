#ifndef MCI_UNIFORMNTYPEVECMOVE_HPP
#define MCI_UNIFORMNTYPEVECMOVE_HPP

#include "mci/TrialMoveInterface.hpp"

#include <stdexcept>
#include <algorithm>

namespace mci
{
    // Uniform single-vector move that can be used when you have different
    // types of particles/vectors, and want to calibrate and use different
    // step size per particle type. Vectors of a common type are expected
    // to be found in a single block within the position arrays. Pass the
    // typeEnds array containing the end-indices of every type block. We
    // follow the convention that "end" means one behind the last element.
    class UniformNTypeVecMove: public TrialMoveInterface
    {
    private:
        const int _nvecs; // how many vector/particles are considered (calculated as ndim/veclen)
        const int _veclen; // how many indices does one particle/vector have? (i.e. space dimension)
        const int _ntypes; // how many different types of particles do you have?
        int * const _typeEnds; // end-indices of every type in x (i.e. the last index of type i is _typeEnds[i]-1 )
        double * const _stepSizes; // holds the step sizes, one per type
        std::uniform_int_distribution<int> _rdidx; // uniform integer distribution to choose vector index
        std::uniform_real_distribution<double> _rdmov; // uniform double distribution to move vector

    protected:
        TrialMoveInterface * _clone() const override {
            return new UniformNTypeVecMove(_nvecs, _veclen, _ntypes, _typeEnds, _stepSizes);
        }

    public:
        UniformNTypeVecMove(int nvecs, int veclen, int ntypes, const int typeEnds[] /*len ntypes*/, double initStepSize /*scalar init*/):
            TrialMoveInterface(nvecs*veclen, 0), _nvecs(nvecs), _veclen(veclen), _ntypes(ntypes),
            _typeEnds(new int[_ntypes]), _stepSizes(new double[_ntypes]),
            _rdidx(std::uniform_int_distribution<int>(0, _nvecs-1)),
            _rdmov(std::uniform_real_distribution<double>(-1.,1.))
        {
            if (_nvecs < 1) { throw std::invalid_argument("[UniformNTypeVecMove] Number of vectors must be at least 1."); }
            if (_veclen < 1) { throw std::invalid_argument("[UniformNTypeVecMove] Vector length must be at least 1."); }
            if (_ntypes < 1) { throw std::invalid_argument("[UniformNTypeVecMove] Number of types must be at least 1."); }
            for (int i=0; i<_ntypes; ++i) { // we rely on this later
                if (typeEnds[i]%_veclen != 0) {
                    throw std::invalid_argument("[UniformNTypeVecMove] All type end indices must be multiples of vector length.");
                }
            }
            std::copy(typeEnds, typeEnds+_ntypes, _typeEnds);
            std::fill(_stepSizes, _stepSizes+_ntypes, initStepSize);
        }

        UniformNTypeVecMove(int nvecs, int veclen, int ntypes, const int typeEnds[], const double initStepSizes[] /*len ntypes*/):
            UniformNTypeVecMove(nvecs, veclen, ntypes, typeEnds, 0.) // reuse above constructor
        {
            std::copy(initStepSizes, initStepSizes+_ntypes, _stepSizes); // put the proper values in
        }

        ~UniformNTypeVecMove() { delete [] _stepSizes; }

        // Methods required for auto-callibration:
        int getNStepSizes() const override { return _ntypes; }
        void setStepSize(int i, double val) override { _stepSizes[i] = val; }
        double getStepSize(int i) const override { return _stepSizes[i]; }
        double getChangeRate() const override { return 1./_nvecs; } // equivalent to _veclen/_ndim
        void getUsedStepSizes(int /*nchangedX*/, const int changedIdx[], int &nusedSizes, int usedSizeIdx[]) const override
        { // we know that we changed only a single vector
            nusedSizes = 1;
            for (int i=0; i<_ntypes; ++i) {
                if (changedIdx[0] < _typeEnds[i]) {
                    usedSizeIdx[0] = i;
                    return;
                }
            }
            // this method is not performance-critical, so let's check for this
            throw std::runtime_error("[UniformNTypeVecMove::getUsedStepSizes] End of method reached, without result.");
        }

        void protoFunction(const double[], double[]) override {} // not needed

        void onAcceptance(const SamplingFunctionContainer &, double[]) override {} // not needed

        double trialMove(double xnew[], int &nchanged, int changedIdx[], const double[], double[]) override
        {
            // determine vector to change and its type
            const int vidx = _rdidx( *(this->getRGen()) );
            const int xidx = vidx*_veclen; // first x index to change
            int tidx = 0; // type index
            while (tidx<_ntypes) { // This executes fast, don't worry ;-)
                if (xidx < _typeEnds[tidx]) {
                    break;
                }
                ++tidx;
            }

            // do step
            for (int i=0; i<_veclen; ++i) {
                xnew[xidx + i] += _stepSizes[tidx] * _rdmov( *(this->getRGen()) );
                changedIdx[i] = xidx + i;
            }
            nchanged = _veclen; // how many indices we changed

            return 1.; // uniform -> no move acceptance factor
        }
    };

} // namespace mci

#endif
