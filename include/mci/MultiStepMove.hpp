#ifndef MCI_MULTISTEPMOVE_HPP
#define MCI_MULTISTEPMOVE_HPP

#include "mci/SRRDVecMove.hpp"
#include "mci/SamplingFunctionContainer.hpp"
#include "mci/TrialMoveInterface.hpp"

#include <numeric>

namespace mci
{
    // A TrialMove that contains own sampling functions (contained in pdfcont), together with an own trialMove,
    // to move through nsteps configurations before the end-result is returned to MCI and used to compute the
    // actual sampling function. Finally, the whole result is accepted or rejected.
    // This technique can be useful if the true PDF is very expensive to calculate, but a cheap approximate PDF
    // is available to move through configurations quickly while maintaining high acceptance ratio.
    //
    // Both the sub-PDFs and sub-trialMove can be set in a similar fashion to main MCI.
    // NOTE 1: If no PDF is added, the sub-sampling will accept every step (i.e. constant PDF)!
    // NOTE 2: For simplicity and speed, the sub-sampling itself does not consider any domain boundaries, but MCI
    // will apply them to the end result. Make sure that your sub-PDF does not rely on already applied PBC. This
    // should be a non-issue when only distances calculated via minimum-image convention enter the sub-PDF.
    class MultiStepMove: public TrialMoveInterface
    {
    protected:
        int _nsteps; // how many sub-sampling steps to do
        double * const _origX; // used to backup the original xold
        std::uniform_real_distribution<double> _rd; // used for own accept/reject
        std::unique_ptr<TrialMoveInterface> _trialMove; // the contained sub-move (init: uniform all-move)
        SamplingFunctionContainer _pdfcont; // sampling function container (init: empty)

        TrialMoveInterface * _clone() const final {
            auto * ret = new MultiStepMove(_ndim, _nsteps);
            ret->setTrialMove(*_trialMove);
            for (int i=0; i<_pdfcont.size(); ++i) {
                ret->addSamplingFunction( _pdfcont.getSamplingFunction(i) );
            }
            return ret;
        }

        void _newToOld(const WalkerState&/*wlk*/) final {}; // not needed
        void _oldToNew() final {}; // not needed

    public:
        MultiStepMove(int ndim, int nsteps):
            TrialMoveInterface(ndim, 0), _nsteps(nsteps), _origX(new double[ndim]),
            _rd( std::uniform_real_distribution<double>(0.,1.) ),
            _trialMove( new UniformVecMove(ndim, 1, 0.05) ) /*default to uniform single-move*/
        {}

        explicit MultiStepMove(int ndim): MultiStepMove(ndim, ndim) {} // default to ndim steps (with single index moves)
        ~MultiStepMove() final { delete [] _origX; }

        void setNSteps(int nsteps) { _nsteps=nsteps; }
        void setTrialMove(const TrialMoveInterface &tmove); // pass an existing move to be cloned
        void addSamplingFunction(const SamplingFunctionInterface &pdf); // add a sampling function (we make clone)
        void clearSamplingFunctions();

        int getNSteps() const { return _nsteps; }
        TrialMoveInterface & getTrialMove() const { return *_trialMove; }
        SamplingFunctionInterface & getSamplingFunction(int i) const { return _pdfcont.getSamplingFunction(i); }

        // Here we simply wrap the contained trialMove
        int getNStepSizes() const final { return _trialMove->getNStepSizes(); }
        double getStepSize(int i) const final { return _trialMove->getStepSize(i); }
        void setStepSize(int i, double val) final {_trialMove->setStepSize(i, val); }
        double getChangeRate() const final { return std::min(1., _trialMove->getChangeRate()*_nsteps); } // notice we multiply by nsteps here
        int getStepSizeIndex(int xidx) const final { return _trialMove->getStepSizeIndex(xidx); }

        // Methods used during sampling:
        void protoFunction(const double/*in*/[], double/*protov*/[]) final {} // not needed
        double trialMove(WalkerState &wlk, const double protoold[], double protonew[]) final;
    };

}; // namespace mci

#endif
