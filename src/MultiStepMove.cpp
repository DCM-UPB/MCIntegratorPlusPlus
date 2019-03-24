#include "mci/MultiStepMove.hpp"

#include <algorithm>

namespace mci
{

    double MultiStepMove::trialMove(WalkerState &wlk, const double/*pold*/[], double/*pnew*/[]) // perform mini-MC
    {
        const bool orig_needsObs = wlk.needsObs; // store original need Obs flag
        std::copy(wlk.xold, wlk.xold+_ndim, _origX); // make backup of xold
        wlk.needsObs = false; // we don't do obs here

        _pdfcont.initializeProtoValues(wlk.xold); // initialize the sub-pdf at x
        const double oldPDF = _pdfcont.getOldSamplingFunction(); // remember this for later (faster)

        _trialMove->bindRGen(*_rgen); // we bind rgen here
        _trialMove->initializeProtoValues(wlk.xold); // initialize the sub-move

        for (int i=0; i<_nsteps; ++i) {
            // propose a new position x and get move acceptance
            const double moveAcc = _trialMove->computeTrialMove(wlk);
            // find the corresponding sampling function acceptance
            const double pdfAcc = _pdfcont.computeAcceptance(wlk);
            // determine if the proposed x is accepted or not
            wlk.accepted = ( _rd(*_rgen) <= pdfAcc*moveAcc );
            // set state according to result
            if (wlk.accepted) {
                _pdfcont.newToOld();
                _trialMove->newToOld();
                wlk.newToOld();
            } else { // rejected
                _pdfcont.oldToNew();
                _trialMove->oldToNew();
                wlk.oldToNew();
            }
        }
        const double newPDF = _pdfcont.getOldSamplingFunction(); // compute final PDF value

        // reset wlk to proper state
        std::copy(_origX, _origX+_ndim, wlk.xold); // set xold to original xold (xnew stays as is)
        wlk.nchanged = _ndim; // most indices should have changed, so let's treat this like an all-move
        wlk.needsObs = orig_needsObs;
        wlk.accepted = false; // reset this, to be sure

        // return move acceptance (inverse of "normal" acceptace), for detailed balance
        return oldPDF/newPDF;
    }


    // --- Trial Move Setter

    void MultiStepMove::setTrialMove(const TrialMoveInterface &tmove)
    {
        if (tmove.getNDim() != _ndim) {
            throw std::invalid_argument("[MultiStepMove::setTrialMove] Passed trial move's number of inputs is not equal to number of walkers.");
        }
        _trialMove = tmove.clone(); // unique ptr, old move gets freed automatically
        // we bind rgen later
    }


    // --- Sampling function

    void MultiStepMove::clearSamplingFunctions()
    {
        _pdfcont.clear();
    }

    void MultiStepMove::addSamplingFunction(const SamplingFunctionInterface &pdf)
    {
        if (pdf.getNDim() != _ndim) {
            throw std::invalid_argument("[MultiStepMove::addSamplingFunction] Passed sampling function's number of inputs is not equal to number of walkers.");
        }
        _pdfcont.addSamplingFunction( pdf.clone() );
    }

}  // namespace mci
