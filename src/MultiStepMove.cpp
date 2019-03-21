#include "mci/MultiStepMove.hpp"

#include <algorithm>

namespace mci
{

    double MultiStepMove::trialMove(WalkerState &wlk, const double/*pold*/[], double/*pnew*/[]) // perform mini-MC
    {
        std::copy(wlk.xold, wlk.xold+_ndim, _origX); // make backup of xold
        _pdfcont.initializeProtoValues(wlk.xold); // initialize the sub-pdf at x
        const double oldPDF = _pdfcont.getOldSamplingFunction(); // remember this for later (faster)

        _trialMove->bindRGen(*_rgen); // we bind rgen late
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
                wlk.newToOld();
                _pdfcont.newToOld();
                _trialMove->newToOld();
            } else { // rejected
                wlk.oldToNew();
                _pdfcont.oldToNew();
                _trialMove->oldToNew();
            }
        }
        const double newPDF = _pdfcont.getOldSamplingFunction(); // compute final PDF value

        // reset wlk to proper state
        wlk.accepted = false; // reset this, to be sure
        std::copy(_origX, _origX+_ndim, wlk.xold); // set xold to original xold
        wlk.nchanged = _ndim; // most indices should have changed, so let's treat this like an all-move

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
