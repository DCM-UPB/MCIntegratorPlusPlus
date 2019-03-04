#ifndef MCI_MCIOBSERVABLECONTAINER_HPP
#define MCI_MCIOBSERVABLECONTAINER_HPP


#include "mci/MCIAccumulatorInterface.hpp"

struct MCIObservableContainer
{
    MCIAccumulatorInterface * const accu;
    const bool flag_error; // calculation of errors desired?

    //MCIEstimatorInterface * const estim; // not yet used

    MCIObservableContainer(MCIAccumulatorInterface * accumulator, bool needsError):
                           accu(accumulator), flag_error(needsError)
                           {}
};


#endif
