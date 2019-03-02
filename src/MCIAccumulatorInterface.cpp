#include "mci/MCIAccumulatorInterface.hpp"

#include <stdexcept>

MCIAccumulatorInterface::MCIAccumulatorInterface(MCIObservableFunctionInterface * obs, const int nskip):
    _obs(obs), _nobs(obs->getNObs()), _nskip(nskip),
    _nsteps(0), _stepidx(0), _accuidx(0), _skipidx(0), _flag_eval(true), _data(nullptr)
{
    if (_nskip < 1) { throw std::invalid_argument("[MCIAccumulatorInterface] Provided number of steps per evaluation was < 1 ."); }
}


void MCIAccumulatorInterface::allocate(const int nsteps)
{
    if (nsteps < 1) { throw std::invalid_argument("[MCIAccumulatorInterface::allocate] Provided number of MC steps was < 1 ."); }

    this->deallocate(); // to be safe

    _nsteps = nsteps;
    this->_allocate(); // call child allocate
}


void MCIAccumulatorInterface::accumulate(const double * in, const bool flagacc)
{
    if (flagacc) { _flag_eval = true; } // there was a change

    if (_skipidx == 0) { // accumulate observables
        if (_flag_eval) {
            _obs->computeObservables(in);
            _flag_eval = false;
        }
        this->_accumulate(); // call child storage implementation

        ++_accuidx;
    }

    if (++_skipidx == _nskip) {
        _skipidx = 0; // next step will be evaluated
    }

    ++_stepidx;
}


void MCIAccumulatorInterface::reset()
{
    _stepidx = 0;
    _accuidx = 0;
    _skipidx = 0;
    _flag_eval = true;
    this->_reset(); // call child reset
}


void MCIAccumulatorInterface::deallocate()
{
    this->reset();
    this->_deallocate(); // call child deallocate
}
