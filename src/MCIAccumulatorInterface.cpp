#include "mci/MCIAccumulatorInterface.hpp"

#include <stdexcept>

MCIAccumulatorInterface::MCIAccumulatorInterface(MCIObservableFunctionInterface * obs, const int nskip):
    _obs(obs), _obs_values(obs->getValues()), _nobs(obs->getNObs()), _nskip(nskip),
    _nsteps(0), _stepidx(0), _skipidx(0), _flag_eval(true), _flag_final(false), _data(nullptr)
{
    if (_nskip < 1) { throw std::invalid_argument("[MCIAccumulatorInterface] Provided number of steps per evaluation was < 1 ."); }
}


void MCIAccumulatorInterface::allocate(const int nsteps)
{
    this->deallocate(); // for safety

    if (nsteps < 1) { throw std::invalid_argument("[MCIAccumulatorInterface::allocate] Provided number of MC steps was < 1 ."); }

    _nsteps = nsteps;
    this->_allocate(); // call child allocate
}


void MCIAccumulatorInterface::accumulate(const double in[], const bool flagacc)
{
    if (_stepidx == _nsteps) { throw std::runtime_error("[MCIAccumulatorInterface::accumulate] Number of calls to accumulate exceed the allocation."); }

    if (flagacc) { _flag_eval = true; } // there was a change (so we need to evaluate obs on next skipidx==0)

    if (_skipidx == 0) { // accumulate observables
        if (_flag_eval) {
            _obs->computeValues(in);
            _flag_eval = false;
        }
        this->_accumulate(); // call child storage implementation
    }

    ++_stepidx;
    if (++_skipidx == _nskip) {
        _skipidx = 0; // next step will be evaluated
    }
}


void MCIAccumulatorInterface::finalize()
{
    if (_stepidx != _nsteps) { throw std::runtime_error("[MCIAccumulatorInterface::finalize] Finalize was called before all steps were accumulated."); }
    if (!_flag_final) { this->_finalize(); } // call child finalize
    _flag_final = true;
}


void MCIAccumulatorInterface::reset()
{
    this->_reset(); // call child reset

    // do base reset
    _stepidx = 0;
    _skipidx = 0;
    _flag_eval = true;
    _flag_final = false;
}


void MCIAccumulatorInterface::deallocate()
{
    this->reset(); // achieve clean state

    this->_deallocate(); // call child deallocate
    _nsteps = 0; // base set that one on allocate
}
