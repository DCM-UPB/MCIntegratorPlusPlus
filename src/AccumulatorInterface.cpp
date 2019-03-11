#include "mci/AccumulatorInterface.hpp"

#include <stdexcept>

namespace mci
{

    AccumulatorInterface::AccumulatorInterface(const ObservableFunctionInterface &obs, const int nskip):
        _obs(obs.clone()), _nobs(_obs->getNObs()), _xndim(_obs->getNDim()), _nskip(nskip),
        _obs_values(new double[_nobs]), _flags_xchanged(new bool[_xndim]),
        _nsteps(0), _data(nullptr)
    {
        if (nskip < 1) { throw std::invalid_argument("[AccumulatorInterface] Provided number of steps per evaluation was < 1 ."); }
        this->_init();
    }

    AccumulatorInterface::~AccumulatorInterface()
    {
        delete [] _flags_xchanged;
        delete [] _obs_values;
    }

    void AccumulatorInterface::_init() // reset base variables (except nsteps/_data)
    {
        _stepidx = 0;
        _skipidx = 0;
        _flag_final = false;

        _nchanged = _xndim; // on the first step we always need to evaluate fully
        std::fill(_flags_xchanged, _flags_xchanged+_xndim, true);
    }


    void AccumulatorInterface::allocate(const int nsteps)
    {
        this->deallocate(); // for safety, also calls reset

        if (nsteps < 1) { throw std::invalid_argument("[AccumulatorInterface::allocate] Provided number of MC steps was < 1 ."); }

        _nsteps = nsteps;
        this->_allocate(); // call child allocate
    }


    void AccumulatorInterface::accumulate(const double x[], const int nchanged, const int changedIdx[])
    {
        if (_stepidx == _nsteps) { throw std::runtime_error("[AccumulatorInterface::accumulate] Number of calls to accumulate exceed the allocation."); }

        if (_nchanged<_xndim) { // if internal change counter is already full, we should skip this
            if (nchanged<_xndim) { // similarly for the step's nchange
                for (int i=0; i<nchanged; ++i) { // if nchange>0 (accepted step), we need to evaluate obs on next skipidx==0)
                    if (!_flags_xchanged[changedIdx[i]]) {
                        _flags_xchanged[changedIdx[i]] = true;
                        ++_nchanged; // increase internal change counter
                    }
                }
            }
            else { // all-particle move case
                _nchanged = _xndim; // when _nchanged>=_xndim, the flags get ignored so we don't need to set them
            }
        }

        if (_skipidx == 0) { // accumulate observables
            if (_nchanged>=_xndim) { // call full obs compute
                _obs->observableFunction(x, _obs_values);
            } else { // call optimized recompute
                _obs->observableFunction(x, _nchanged, _flags_xchanged, _obs_values);
            }
            this->_accumulate(); // call child storage implementation

            std::fill(_flags_xchanged, _flags_xchanged+_xndim, false);
            _nchanged = 0;
        }

        ++_stepidx;
        if (++_skipidx == _nskip) {
            _skipidx = 0; // next step will be evaluated
        }
    }


    void AccumulatorInterface::finalize()
    {
        if (_stepidx != _nsteps) { throw std::runtime_error("[AccumulatorInterface::finalize] Finalize was called before all steps were accumulated."); }
        if (!_flag_final) { this->_finalize(); } // call child finalize
        _flag_final = true;
    }


    void AccumulatorInterface::reset()
    {
        this->_reset(); // call child reset
        this->_init(); // call base init/reset
    }


    void AccumulatorInterface::deallocate()
    {
        this->reset(); // achieve clean state
        this->_deallocate(); // call child deallocate

        _nsteps = 0;
    }

}  // namespace mci
