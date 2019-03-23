#include "mci/AccumulatorInterface.hpp"

#include <stdexcept>

namespace mci
{

    AccumulatorInterface::AccumulatorInterface(std::unique_ptr<ObservableFunctionInterface> obs, const int nskip):
        _obs(std::move(obs)), _updobs(dynamic_cast<UpdateableObservableInterface *>(_obs.get())), _flag_updobs(_updobs!=nullptr),
        _nobs(_obs->getNObs()), _xndim(_obs->getNDim()), _nskip(nskip),
        _obs_values(new double[_nobs]), _flags_xchanged(_flag_updobs ? new bool[_xndim] : nullptr),
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

    std::unique_ptr<ObservableFunctionInterface> AccumulatorInterface::removeObs()
    {
        return std::move(_obs); // move away the obs (NOW THE OBJECT IS INVALID; DELETE IT)
    }

    void AccumulatorInterface::_init() // reset base variables (except nsteps/_data)
    {
        _stepidx = 0;
        _skipidx = _nskip-1; // first step should not be skipped, so we prepare ++_skipidx == _nskip
        _flag_final = false;

        _nchanged = _xndim; // on the first step we always need to evaluate fully
        if (_flag_updobs) { std::fill(_flags_xchanged, _flags_xchanged+_xndim, true); }
    }


    void AccumulatorInterface::allocate(const int64_t nsteps)
    {
        this->deallocate(); // for safety, also calls reset

        if (nsteps < 1) { throw std::invalid_argument("[AccumulatorInterface::allocate] Provided number of MC steps was < 1 ."); }

        _nsteps = nsteps;
        this->_allocate(); // call child allocate
    }


    void AccumulatorInterface::accumulate(const WalkerState &wlk)
    {
        if (_stepidx >= _nsteps) { throw std::runtime_error("[AccumulatorInterface::accumulate] Number of calls to accumulate exceed the allocation."); }

        if (_nchanged<_xndim && wlk.accepted) { // we need to record changes
            if (_flag_updobs && wlk.nchanged<_xndim) { // track changes by index
                for (int i=0; i<wlk.nchanged; ++i) { // if nchange>0 (accepted step), we need to evaluate obs on next skipidx==0)
                    if (!_flags_xchanged[wlk.changedIdx[i]]) {
                        _flags_xchanged[wlk.changedIdx[i]] = true;
                        ++_nchanged; // increase internal change counter
                    }
                }
            }
            else { // all-particle move case or no tracking
                _nchanged = _xndim; // note: if _nchanged>=_xndim, the flags get ignored, so no need to set them
            }
        }

        if (++_skipidx == _nskip) { // accumulate observables
            _skipidx = 0;

            if (_nchanged > 0) { // we need to compute new obs
                if (_flag_updobs && _nchanged<_xndim) { // call optimized recompute
                    _updobs->updatedObservable(wlk.xnew, _nchanged, _flags_xchanged, _obs_values);
                } else { // call full obs compute
                    _obs->observableFunction(wlk.xnew, _obs_values);
                }
                if (_flag_updobs) { std::fill(_flags_xchanged, _flags_xchanged+_xndim, false); }
                _nchanged = 0;
            }

            this->_accumulate(); // call child storage implementation
        }
        ++_stepidx;
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
