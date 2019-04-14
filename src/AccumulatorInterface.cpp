#include "mci/AccumulatorInterface.hpp"

namespace mci
{

AccumulatorInterface::AccumulatorInterface(ObservableFunctionInterface &obs, const int nskip):
        _obs(obs), _flag_updobs(_obs.isUpdateable()), _nobs(_obs.getNObs()), _xndim(_obs.getNDim()),
        _nskip(nskip), _obs_values(new double[_nobs]), _flags_xchanged(_flag_updobs ? new bool[_xndim] : nullptr),
        _nsteps(0), _data(nullptr)
{
    if (nskip < 1) { throw std::invalid_argument("[AccumulatorInterface] Provided number of steps per evaluation was < 1 ."); }
    this->_init();
}

AccumulatorInterface::~AccumulatorInterface()
{
    delete[] _flags_xchanged;
    delete[] _obs_values;
}

void AccumulatorInterface::_init() // reset base variables (except nsteps/_data)
{
    _stepidx = 0;
    _skipidx = _nskip - 1; // first step should not be skipped, so we prepare ++_skipidx == _nskip
    _flag_final = false;

    _nchanged = _xndim; // on the first step we always need to evaluate fully
    if (_flag_updobs) { std::fill(_flags_xchanged, _flags_xchanged + _xndim, true); }
}

void AccumulatorInterface::_processOld(const WalkerState &wlk, const SamplingFunctionContainer &pdfcont)
{
    // this is used when both !wlk.accepted and _nchanged==0
    if (++_skipidx == _nskip) { // accumulate observables
        _skipidx = 0;
        this->_accumulate(); // call child storage implementation
    }
}

void AccumulatorInterface::_processFull(const WalkerState &wlk, const SamplingFunctionContainer &pdfcont)
{
    // this is used when something changed (wlk.accepted || _nchanged>0) and obs is not updateable
    _nchanged = _xndim; // remember change even when we skip
    if (++_skipidx == _nskip) { // accumulate observables
        _skipidx = 0;

        // call full obs compute
        _obs.observableFunction(wlk.xnew, pdfcont, _obs_values);
        _nchanged = 0;

        this->_accumulate(); // call child storage implementation
    }
}

void AccumulatorInterface::_processSelective(const WalkerState &wlk, const SamplingFunctionContainer &pdfcont)
{   // this is used when something changed (wlk.accepted || _nchanged>0) and obs is updateable
    if (_nchanged < _xndim && wlk.accepted) { // we need to record changes
        if (wlk.nchanged < _xndim) { // track changes by index
            for (int i = 0; i < wlk.nchanged; ++i) {
                if (!_flags_xchanged[wlk.changedIdx[i]]) {
                    _flags_xchanged[wlk.changedIdx[i]] = true;
                    ++_nchanged; // increase internal change counter
                }
            }
        }
        else { // all-particle move case
            _nchanged = _xndim; // note: if _nchanged>=_xndim, the flags get ignored, so no need to set them
        }
    }

    if (++_skipidx == _nskip) { // accumulate observables
        _skipidx = 0;

        if (_nchanged < _xndim) { // call optimized recompute
            _obs.updatedObservable(wlk.xnew, _nchanged, _flags_xchanged, pdfcont, _obs_values);
        }
        else { // call full obs compute
            _obs.observableFunction(wlk.xnew, pdfcont, _obs_values);
        }
        std::fill(_flags_xchanged, _flags_xchanged + _xndim, false);
        _nchanged = 0;

        this->_accumulate(); // call child storage implementation
    }
}


void AccumulatorInterface::allocate(const int64_t nsteps)
{
    this->deallocate(); // for safety, also calls reset

    if (nsteps < 1) { throw std::invalid_argument("[AccumulatorInterface::allocate] Provided number of MC steps was < 1 ."); }

    _nsteps = nsteps;
    this->_allocate(); // call child allocate
}


void AccumulatorInterface::accumulate(const WalkerState &wlk, const SamplingFunctionContainer &pdfcont)
{
    if (wlk.accepted || _nchanged > 0) {
        if (_flag_updobs) {
            this->_processSelective(wlk, pdfcont);
        }
        else {
            this->_processFull(wlk, pdfcont);
        }
    }
    else {
        this->_processOld(wlk, pdfcont);
    }

    ++_stepidx;
}


void AccumulatorInterface::finalize()
{
    if (_stepidx != _nsteps) { throw std::runtime_error("[AccumulatorInterface::finalize] Finalize was called, but number of accumulated steps do not match the planned amount."); }
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
