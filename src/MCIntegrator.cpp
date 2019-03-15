#include "mci/MCIntegrator.hpp"

#include "mci/BlockAccumulator.hpp"
#include "mci/FullAccumulator.hpp"
#include "mci/SimpleAccumulator.hpp"
#include "mci/Estimators.hpp"
#include "mci/Factories.hpp"
#include "mci/UniformAllMove.hpp"
#include "mci/UniformVecMove.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace mci
{

    //  --- Integrate
    void MCI::integrate(const int Nmc, double average[], double error[], const bool doFindMRT2step, const bool doDecorrelation)
    {
        if ( !_pdfcont.empty() ) {
            //find the optimal mrt2 step
            if (doFindMRT2step) { this->findMRT2Step(); }
            // take care to do the initial decorrelation of the walker
            if (doDecorrelation) { this->initialDecorrelation(); }
        }

        // allocation of the accumulators where the data will be stored
        _obscont.allocate(Nmc);

        //sample the observables
        if (_flagobsfile) { _obsfile.open(_pathobsfile); }
        if (_flagwlkfile) { _wlkfile.open(_pathwlkfile); }
        _flagMC = true;
        this->sample(Nmc, _obscont); // let sample accumulate data
        _flagMC = false;
        if (_flagobsfile) { _obsfile.close(); }
        if (_flagwlkfile) { _wlkfile.close(); }

        // estimate average and standard deviation
        _obscont.estimate(average, error);

        // if we sampled uniformly, scale results by volume
        if (_pdfcont.empty()) {
            for (int i=0; i<_obscont.getNObsDim(); ++i) {
                average[i] *=_vol;
                error[i] *=_vol;
            }
        }

        // deallocate
        _obscont.deallocate();
    }



    // --- "High-level" internal methods

    void MCI::initialDecorrelation()
    {
        if (_NdecorrelationSteps < 0) {
            // automatic equilibration of contained observables with flag_equil = true

            //create the temporary observable container to be used
            ObservableContainer obs_equil;
            for (int i=0; i<_obscont.getNObs(); ++i) {
                if (_obscont.getFlagEquil(i)) {
                    obs_equil.addObservable(std::unique_ptr<AccumulatorInterface>( new FullAccumulator(_obscont.getObservableFunction(i).clone(), 1) ),
                                            mci::CorrelatedEstimator, true);
                }
            }
            const int MIN_NMC=100;
            const int nobsdim = obs_equil.getNObsDim();
            // allocate memory for observables
            obs_equil.allocate(MIN_NMC);

            //do a first estimate of the observables
            this->sample(MIN_NMC, obs_equil);
            auto * oldestimate = new double[nobsdim];
            auto * olderrestim = new double[nobsdim];
            obs_equil.estimate(oldestimate, olderrestim);

            //start a loop which will stop when the observables are stabilized
            bool flag_loop=true;
            auto * newestimate = new double[nobsdim];
            auto * newerrestim = new double[nobsdim];
            while ( flag_loop ) {
                flag_loop = false;
                this->sample(MIN_NMC, obs_equil);
                obs_equil.estimate(newestimate, newerrestim);

                for (int i=0; i<nobsdim; ++i) {
                    if ( fabs( oldestimate[i] - newestimate[i] ) > 2*sqrt( olderrestim[i]*olderrestim[i] + newerrestim[i]*newerrestim[i] ) ) {
                        flag_loop=true; // if any difference is too large, continue the loop
                        break; // break the inner for loop
                    }
                }
                // swap array pointers
                std::swap(oldestimate, newestimate);
                std::swap(olderrestim, newerrestim);
            }

            //memory deallocation
            delete [] newerrestim;
            delete [] newestimate;
            delete [] olderrestim;
            delete [] oldestimate;
            obs_equil.clear();
        }
        else if (_NdecorrelationSteps > 0) {
            this->sample(_NdecorrelationSteps);
        }
    }


    void MCI::findMRT2Step()
    {
        if (!_trialMove->hasStepSizes()) { return; } // in the odd case that our mover has no adjustable step sizes

        //constants
        const int MIN_STAT=200;  //minimum statistic: number of M(RT)^2 steps done before deciding if the step must be increased or decreased
        const int MIN_CONS=5;   //minimum consecutive: minimum number of consecutive loops without need of changing mrt2step
        const double TOLERANCE=0.05;  //tolerance: tolerance for the acceptance rate
        const int MAX_NUM_ATTEMPTS=50;  //maximum number of attempts: maximum number of time that the main loop can be executed
        const double SMALLEST_ACCEPTABLE_DOUBLE=1.e-50;

        //initialize index
        int cons_count = 0;  //number of consecutive loops without need of changing mrt2step
        int counter = 0;  //counter of loops
        while ( ( _NfindMRT2Iterations < 0 && cons_count < MIN_CONS ) || counter < _NfindMRT2Iterations ) {
            const int nsizes = _trialMove->getNStepSizes();
            const double minboxlen = this->getMinBoxLen();

            //do MIN_STAT M(RT)^2 steps
            this->sample(MIN_STAT);

            //increase or decrease mrt2step depending on the acceptance rate
            const double rate = this->getAcceptanceRate();
            if ( fabs(rate-_targetaccrate) < TOLERANCE ) {
                //mrt2step was ok
                cons_count++;
            }
            else {
                //need to change mrt2step
                cons_count=0;
                const double fact = std::min(2., std::max(0.5, rate/_targetaccrate) );
                for (int j=0; j<nsizes; ++j) {
                    _trialMove->scaleStepSize(j, fact);
                }

                // sanity checks
                for (int j=0; j<nsizes; ++j) { //mrt2step = Infinity
                    if ( _trialMove->getStepSize(j) > minboxlen ) {
                        _trialMove->setStepSize(j, minboxlen);
                    }
                }
                for (int j=0; j<nsizes; ++j) { //mrt2step ~ 0
                    if ( _trialMove->getStepSize(j) < SMALLEST_ACCEPTABLE_DOUBLE ) {
                        _trialMove->setStepSize(j, SMALLEST_ACCEPTABLE_DOUBLE);
                    }
                }
            }
            counter++;

            if ( _NfindMRT2Iterations < 0 && counter >= MAX_NUM_ATTEMPTS ) {
                std::cout << "Warning [MCI::findMRT2Step]: Max number of attempts reached without convergence." << std::endl;
                break;
            }
        }
    }

    // --- Sampling

    void MCI::sample(const int npoints) // sample without taking observables, file io or callbacks
    {
        // Initialize
        this->initializeSampling(nullptr);

        //run the main loop for sampling
        const bool flagpdf = !(_pdfcont.empty());
        int changedIdx[_ndim]; // index array for tracking changed indices
        std::iota(changedIdx, changedIdx+_ndim, 0); // init with 0...ndim-1 (just because we can)
        for (_ridx=0; _ridx<npoints; ++_ridx) {
            if (flagpdf) { // use sampling function
                this->doStepMRT2(changedIdx);
            }
            else {
                this->newRandomX();
                ++_acc; // "accept" move
            }
        }
    }

    void MCI::sample(const int npoints, ObservableContainer &container)
    {
        // Initialize
        this->initializeSampling(&container);

        //run the main loop for sampling
        const bool flagpdf = !(_pdfcont.empty());
        int changedIdx[_ndim]; // index array for tracking changed indices
        std::iota(changedIdx, changedIdx+_ndim, 0); // init with 0...ndim-1 (just because we can)
        for (_ridx=0; _ridx<npoints; ++_ridx) {
            // do MC step
            int nchanged;
            if (flagpdf) { // use sampling function
                nchanged = this->doStepMRT2(changedIdx);
            }
            else {
                this->newRandomX();
                ++_acc; // "accept" move
                nchanged = _ndim;
            }
            // call the callbacks
            this->callBackOnMove(_xold, nchanged>0);

            // accumulate observables
            container.accumulate(_xold, nchanged, changedIdx);

            // file output
            if (_flagMC && _flagobsfile) { this->storeObservables(); } // store obs on file
            if (_flagMC && _flagwlkfile) { this->storeWalkerPositions(); } // store walkers on file
        }

        // finalize data
        container.finalize();
    }


    // --- Walking

    int MCI::doStepMRT2(int changedIdx[])
    {
        // propose a new position x
        int nchanged;
        const double moveAcc = _trialMove->computeTrialMove(_xnew, nchanged, changedIdx);
        applyPBC(_xnew);

            /*this->computeNewX();
        const double moveAcc = 1.;
        nchanged = _ndim;*/

        // find the corresponding sampling function value
        const double pdfAcc = (nchanged < _ndim) ?
                                           _pdfcont.computeAcceptance(_xold, _xnew, nchanged, changedIdx)
                                           : _pdfcont.computeAcceptance(_xnew);

        //determine if the proposed x is accepted or not
        nchanged = ( _rd(_rgen) <= pdfAcc * moveAcc /* maybe it must be / */ ) ? nchanged : 0;

        //update some values according to the acceptance of the mrt2 step
        if ( nchanged > 0 ) {
            //accepted
            this->acceptX();
        } else {
            //rejected
            this->rejectX();
        }

        return nchanged;
    }

    void MCI::acceptX()
    {
        std::copy(_xnew, _xnew+_ndim, _xold);
        _pdfcont.newToOld();
        _trialMove->newToOld();
        _trialMove->callOnAcceptance(_pdfcont);
        ++_acc;
    }

    void MCI::rejectX()
    {
        std::copy(_xold, _xold+_ndim, _xnew);
        _pdfcont.oldToNew();
        _trialMove->oldToNew();
        ++_rej;
    }

    void MCI::newRandomX()
    {
        //set xold to new random values (within the irange)
        for (int i=0; i<_ndim; ++i) {
            _xold[i] = _lbound[i] + ( _ubound[i] - _lbound[i] ) * _rd(_rgen);
        }
    }

    void MCI::initializeSampling(ObservableContainer * obsCont)
    {
        // reset running counters
        _acc = 0;
        _rej = 0;
        _ridx = 0;

        // init xnew and all protovalues
        std::copy(_xold, _xold+_ndim, _xnew);
        _pdfcont.initializeProtoValues(_xold); // initialize the pdf at x
        _trialMove->initializeProtoValues(_xold); // initialize the trial mover

        // call callbacks and initialize passed obs container
        this->callBackOnMove(_xold, true); // first call of the call-back functions
        if (obsCont != nullptr) { obsCont->reset(); }// reset observable data (passed per argument!)
    }

    /* // deprecated
    void MCI::computeNewX()
    {
        for (int i=0; i<_ndim; ++i) {
            _xnew[i] = _xold[i] + _mrt2step[i] * (_rd(_rgen)-0.5);
        }
        applyPBC(_xnew);
    }
    */


    // --- File Output

    void MCI::storeObservablesOnFile(const std::string &filepath, const int freq)
    {
        _pathobsfile = filepath;
        _freqobsfile = freq;
        _flagobsfile = true;
    }

    void MCI::storeObservables()
    {
        if ( _ridx%_freqobsfile == 0 ) {
            _obscont.printObsValues(_obsfile);
            _obsfile << std::endl;
        }
    }


    void MCI::storeWalkerPositionsOnFile(const std::string &filepath, const int freq)
    {
        _pathwlkfile = filepath;
        _freqwlkfile = freq;
        _flagwlkfile = true;
    }

    void MCI::storeWalkerPositions()
    {
        if ( _ridx%_freqwlkfile == 0 ) {
            _wlkfile << _ridx;
            for (int j=0; j<_ndim; ++j) {
                _wlkfile << "   " << _xold[j] ;
            }
            _wlkfile << std::endl;
        }
    }


    // --- Observables

    void MCI::clearObservables()
    {
        _obscont.clear();
    }

    void MCI::addObservable(const ObservableFunctionInterface &obs, int blocksize, int nskip, const bool flag_equil, const bool flag_correlated)
    {
        // sanity
        if (obs.getNDim() != _ndim) {
            throw std::invalid_argument("[MCI::addObservable] Passed observable function's number of inputs is not equal to MCI's number of walkers.");
        }
        blocksize = std::max(0, blocksize);
        nskip = std::max(1, nskip);
        const bool flag_error = (blocksize > 0); // will we calculate errors?
        if (flag_equil && blocksize==0) {
            throw std::invalid_argument("[MCI::addObservable] Requested automatic observable equilibration requires blocksize > 0.");
        }

        // add accumulator&estimator from factory functions
        _obscont.addObservable(createAccumulator(obs, blocksize, nskip), createEstimator(flag_correlated, flag_error), flag_equil);
    }


    // --- Sampling functions

    void MCI::clearSamplingFunctions()
    {
        _pdfcont.clear();
    }

    void MCI::addSamplingFunction(const SamplingFunctionInterface &mcisf)
    {
        if (mcisf.getNDim() != _ndim) {
            throw std::invalid_argument("[MCI::addObservable] Passed sampling function's number of inputs is not equal to MCI's number of walkers.");
        }
        _pdfcont.addSamplingFunction( mcisf.clone() );
    }


    // --- Callbacks

    void MCI::clearCallBackOnMove() {
        _cbacks.clear();
    }

    void MCI::addCallBackOnMove(const CallBackOnMoveInterface &cback)
    {
        if (cback.getNDim() != _ndim) {
            throw std::invalid_argument("[MCI::addObservable] Passed callback function's number of inputs is not equal to MCI's number of walkers.");
        }
        _cbacks.emplace_back( std::unique_ptr<CallBackOnMoveInterface>(cback.clone()) ); // we add unique clone
    }

    void MCI::callBackOnMove(const double x[], const bool accepted)
    {
        for (auto & cback : _cbacks){
            cback->callBackFunction(x, accepted);
        }
    }

    // setters

    void MCI::setTargetAcceptanceRate(const double targetaccrate)
    {
        _targetaccrate = targetaccrate;
    }

    void MCI::setMRT2Step(const double mrt2step)
    {
        for (int i=0; i<_trialMove->getNStepSizes(); ++i) {
            _trialMove->setStepSize(i, mrt2step);
        }
    }

    void MCI::setMRT2Step(const int i, const double mrt2step)
    {
        if (i<_trialMove->getNStepSizes()) {
        _trialMove->setStepSize(i, mrt2step);
        }
        else {
            std::cout << "[MCI::setMRT2Step] Warning: Tried to set non-existing MRT2step index." << std::endl;
        }
    }

    void MCI::setMRT2Step(const double mrt2step[])
    {
        for (int i=0; i<_trialMove->getNStepSizes(); ++i) {
            _trialMove->setStepSize(i, mrt2step[i]);
        }
    }


    void MCI::setX(const double x[])
    {
        std::copy(x, x+_ndim, _xold);
        applyPBC(_xold);
    }

    void MCI::setSeed(const uint_fast64_t seed) // fastest unsigned integer which is at least 64 bit (as expected by rgen)
    {
        _rgen.seed(seed);
    }

    // Domain
    void MCI::updateVolume()
    {
        // Set the integration volume
        _vol=1.;
        for (int i=0; i<_ndim; ++i) {
            _vol = _vol*( _ubound[i] - _lbound[i] );
        }
    }

    double MCI::getMinBoxLen() const {
        // find minimal box length
        double minboxlen = _ubound[0] - _lbound[0];
        for (int i=1; i<_ndim; ++i) {
            if ( (_ubound[i] - _lbound[i]) < minboxlen ) {
                minboxlen = _ubound[i] - _lbound[i];
            }
        }
        return minboxlen;
    }

    void MCI::applyPBC(double v[]) const
    {
        for (int i=0; i<_ndim; ++i) {
            while ( v[i] < _lbound[i] ) {
                v[i] += _ubound[i] - _lbound[i];
            }
            while ( v[i] > _ubound[i] ) {
                v[i] -= _ubound[i] - _lbound[i];
            }
        }
    }

    void MCI::setIRange(const double lbound, const double ubound)
    {
        // Set irange and apply PBC to the initial walker position _x
        std::fill(_lbound, _lbound+_ndim, lbound);
        std::fill(_ubound, _ubound+_ndim, ubound);
        updateVolume();

        applyPBC(_xold);
    }

    void MCI::setIRange(const double lbound[], const double ubound[])
    {
        // Set irange and apply PBC to the initial walker position _x
        std::copy(lbound, lbound+_ndim, _lbound);
        std::copy(ubound, ubound+_ndim, _ubound);
        updateVolume();

        applyPBC(_xold);
    }


    //   --- Constructor and Destructor

    MCI::MCI(const int ndim)
    {
        // _ndim
        _ndim = ndim;
        // _lbound and _ubound
        _lbound = new double[_ndim];
        _ubound = new double[_ndim];
        std::fill(_lbound, _lbound+_ndim, -0.49*std::numeric_limits<double>::max()); // set it so that ubound-lbound can still be computed
        std::fill(_ubound, _ubound+_ndim, 0.49*std::numeric_limits<double>::max());
        // _vol (will only be relevant without sampling function)
        _vol=0.;

        // _x
        _xold = new double[_ndim];
        std::fill(_xold, _xold+_ndim, 0.);
        _xnew = new double[_ndim];
        std::fill(_xnew, _xnew+_ndim, 0.);

        // _mrt2step
        //_mrt2step = new double[_ndim];
        //std::fill(_mrt2step, _mrt2step+_ndim, INITIAL_STEP);

        // initialize random generator
        _rgen = std::mt19937_64(_rdev());
        _rd = std::uniform_real_distribution<double>(0.,1.);

        // trial move
        _trialMove = std::unique_ptr<TrialMoveInterface>(new UniformAllMove(_ndim, INITIAL_STEP));
        _trialMove->bindRGen(_rgen);

        // other controls, defaulting to auto behavior
        _NfindMRT2Iterations = -1;
        _NdecorrelationSteps = -1;
        // initialize file flags
        _flagwlkfile=false;
        _flagobsfile=false;

        //initialize the running counters
        _ridx=0;
        _acc=0;
        _rej=0;
        // initialize all the other variables
        _targetaccrate=0.5;
        _flagMC=false;
    }

    MCI::~MCI()
    {
        // _mrt2step
        //delete [] _mrt2step;

        // _xold and _xnew
        delete [] _xnew;
        delete [] _xold;

        // lbound and ubound
        delete [] _ubound;
        delete [] _lbound;
    }

}  // namespace mci
