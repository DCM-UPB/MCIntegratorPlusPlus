#include "mci/MCIntegrator.hpp"

#include "mci/BlockAccumulator.hpp"
#include "mci/Estimators.hpp"
#include "mci/FullAccumulator.hpp"
#include "mci/SimpleAccumulator.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace mci
{

    //   --- Integrate
    void MCI::integrate(const int Nmc, double average[], double error[], const bool doFindMRT2step, const bool doDecorrelation)
    {
        if ( _flagpdf ) {
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
        this->sample(Nmc, &_obscont); // let sample accumulate data
        _flagMC = false;
        if (_flagobsfile) { _obsfile.close(); }
        if (_flagwlkfile) { _wlkfile.close(); }

        // estimate average and standard deviation
        _obscont.estimate(average, error);

        // if we sampled uniformly, scale results by volume
        if (!_flagpdf) {
            for (int i=0; i<_obscont.getNObsDim(); ++i) {
                average[i] *=_vol;
                error[i] *=_vol;
            }
        }

        // deallocate
        _obscont.deallocate();
    }



    //   --- Internal methods


    void MCI::storeObservables()
    {
        if ( _ridx%_freqobsfile == 0 ) {
            _obscont.printObsValues(_obsfile);
            _obsfile << std::endl;
        }
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


    void MCI::sample(const int npoints, ObservableContainer * container)
    {
        const bool flagobs = (container != nullptr); // should we sample observables ?

        // reset acceptance and rejection
        this->resetAccRejCounters();

        // reset observable data (to be sure)
        if (flagobs) { container->reset(); }

        // first call of the call-back functions
        if (_flagpdf) {
            for (auto & cback : _cbacks){
                cback->callBackFunction(_xold, true);
            }
        }

        // initialize the pdf at x
        computeOldSamplingFunction();

        //run the main loop for sampling
        bool flags_xchanged[_ndim]; // will remember which x elements changed
        for (_ridx=0; _ridx<npoints; ++_ridx) {
            std::fill(flags_xchanged, flags_xchanged+_ndim, false); // reset flags

            bool flagacc;
            if (_flagpdf) { // use sampling function
                flagacc = this->doStepMRT2(flags_xchanged);
            }
            else {
                this->newRandomX();
                ++_acc; // "accept" move
                flagacc = true;
                std::fill(flags_xchanged, flags_xchanged+_ndim, true);
            }

            if (flagobs) {
                container->accumulate(_xold, flagacc, flags_xchanged); // accumulate observables
                if (_flagMC && _flagobsfile) { this->storeObservables(); } // store obs on file
            }

            if (_flagMC && _flagwlkfile) { this->storeWalkerPositions(); } // store walkers on file
        }

        // finalize data
        if (flagobs) { container->finalize(); }
    }


    void MCI::initialDecorrelation()
    {
        if (_NdecorrelationSteps < 0) {
            // automatic equilibration of contained observables with flag_equil = true

            //create the temporary observable container to be used
            ObservableContainer obs_equil;
            for (int i=0; i<_obscont.getNObs(); ++i) {
                if (_obscont.getFlagEquil(i)) {
                    obs_equil.addObservable(std::unique_ptr<AccumulatorInterface>( new FullAccumulator(*_obscont.getObservableFunction(i).clone(), 1) ),
                                            mci::CorrelatedEstimator, true);
                }
            }
            const int MIN_NMC=100;
            const int nobsdim = obs_equil.getNObsDim();
            // allocate memory for observables
            obs_equil.allocate(MIN_NMC);

            //do a first estimate of the observables
            this->sample(MIN_NMC, &obs_equil);
            auto * oldestimate = new double[nobsdim];
            auto * olderrestim = new double[nobsdim];
            obs_equil.estimate(oldestimate, olderrestim);

            //start a loop which will stop when the observables are stabilized
            bool flag_loop=true;
            auto * newestimate = new double[nobsdim];
            auto * newerrestim = new double[nobsdim];
            while ( flag_loop ) {
                flag_loop = false;
                this->sample(MIN_NMC, &obs_equil);
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
            this->sample(_NdecorrelationSteps, nullptr);
        }
    }


    void MCI::findMRT2Step()
    {
        //constants
        const int MIN_STAT=200;  //minimum statistic: number of M(RT)^2 steps done before deciding if the step must be increased or decreased
        const int MIN_CONS=5;   //minimum consecutive: minimum number of consecutive loops without need of changing mrt2step
        const double TOLERANCE=0.05;  //tolerance: tolerance for the acceptance rate
        const int MAX_NUM_ATTEMPTS=50;  //maximum number of attempts: maximum number of time that the main loop can be executed
        const double SMALLEST_ACCEPTABLE_DOUBLE=1.e-50;

        //initialize index
        int cons_count = 0;  //number of consecutive loops without need of changing mrt2step
        int counter = 0;  //counter of loops
        double fact;
        while ( ( _NfindMRT2Iterations < 0 && cons_count < MIN_CONS ) || counter < _NfindMRT2Iterations ) {
            //do MIN_STAT M(RT)^2 steps
            this->sample(MIN_STAT, nullptr);

            //increase or decrease mrt2step depending on the acceptance rate
            double rate = this->getAcceptanceRate();
            if ( fabs(rate-_targetaccrate) < TOLERANCE ) {
                //mrt2step was ok
                cons_count++;
            }
            else {
                //need to change mrt2step
                cons_count=0;
                fact = std::min(2., std::max(0.5, rate/_targetaccrate) );
                for (int j=0; j<_ndim; ++j) {
                    _mrt2step[j] *= fact;
                }

                // sanity checks
                for (int j=0; j<_ndim; ++j) { //mrt2step = Infinity
                    if ( _mrt2step[j] > ( _ubound[j] - _lbound[j] ) ) {
                        _mrt2step[j] = ( _ubound[j] - _lbound[j] );
                    }
                }
                for (int j=0; j<_ndim; ++j) { //mrt2step ~ 0
                    if ( _mrt2step[j] < SMALLEST_ACCEPTABLE_DOUBLE ) {
                        _mrt2step[j] = SMALLEST_ACCEPTABLE_DOUBLE;
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


    bool MCI::doStepMRT2(bool * flags_xchanged)
    {
        // propose a new position x
        this->computeNewX();

        // find the corresponding sampling function value
        this->computeNewSamplingFunction();

        //determine if the proposed x is accepted or not
        const bool flagacc = ( _rd(_rgen) <= this->computeAcceptance() );

        //update some values according to the acceptance of the mrt2 step
        if ( flagacc ) {
            std::fill(flags_xchanged, flags_xchanged+_ndim, true); // currently we do all-particle steps
            //accepted
            _acc++;
            //update the walker position x
            this->updateX();
            //update the sampling function values pdfx
            this->updateSamplingFunction();
            //if there are some call back functions, invoke them
            for (auto & cback : _cbacks){
                cback->callBackFunction(_xold, _flagMC);
            }
        } else {
            //rejected
            _rej++;
        }

        return flagacc;
    }


    void MCI::updateVolume()
    {
        // Set the integration volume
        _vol=1.;
        for (int i=0; i<_ndim; ++i) {
            _vol = _vol*( _ubound[i] - _lbound[i] );
        }
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


    void MCI::updateX()
    {
        std::swap(_xold, _xnew);
    }


    void MCI::newRandomX()
    {
        //set xold to new random values (within the irange)
        for (int i=0; i<_ndim; ++i) {
            _xold[i] = _lbound[i] + ( _ubound[i] - _lbound[i] ) * _rd(_rgen);
        }
    }


    void MCI::resetAccRejCounters()
    {
        _acc = 0;
        _rej = 0;
    }


    void MCI::computeNewX()
    {
        for (int i=0; i<_ndim; ++i) {
            _xnew[i] = _xold[i] + _mrt2step[i] * (_rd(_rgen)-0.5);
        }
        applyPBC(_xnew);
    }


    void MCI::updateSamplingFunction()
    {
        for (auto & sf : _pdfs) {
            sf->newToOld();
        }
    }


    double MCI::computeAcceptance() const
    {
        double acceptance=1.;
        for (auto & sf : _pdfs) {
            acceptance*=sf->getAcceptance();
        }
        return acceptance;
    }


    void MCI::computeOldSamplingFunction()
    {
        for (auto & sf : _pdfs) {
            sf->computeNewSamplingFunction(_xold);
            sf->newToOld();
        }
    }


    void MCI::computeNewSamplingFunction()
    {
        for (auto & sf : _pdfs) {
            sf->computeNewSamplingFunction(_xnew);
        }
    }


    //   --- Setters


    void MCI::storeObservablesOnFile(const std::string &filepath, const int freq)
    {
        _pathobsfile = filepath;
        _freqobsfile = freq;
        _flagobsfile = true;
    }


    void MCI::storeWalkerPositionsOnFile(const std::string &filepath, const int freq)
    {
        _pathwlkfile = filepath;
        _freqwlkfile = freq;
        _flagwlkfile = true;
    }


    void MCI::clearCallBackOnMove(){
        _cbacks.clear();
    }


    void MCI::addCallBackOnMove(const CallBackOnMoveInterface &cback){
        _cbacks.emplace_back( std::unique_ptr<CallBackOnMoveInterface>(cback.clone()) ); // we add unique clone
    }


    void MCI::clearObservables()
    {
        _obscont.clear();
    }


    void MCI::addObservable(const ObservableFunctionInterface &obs, int blocksize, int nskip, const bool flag_equil, const bool flag_correlated)
    {
        // sanity
        blocksize = std::max(0, blocksize);
        nskip = std::max(1, nskip);
        if (flag_equil && blocksize==0) {
            throw std::invalid_argument("[MCI::addObservable] Requested automatic observable equilibration requires blocksize > 0.");
        }
        if (flag_equil && blocksize==0) {
            throw std::invalid_argument("[MCI::addObservable] Requested correlated error estimation requires blocksize > 0.");
        }

        // we need to select these two
        std::unique_ptr<AccumulatorInterface> accu;
        std::function< void(int /*nstore*/, int /*nobs*/, const double [] /*data*/, double [] /*avg*/, double [] /*error*/) > estim;

        if (blocksize == 0) {
            accu = std::unique_ptr<AccumulatorInterface>( new SimpleAccumulator(obs, nskip) );
            // data is already the average, so this estimator just copies the average and fills error with 0
            estim = [](int /*unused*/, int nobs, const double data[], double avg[], double err[]) {
                        std::copy(data, data+nobs, avg);
                        std::fill(err, err+nobs, 0.);
                    };
        } else {
            if (blocksize == 1) {
                accu = std::unique_ptr<AccumulatorInterface>( new FullAccumulator(obs, nskip) );
            } else {
                accu = std::unique_ptr<AccumulatorInterface>( new BlockAccumulator(obs, nskip, blocksize) );
            }
            estim = flag_correlated ? mci::CorrelatedEstimator : mci::UncorrelatedEstimator ;
        }

        // append to container
        _obscont.addObservable(std::move(accu), estim, flag_equil);
    }


    void MCI::clearSamplingFunctions()
    {
        _pdfs.clear();
        _flagpdf = false;
    }


    void MCI::addSamplingFunction(const SamplingFunctionInterface &mcisf)
    {
        _pdfs.emplace_back( std::unique_ptr<SamplingFunctionInterface>(mcisf.clone()) ); // we add unique clone
        _flagpdf = true;
    }


    void MCI::setTargetAcceptanceRate(const double targetaccrate)
    {
        _targetaccrate = targetaccrate;
    }


    void MCI::setMRT2Step(const double mrt2step[])
    {
        std::copy(mrt2step, mrt2step+_ndim, _mrt2step);
    }


    void MCI::setX(const double x[])
    {
        std::copy(x, x+_ndim, _xold);
        applyPBC(_xold);
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

    void MCI::setSeed(const uint_fast64_t seed) // fastest unsigned integer which is at least 64 bit (as expected by rgen)
    {
        _rgen.seed(seed);
    }


    //   --- Constructor and Destructor

    MCI::MCI(const int ndim)
    {
        // _ndim
        _ndim = ndim;
        // _lbound and _ubound
        _lbound = new double[_ndim];
        _ubound = new double[_ndim];
        std::fill(_lbound, _lbound+_ndim, -std::numeric_limits<double>::max());
        std::fill(_ubound, _ubound+_ndim, std::numeric_limits<double>::max());
        // _vol (will only be relevant without sampling function)
        _vol=0.;

        // _x
        _xold = new double[_ndim];
        std::fill(_xold, _xold+_ndim, 0.);
        _xnew = new double[_ndim];
        std::fill(_xnew, _xnew+_ndim, 0.);

        // _mrt2step
        _mrt2step = new double[_ndim];
        std::fill(_mrt2step, _mrt2step+_ndim, INITIAL_STEP);

        // other controls, defaulting to auto behavior
        _NfindMRT2Iterations = -1;
        _NdecorrelationSteps = -1;

        // probability density function
        _flagpdf = false;
        // initialize file flags
        _flagwlkfile=false;
        _flagobsfile=false;
        // initialize random generator
        _rgen = std::mt19937_64(_rdev());
        _rd = std::uniform_real_distribution<double>(0.,1.);
        //initialize the running index
        _ridx=0;
        //initialize the acceptance counters
        _acc=0;
        _rej=0;
        // initialize all the other variables
        _targetaccrate=0.5;
        _flagMC=false;
    }

    MCI::~MCI()
    {
        // _mrt2step
        delete [] _mrt2step;

        // _xold and _xnew
        delete [] _xnew;
        delete [] _xold;

        // lbound and ubound
        delete [] _ubound;
        delete [] _lbound;
    }

}  // namespace mci
