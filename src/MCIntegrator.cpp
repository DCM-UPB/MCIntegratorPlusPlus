#include "mci/MCIntegrator.hpp"

#include "mci/BlockAccumulator.hpp"
#include "mci/Estimators.hpp"
#include "mci/Factories.hpp"
#include "mci/FullAccumulator.hpp"
#include "mci/SRRDAllMove.hpp"
#include "mci/SRRDVecMove.hpp"
#include "mci/SimpleAccumulator.hpp"
#include "mci/UnboundDomain.hpp"
#include "mci/OrthoPeriodicDomain.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace mci
{

    //  --- Integrate

    void MCI::integrate(const int64_t Nmc, double average[], double error[], const bool doFindMRT2step, const bool doDecorrelation)
    {
        if ( _pdfcont.empty() && !_domain->isFinite() ) {
            throw std::domain_error("[MCI::integrate] Integrating over an infinite domain requires a sampling function.");
        }

        if ( !_pdfcont.empty() ) {
            //find the optimal mrt2 step
            if (doFindMRT2step) { this->findMRT2Step(); }
            // take care to do the initial decorrelation of the walker
            if (doDecorrelation) { this->initialDecorrelation(); }
        }

        if (Nmc>0) {
            // allocation of the accumulators where the data will be stored
            _obscont.allocate(Nmc);

            //sample the observables
            if (_flagobsfile) { _obsfile.open(_pathobsfile); }
            if (_flagwlkfile) { _wlkfile.open(_pathwlkfile); }
            this->sample(Nmc, _obscont, true); // let sample accumulate data
            if (_flagobsfile) { _obsfile.close(); }
            if (_flagwlkfile) { _wlkfile.close(); }

            // estimate average and standard deviation
            _obscont.estimate(average, error);

            // if we sampled randomly, scale results by volume
            if (_pdfcont.empty()) {
                const double vol = _domain->getVolume();
                for (int i=0; i<_obscont.getNObsDim(); ++i) {
                    average[i] *= vol;
                    error[i] *= vol;
                }
            }

            // deallocate
            _obscont.deallocate();
        }
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
            const int64_t MIN_NMC=100;
            const int nobsdim = obs_equil.getNObsDim();
            // allocate memory for observables
            obs_equil.allocate(MIN_NMC);

            //do a first estimate of the observables
            this->sample(MIN_NMC, obs_equil, false);
            auto * oldestimate = new double[nobsdim];
            auto * olderrestim = new double[nobsdim];
            obs_equil.estimate(oldestimate, olderrestim);

            //start a loop which will stop when the observables are stabilized
            bool flag_loop=true;
            auto * newestimate = new double[nobsdim];
            auto * newerrestim = new double[nobsdim];
            while ( flag_loop ) {
                flag_loop = false;
                this->sample(MIN_NMC, obs_equil, false);
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
        const int64_t MIN_STAT=200*_trialMove->getNStepSizes(); // minimum statistic: number of M(RT)^2 steps done to decide on step size change
        const int MIN_CONS=5;   //minimum consecutive: minimum number of consecutive loops without need of changing mrt2step
        const double TOLERANCE=0.05;  //tolerance: tolerance for the acceptance rate
        const int MAX_NUM_ATTEMPTS=50;  //maximum number of attempts: maximum number of time that the main loop can be executed
        const double SMALLEST_ACCEPTABLE_DOUBLE=1.e-50;

        //initialize index
        int cons_count = 0;  //number of consecutive loops without need of changing mrt2step
        int counter = 0;  //counter of loops
        while ( ( _NfindMRT2Iterations < 0 && cons_count < MIN_CONS ) || counter < _NfindMRT2Iterations ) {
            const int nsizes = _trialMove->getNStepSizes();
            const double maxstepsize = _domain->getMaxStepSize();

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
                    if ( _trialMove->getStepSize(j) > maxstepsize ) {
                        _trialMove->setStepSize(j, maxstepsize);
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

    void MCI::initializeSampling(ObservableContainer * obsCont)
    {
        // reset running counters
        _acc = 0;
        _rej = 0;
        _ridx = 0;

        // init xnew and all protovalues
        _wlkstate.initialize();
        _pdfcont.initializeProtoValues(_wlkstate.xold); // initialize the pdf at x
        _trialMove->initializeProtoValues(_wlkstate.xold); // initialize the trial mover

        // init rest
        this->callBackOnMove(); // first call of the call-back functions
        if (obsCont != nullptr) { // optional passed observable container
            obsCont->reset(); // reset observable accumulators
        }
    }

    void MCI::sample(const int64_t npoints) // sample without taking observables or printing to file
    {
        // Initialize
        this->initializeSampling(nullptr);

        //run the main loop for sampling
        const bool flagpdf = _pdfcont.hasPDF();
        for (_ridx=0; _ridx<npoints; ++_ridx) {
            if (flagpdf) { // use sampling function
                this->doStepMRT2();
            }
            else { // sample randomly
                this->doStepRandom();
            }
        }
    }

    void MCI::sample(const int64_t npoints, ObservableContainer &container, const bool flagMC)
    {
        // Initialize
        this->initializeSampling(&container);

        //run the main loop for sampling
        const bool flagpdf = _pdfcont.hasPDF();
        for (_ridx=0; _ridx<npoints; ++_ridx) {
            // do MC step
            if (flagpdf) { // use sampling function
                this->doStepMRT2();
            }
            else { // sample randomly
                this->doStepRandom();
            }

            // accumulate obs
            container.accumulate(_wlkstate);

            // file output
            if (flagMC && _flagobsfile) { this->storeObservables(); } // store obs on file
            if (flagMC && _flagwlkfile) { this->storeWalkerPositions(); } // store walkers on file
        }

        // finalize data
        container.finalize();
    }


    // --- Walking

    void MCI::doStepMRT2()
    {
        // propose a new position x and get move acceptance
        _wlkstate.nchanged=0; // better don't rely on this!
        const double moveAcc = _trialMove->computeTrialMove(_wlkstate);

        // apply PBC update
        if (_wlkstate.nchanged < _ndim) {
            _domain->applyDomain(_wlkstate); // selective update
        } else {
            _domain->applyDomain(_wlkstate.xnew);
        }

        // find the corresponding sampling function acceptance
        const double pdfAcc = _pdfcont.computeAcceptance(_wlkstate);

        // determine if the proposed x is accepted or not
        _wlkstate.accepted = ( _rd(_rgen) <= pdfAcc*moveAcc );
        _wlkstate.accepted ? ++_acc : ++_rej; // increase counters

        // call callbacks
        this->callBackOnMove();

        // set state according to result
        if (_wlkstate.accepted) {
            _wlkstate.newToOld();
            _pdfcont.newToOld();
            _trialMove->newToOld();
            _trialMove->callOnAcceptance(_pdfcont); // expects pdfcont swapped already
        } else { // rejected
            _wlkstate.oldToNew();
            _pdfcont.oldToNew();
            _trialMove->oldToNew();
        }
    }

    void MCI::doStepRandom()
    {
        // set xnew to new random values (within the irange)
        // use fixed uniform all-moves on maxStep*[-0.5, 0.5],
        // which means effectively uncorrelated samples
        const double stepSize = _domain->getMaxStepSize();
        for (int i=0; i<_ndim; ++i) {
            _wlkstate.xnew[i] += stepSize*(_rd(_rgen) - 0.5);
        }
        _wlkstate.nchanged = _ndim;
        _domain->applyDomain(_wlkstate.xnew);

        // "accept" move
        _wlkstate.accepted = true;
        ++_acc;

        // rest
        this->callBackOnMove(); // call callbacks
        _wlkstate.newToOld(); // to mimic doStepMRT2()
    }

    // --- Domain

    void MCI::setDomain(const DomainInterface &domain)
    {
        if (domain.ndim != _ndim) {
            throw std::invalid_argument("[MCI::setDomain] Passed domain's number of dimensions is not equal to MCI's number of walkers.");
        }
        _domain = domain.clone();
    }

    void MCI::resetDomain() {
        _domain = std::unique_ptr<DomainInterface>(new UnboundDomain(_ndim));
    }

    void MCI::setIRange(const double lbound, const double ubound)
    {
        _domain = std::unique_ptr<DomainInterface>(new OrthoPeriodicDomain(_ndim, lbound, ubound));
        _domain->applyDomain(_wlkstate.xold);
    }

    void MCI::setIRange(const double lbounds[], const double ubounds[])
    {
        _domain = std::unique_ptr<DomainInterface>(new OrthoPeriodicDomain(_ndim, lbounds, ubounds));
        _domain->applyDomain(_wlkstate.xold);
    }


    // --- Trial Moves

    void MCI::setTrialMove(const TrialMoveInterface &tmove)
    {
        if (tmove.getNDim() != _ndim) {
            throw std::invalid_argument("[MCI::setTrialMove] Passed trial move's number of inputs is not equal to MCI's number of walkers.");
        }
        _trialMove = tmove.clone(); // unique ptr, old move gets freed automatically
        _trialMove->bindRGen(_rgen);
    }

    void MCI::setTrialMove(MoveType move)
    {
        _trialMove = createMoveDefault(move, _ndim); // use factory default function
        _trialMove->bindRGen(_rgen);
    }

    void MCI::setTrialMove(SRRDType srrd, int veclen, int ntypes, int typeEnds[])
    {
        if (veclen>0) {
            if (_ndim % veclen != 0) {
                throw std::invalid_argument("[MCI::setTrialMove] MCI's number of walkers must be a multiple of passed veclen.");
            }
            _trialMove = createSRRDVecMove(srrd, _ndim/veclen, veclen, ntypes, typeEnds);
        }
        else {
            _trialMove = createSRRDAllMove(srrd, _ndim, ntypes, typeEnds);
        }

        _trialMove->bindRGen(_rgen);
    }


    // --- Observables

    void MCI::clearObservables()
    {
        _obscont.clear();
    }

    void MCI::addObservable(const ObservableFunctionInterface &obs, int blocksize, int nskip, const bool flag_equil, const EstimatorType estimType)
    {
        // sanity
        blocksize = std::max(0, blocksize);
        nskip = std::max(1, nskip);
        if (obs.getNDim() != _ndim) {
            throw std::invalid_argument("[MCI::addObservable] Passed observable function's number of inputs is not equal to MCI's number of walkers.");
        }
        if (flag_equil && estimType == EstimatorType::Noop) {
            throw std::invalid_argument("[MCI::addObservable] Requested automatic observable equilibration requires estimator with error calculation.");
        }

        // add accumulator&estimator from factory functions
        _obscont.addObservable(createAccumulator(obs, blocksize, nskip), createEstimator(estimType), flag_equil);
    }

    void MCI::addObservable(const ObservableFunctionInterface &obs, const int blocksize, const int nskip, const bool flag_equil, const bool flag_correlated)
    {
        // select type
        const bool flag_error = (blocksize > 0); // will we calculate errors?
        const EstimatorType estimType = selectEstimatorType(flag_correlated, flag_error);

        // use addObservable above
        this->addObservable(obs, blocksize, nskip, flag_equil, estimType);
    }


    // --- Sampling functions

    void MCI::clearSamplingFunctions()
    {
        _pdfcont.clear();
    }

    void MCI::addSamplingFunction(const SamplingFunctionInterface &mcisf)
    {
        if (mcisf.getNDim() != _ndim) {
            throw std::invalid_argument("[MCI::addSamplingFunction] Passed sampling function's number of inputs is not equal to MCI's number of walkers.");
        }
        _pdfcont.addSamplingFunction( mcisf.clone() );
    }


    // --- Callbacks

    void MCI::clearCallBacks() {
        _cbacks.clear();
    }

    void MCI::addCallBack(const CallBackOnMoveInterface &cback)
    {
        if (cback.getNDim() != _ndim) {
            throw std::invalid_argument("[MCI::addCallBack] Passed callback function's number of inputs is not equal to MCI's number of walkers.");
        }
        _cbacks.emplace_back( std::unique_ptr<CallBackOnMoveInterface>(cback.clone()) ); // we add unique clone
    }

    void MCI::callBackOnMove()
    {
        for (auto & cback : _cbacks){
            cback->callBackFunction(_wlkstate);
        }
    }


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
                _wlkfile << "   " << _wlkstate.xold[j] ;
            }
            _wlkfile << std::endl;
        }
    }


    // --- Setters

    void MCI::setSeed(const uint_fast64_t seed) // fastest unsigned integer which is at least 64 bit (as expected by rgen)
    {
        _rgen.seed(seed);
    }

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
        std::copy(x, x+_ndim, _wlkstate.xold);
        _domain->applyDomain(_wlkstate.xold);
    }

    void MCI::moveX() // for external user, to manually use trialMove on xold
    {
        _trialMove->computeTrialMove(_wlkstate);
        _wlkstate.newToOld();
        _domain->applyDomain(_wlkstate.xold);
    }

    //   --- Constructor and Destructor

    MCI::MCI(const int ndim): _ndim(ndim), _wlkstate(_ndim)
    {
        // initialize random generator
        _rgen = std::mt19937_64(_rdev()); // passed through to trial moves (for seed consistency)
        _rd = std::uniform_real_distribution<double>(0.,1.); // used for acceptance (and random moves without pdf)

        // domain
        this->resetDomain();

        // default trial move (will be used when sampling with pdf)
        this->setTrialMove(MoveType::All);

        // other controls, defaulting to auto behavior
        _targetaccrate=0.5;
        _NfindMRT2Iterations = -1;
        _NdecorrelationSteps = -1;

        // initialize file flags
        _flagwlkfile=false;
        _flagobsfile=false;

        //initialize the running counters
        _ridx=0;
        _acc=0;
        _rej=0;
    }

}  // namespace mci
