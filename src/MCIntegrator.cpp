#include "mci/MCIntegrator.hpp"

#include "mci/BlockAccumulator.hpp"
#include "mci/Estimators.hpp"
#include "mci/Factories.hpp"
#include "mci/FullAccumulator.hpp"
#include "mci/OrthoPeriodicDomain.hpp"
#include "mci/SRRDAllMove.hpp"
#include "mci/SRRDVecMove.hpp"
#include "mci/SimpleAccumulator.hpp"
#include "mci/UnboundDomain.hpp"

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
        if ( !_pdfcont.hasPDF() && !_domain->isFinite() ) {
            throw std::domain_error("[MCI::integrate] Integrating over an infinite domain requires a sampling function.");
        }

        if ( _pdfcont.hasPDF() ) {
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
            if (!_pdfcont.hasPDF()) {
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

    void MCI::findMRT2Step()
    {
        if (!_trialMove->hasStepSizes()) { return; } // in the odd case that our mover has no adjustable step sizes

        //constants
        const int nStepSizes = _trialMove->getNStepSizes();
        //const double changeRate = _trialMove->getChangeRate();
        const int64_t MIN_STAT=200LL*nStepSizes; // minimum statistic: number of M(RT)^2 steps done to decide on step size change
        const int MIN_CONS=5;   //minimum consecutive: minimum number of consecutive loops without need of changing mrt2step
        const double TOLERANCE=0.05;  //tolerance: tolerance for the acceptance rate
        const double SMALLEST_ACCEPTABLE_DOUBLE=std::numeric_limits<float>::min(); // use smallest float value as limit for double

        // fill temporary vectors
        std::vector<double> dimSizes(_ndim); // vector holding dimension sizes
        _domain->getSizes(dimSizes.data());

        std::vector<int> stepSizeIdx(_ndim); // mapping from x indices to used step size indices
        for (int i=0; i<_ndim; ++i) {
            stepSizeIdx[i] = _trialMove->getStepSizeIndex(i);
        }

        //initialize index
        int cons_count = 0;  //number of consecutive loops without need of changing mrt2step
        int counter = 0;  //counter of loops
        while ( ( _NfindMRT2Iterations < 0 && cons_count < MIN_CONS ) || counter < _NfindMRT2Iterations ) {
            //do MIN_STAT M(RT)^2 steps
            this->sample(MIN_STAT);

            //increase or decrease mrt2step depending on the acceptance rate
            const double rate = this->getAcceptanceRate();
            if ( fabs(rate-_targetaccrate) < TOLERANCE ) {
                ++cons_count; // acceptance was within tolerance
            } else {
                cons_count = 0; // we reset consecutive counter
            }

            const double fact = std::min(2., std::max(0.5, rate/_targetaccrate) );
            _trialMove->scaleStepSizes(fact); // scale move according to ratio

            // keep large step sizes in check
            for (int i=0; i<_ndim; ++i) {
                if (_trialMove->getStepSize(stepSizeIdx[i]) > 0.5*dimSizes[i]) {
                    _trialMove->setStepSize(stepSizeIdx[i], 0.5*dimSizes[i]);
                }
            }
            // keep small step sizes in check
            for (int j=0; j<nStepSizes; ++j) { //mrt2step ~ 0
                if ( _trialMove->getStepSize(j) < SMALLEST_ACCEPTABLE_DOUBLE ) {
                    _trialMove->setStepSize(j, SMALLEST_ACCEPTABLE_DOUBLE);
                }
            }

            ++counter;
            if ( _NfindMRT2Iterations < 0 && counter >= std::abs(_NfindMRT2Iterations) ) {
                std::cout << "Warning [MCI::findMRT2Step]: Max number of attempts reached without convergence." << std::endl;
                break;
            }
        }
    }


    void MCI::initialDecorrelation()
    {
        if (_NdecorrelationSteps < 0) {
            // automatic equilibration of contained observables with flag_equil = true

            //create the temporary observable container to be used
            ObservableContainer obs_equil;
            for (int i=0; i<_obscont.getNObs(); ++i) {
                if (_obscont.getFlagEquil(i)) {
                    obs_equil.addObservable(std::unique_ptr<AccumulatorInterface>( new FullAccumulator(_obscont.getObservableFunction(i).clone(), 1) ),
                                            createEstimator(EstimatorType::Correlated), true);
                }
            }
            const auto MIN_NMC=static_cast<int64_t>( sqrt(10000.*_ndim) ); // a guess on how many steps we need
            const int nobsdim = obs_equil.getNObsDim();
            // allocate memory for observables
            obs_equil.allocate(MIN_NMC);

            //do a first estimate of the observables
            this->sample(MIN_NMC, obs_equil, false);
            double oldestimate[nobsdim];
            double olderrestim[nobsdim];
            obs_equil.estimate(oldestimate, olderrestim);

            //start a loop which will stop when the observables are stabilized
            bool flag_loop=true;
            double newestimate[nobsdim];
            double newerrestim[nobsdim];
            int64_t countNMC = 0;
            while ( flag_loop ) {
                flag_loop = false;
                this->sample(MIN_NMC, obs_equil, false);
                countNMC += MIN_NMC;

                if (countNMC >= std::abs(_NdecorrelationSteps)) {
                    std::cout << "Warning [MCI::initialDecorrelation]: Max number of MC steps reached without equilibration." << std::endl;
                    break;
                }

                obs_equil.estimate(newestimate, newerrestim);
                for (int i=0; i<nobsdim; ++i) {
                    if ( fabs( oldestimate[i] - newestimate[i] ) > 2*sqrt( olderrestim[i]*olderrestim[i] + newerrestim[i]*newerrestim[i] ) ) {
                        flag_loop=true; // if any difference is too large, continue the loop
                        break; // break the inner for loop
                    }
                }
                // copy new to old
                std::copy(newestimate, newestimate+nobsdim, oldestimate);
                std::copy(newerrestim, newerrestim+nobsdim, olderrestim);
            }
            //memory deallocated automatically
        }
        else if (_NdecorrelationSteps > 0) {
            this->sample(_NdecorrelationSteps);
        }
    }


    // --- Sampling

    void MCI::initializeSampling(ObservableContainer * obsCont)
    {
        const bool flag_obs = (obsCont != nullptr);

        // reset running counters
        _acc = 0;
        _rej = 0;
        _ridx = 0;

        // init xnew and all protovalues
        _wlkstate.initialize(flag_obs);
        _pdfcont.initializeProtoValues(_wlkstate.xold); // initialize the pdf at x
        _trialMove->initializeProtoValues(_wlkstate.xold); // initialize the trial mover

        // init rest
        this->callBackOnMove(); // first call of the call-back functions
        if (flag_obs) {
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

    void MCI::doStepMRT2() // do MC step, sampling from _pdfcont
    {
        // propose a new position x and get move acceptance
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
            _pdfcont.newToOld();
            _trialMove->newToOld();
            _wlkstate.newToOld();
        } else { // rejected
            _pdfcont.oldToNew();
            _trialMove->oldToNew();
            _wlkstate.oldToNew();
        }
    }

    void MCI::doStepRandom() // do MC step, sampling randomly (used when _pdfcont is empty)
    {
        // set xnew to new random values within the domain
        for (int i=0; i<_ndim; ++i) { _wlkstate.xnew[i] = _rd(_rgen); } // between 0 and 1
        _domain->scaleToDomain(_wlkstate.xnew); // make it proper coordinates
        _wlkstate.nchanged = _ndim;

        // "accept" move
        _wlkstate.accepted = true;
        ++_acc;

        // rest
        this->callBackOnMove(); // call callbacks
        _wlkstate.newToOld(); // to mimic doStepMRT2()
    }

    // --- Domain

    std::unique_ptr<DomainInterface> MCI::setDomain(std::unique_ptr<DomainInterface> domain)
    {
        if (domain->ndim != _ndim) {
            throw std::invalid_argument("[MCI::setDomain] Passed domain's number of dimensions is not equal to MCI's number of walkers.");
        }
        std::swap(domain, _domain); // swap owned pointers
        _domain->applyDomain(_wlkstate.xold);
        return domain; // return old (will get deleted if not taken)
    }

    std::unique_ptr<DomainInterface> MCI::resetDomain() {
        std::unique_ptr<DomainInterface> domain( new UnboundDomain(_ndim) );
        std::swap(domain, _domain);
        _domain->applyDomain(_wlkstate.xold); // just in case
        return domain; // same as above
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

    std::unique_ptr<TrialMoveInterface> MCI::setTrialMove(std::unique_ptr<TrialMoveInterface> tmove)
    {
        if (tmove->getNDim() != _ndim) {
            throw std::invalid_argument("[MCI::setTrialMove] Passed trial move's number of inputs is not equal to MCI's number of walkers.");
        }
        std::swap(tmove, _trialMove); // unique ptr, old move gets freed automatically
        _trialMove->bindRGen(_rgen);
        return tmove; // deleted if not taken
    }

    std::unique_ptr<TrialMoveInterface> MCI::setTrialMove(MoveType move)
    {
        auto tmove = createMoveDefault(move, _ndim); // use factory default function
        std::swap(tmove, _trialMove);
        _trialMove->bindRGen(_rgen);
        return tmove;
    }

    std::unique_ptr<TrialMoveInterface> MCI::setTrialMove(SRRDType srrd, int veclen, int ntypes, int typeEnds[])
    {
        std::unique_ptr<TrialMoveInterface> tmove;
        if (veclen>0) {
            if (_ndim % veclen != 0) {
                throw std::invalid_argument("[MCI::setTrialMove] MCI's number of walkers must be a multiple of passed veclen.");
            }
            tmove = createSRRDVecMove(srrd, _ndim/veclen, veclen, ntypes, typeEnds);
        }
        else {
            tmove = createSRRDAllMove(srrd, _ndim, ntypes, typeEnds);
        }
        std::swap(tmove, _trialMove);
        _trialMove->bindRGen(_rgen);
        return tmove;
    }


    // --- Observables

    void MCI::addObservable(std::unique_ptr<ObservableFunctionInterface> obs, int blocksize, int nskip, const bool flag_equil, const EstimatorType estimType)
    {
        // sanity
        blocksize = std::max(0, blocksize);
        nskip = std::max(1, nskip);
        if (obs->getNDim() != _ndim) {
            throw std::invalid_argument("[MCI::addObservable] Passed observable function's number of inputs is not equal to MCI's number of walkers.");
        }
        if (flag_equil && estimType == EstimatorType::Noop) {
            throw std::invalid_argument("[MCI::addObservable] Requested automatic observable equilibration requires estimator with error calculation.");
        }

        // add accumulator&estimator from factory functions
        _obscont.addObservable(createAccumulator(std::move(obs), blocksize, nskip), createEstimator(estimType), flag_equil);
    }

    void MCI::addObservable(std::unique_ptr<ObservableFunctionInterface> obs, const int blocksize, const int nskip, const bool flag_equil, const bool flag_correlated)
    {
        // select type
        const bool flag_error = (blocksize > 0); // will we calculate errors?
        const EstimatorType estimType = selectEstimatorType(flag_correlated, flag_error);

        // use addObservable above
        this->addObservable(std::move(obs), blocksize, nskip, flag_equil, estimType);
    }

    std::unique_ptr<ObservableFunctionInterface> MCI::popObservable()
    {
        auto accu = _obscont.pop_back(); // remove accu from container
        auto obs = accu->removeObs(); // extract contained observable, accu will be destroyed after return
        return obs; // return obs
    }

    // --- Sampling functions

    void MCI::addSamplingFunction(std::unique_ptr<SamplingFunctionInterface> pdf)
    {
        if (pdf->getNDim() != _ndim) {
            throw std::invalid_argument("[MCI::addSamplingFunction] Passed sampling function's number of inputs is not equal to MCI's number of walkers.");
        }
        _pdfcont.addSamplingFunction(std::move(pdf)); // we move pdf into pdfcont
    }


    // --- Callbacks

    void MCI::addCallBack(std::unique_ptr<CallBackOnMoveInterface> cback)
    {
        if (cback->getNDim() != _ndim) {
            throw std::invalid_argument("[MCI::addCallBack] Passed callback function's number of inputs is not equal to MCI's number of walkers.");
        }
        _cbacks.emplace_back( std::move(cback) ); // we move cback ptr into a vector
    }

    std::unique_ptr<CallBackOnMoveInterface> MCI::popCallBack()
    {
        auto cback = std::move(_cbacks.back()); // move last element out of vector
        _cbacks.pop_back(); // resize vec
        return cback;
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
            _obsfile << _ridx;
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

    void MCI::setX(const int i, const double val)
    {
        _wlkstate.xold[i] = val;
        _domain->applyDomain(_wlkstate.xold);
    }

    void MCI::setX(const double x[])
    {
        std::copy(x, x+_ndim, _wlkstate.xold);
        _domain->applyDomain(_wlkstate.xold);
    }

    void MCI::moveX() // for external user, to manually use trialMove on xold
    {
        _wlkstate.oldToNew(); // our trial moves expect proper xnew
        _trialMove->computeTrialMove(_wlkstate);
        _domain->applyDomain(_wlkstate);
        _wlkstate.newToOld(); // but here we want to set xold
    }

    void MCI::newRandomX() // also meant for the user
    {
        for (int i=0; i<_ndim; ++i) { _wlkstate.xnew[i] = _rd(_rgen); } // draw random numbers between 0 and 1
        _domain->scaleToDomain(_wlkstate.xnew); // shift/scale to proper domain coordinates
        _wlkstate.newToOld();
    }

    //   --- Constructor and Destructor

    MCI::MCI(const int ndim): _ndim(ndim), _wlkstate(_ndim, false)
    {
        // initialize random generator
        _rgen = std::mt19937_64(_rdev()); // passed through to trial moves (for seed consistency)
        _rd = std::uniform_real_distribution<double>(0.,1.); // used for acceptance (and random moves without pdf)

        // domain
        this->resetDomain(); // default to unbound domain

        // trial move (will be used when sampling with pdf)
        this->setTrialMove(MoveType::All); // default to uniform all-move

        // other controls, defaulting to auto behavior
        _targetaccrate=0.5;
        _NfindMRT2Iterations = -50; // default to max 50 auto-iterations
        _NdecorrelationSteps = -10000; // default to max 10k auto-steps

        // initialize file flags
        _flagwlkfile=false;
        _flagobsfile=false;

        //initialize the running counters
        _ridx=0;
        _acc=0;
        _rej=0;
    }

}  // namespace mci
