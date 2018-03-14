#include "MCIntegrator.hpp"
#include "Estimators.hpp"

#include <limits>
#include <iostream>
#include <algorithm>
#include <vector>
#include <random>


//   --- Integrate

void MCI::integrate(const long &Nmc, double * average, double * error, bool findMRT2step, bool initialdecorrelation)
{
    long i;

    if ( _flagpdf )
    {
        //find the optimal mrt2 step
        if (findMRT2step) this->findMRT2Step();
        // take care to do the initial decorrelation of the walker
        if (initialdecorrelation) this->initialDecorrelation();
    }

    //allocation of the array where the data will be stored
    _datax = new double*[Nmc];
    for (i=0; i<Nmc; ++i){ *(_datax+i) = new double[_nobsdim]; }

    //sample the observables
    if (_flagobsfile) _obsfile.open(_pathobsfile);
    if (_flagwlkfile) _wlkfile.open(_pathwlkfile);
    _flagMC = true;
    this->sample(Nmc, true);
    _flagMC = false;
    if (_flagobsfile) _obsfile.close();
    if (_flagwlkfile) _wlkfile.close();

    //estimate average and standard deviation
    if ( _flagpdf )
        {
            mci::MultiDimCorrelatedEstimator(Nmc, _nobsdim, _datax, average, error);
        }
    else
        {
            mci::MultiDimUncorrelatedEstimator(Nmc, _nobsdim, _datax, average, error);
            for (i=0; i<_nobsdim; ++i)
                {
                    *(average+i) *=_vol; *(error+i) *=_vol;
                }
        }

    //deallocation of the data array
    for (int i=0; i<Nmc; ++i){ delete [] *(_datax+i); }
    delete [] _datax;
}



//   --- Internal methods


void MCI::storeObservables()
{
    if ( _ridx%_freqobsfile == 0 )
        {
            _obsfile << _ridx;
            for (int j=0; j<_nobsdim; ++j)
                {
                    _obsfile << "   " << *(*(_datax+_ridx)+j) ;
                }
            _obsfile << std::endl;
        }
}


void MCI::storeWalkerPositions()
{
    if ( _ridx%_freqwlkfile == 0 )
        {
            _wlkfile << _ridx;
            for (int j=0; j<_ndim; ++j)
                {
                    _wlkfile << "   " << *(_xold+j) ;
                }
            _wlkfile << std::endl;
        }
}


void MCI::initialDecorrelation()
{
    long i;
    //constants
    const long MIN_NMC=100;
    //allocate the data array that will be used
    _datax = new double*[MIN_NMC];
    for (i=0; i<MIN_NMC; ++i) {*(_datax+i) = new double[_nobsdim];   }
    //do a first estimate of the observables
    this->sample(MIN_NMC, true);
    double * oldestimate = new double[_nobsdim];
    double * olderrestim = new double[_nobsdim];
    mci::MultiDimCorrelatedEstimator(MIN_NMC, _nobsdim, _datax, oldestimate, olderrestim);
    //start a loop which will stop when the observables are stabilized
    bool flag_loop=true;
    double * newestimate = new double[_nobsdim];
    double * newerrestim = new double[_nobsdim];
    double * foo;
    while ( flag_loop )
        {
            this->sample(MIN_NMC, true);
            mci::MultiDimCorrelatedEstimator(MIN_NMC, _nobsdim, _datax, newestimate, newerrestim);
            for (i=0; i<_nobsdim; ++i)
                {
                    if ( abs( *(oldestimate+i) - *(newestimate+i) ) <= *(olderrestim+i) + *(newerrestim+i) ) flag_loop=false;
                    if (flag_loop) *(oldestimate+i) = *(newestimate+i);
                }
            foo=oldestimate;
            oldestimate=newestimate;
            newestimate=foo;
        }
    //memory deallocation
    delete [] newestimate;
    delete [] newerrestim;
    delete [] oldestimate;
    delete [] olderrestim;
    for (i=0; i<MIN_NMC; ++i){ delete[] *(_datax+i); }
    delete [] _datax;
}


void MCI::sample(const long &npoints, const bool &flagobs)
{
    int i;
    //initialize the running index
    _ridx=0;
    //initialize the pdf at x
    computeOldSamplingFunction();
    //initialize the observables values at x
    if (flagobs)
        {
            this->computeObservables();
        }
    //reset acceptance and rejection
    this->resetAccRejCounters();
    //first call of the call-back functions
    if (flagobs){
        for (MCICallBackOnAcceptanceInterface * cback : _cback){
            cback->callBackFunction(_xold, true);
        }
    }
    //start the main loop for sampling
    if ( _flagpdf )
        {
            bool flagacc;
            if ( flagobs )
                {
                    for (i=0; i<npoints; ++i)
                        {
                            this->doStepMRT2(&flagacc);
                            if (flagacc)
                                {
                                    this->computeObservables();
                                    this->saveObservables();
                                } else
                                {
                                    this->saveObservables();
                                }
                            if (_flagMC){
                                if (_flagobsfile) this->storeObservables();
                                if (_flagwlkfile) this->storeWalkerPositions();
                            }
                            _ridx++;
                        }
                } else
                {
                    for (i=0; i<npoints; ++i)
                        {
                            this->doStepMRT2(&flagacc);
                            _ridx++;
                        }
                }
        }
    else
        {
            if ( flagobs )
                {
                    for (i=0; i<npoints; ++i)
                        {
                            this->newRandomX();
                            this->computeObservables();
                            this->saveObservables();
                            if (_flagMC){
                                if (_flagobsfile) this->storeObservables();
                                if (_flagwlkfile) this->storeWalkerPositions();
                            }
                            _ridx++;
                        }
                } else
                {
                    for (i=0; i<npoints; ++i)
                        {
                            this->newRandomX();
                            _ridx++;
                        }
                }
        }
}


void MCI::findMRT2Step()
{
    int j;
    //constants
    const int MIN_STAT=100;  //minimum statistic: number of M(RT)^2 steps done before deciding if the step must be increased or decreased
    const int MIN_CONS=10;  //minimum consecutive: minimum number of consecutive loops without need of changing mrt2step
    const double TOLERANCE=0.05;  //tolerance: tolerance for the acceptance rate
    const int MAX_NUM_ATTEMPTS=100;  //maximum number of attempts: maximum number of time that the main loop can be executed
    const double SMALLEST_ACCEPTABLE_DOUBLE=1.e-50;

    //initialize index
    int cons_count = 0;  //number of consecutive loops without need of changing mrt2step
    int counter = 0;  //counter of loops
    double fact;
    while ( cons_count < MIN_CONS )
        {
            counter++;
            //do MIN_STAT M(RT)^2 steps
            this->sample(MIN_STAT,false);
            //increase or decrease mrt2step depending on the acceptance rate
            double rate = this->getAcceptanceRate();
            if ( rate > _targetaccrate+TOLERANCE )
                {
                    //need to increase mrt2step
                    cons_count=0;
                    fact = std::min(2.,rate/_targetaccrate);
                    for (j=0; j<_ndim; ++j)
                        {
                            *(_mrt2step+j) *= fact;
                        }
                } else if ( rate < _targetaccrate-TOLERANCE)
                {
                    //need to decrease mrt2step
                    cons_count=0;
                    fact = std::max(0.5, rate/_targetaccrate);
                    for (j=0; j<_ndim; ++j)
                        {
                            *(_mrt2step+j) *= fact;
                        }
                } else
                {
                    //mrt2step was ok
                    cons_count++;
                }
            //check if the loop is running since too long
            if ( counter > MAX_NUM_ATTEMPTS ) cons_count=MIN_CONS;
            //mrt2step = Infinity
            for (j=0; j<_ndim; ++j)
                {
                    if ( *(_mrt2step+j) - ( *(*(_irange+j)+1) - *(*(_irange+j)) ) > 0. )
                        {
                            cons_count = MIN_CONS; //make the main loop terminate
                            *(_mrt2step+j) = ( *(*(_irange+j)+1) - *(*(_irange+j)) );
                        }
                }
            //mrt2step ~ 0
            for (j=0; j<_ndim; ++j)
                {
                    if ( *(_mrt2step+j) < SMALLEST_ACCEPTABLE_DOUBLE )
                        {
                            cons_count = MIN_CONS; //make the main loop terminate
                            *(_mrt2step+j) = SMALLEST_ACCEPTABLE_DOUBLE;
                        }
                }
        }
}


void MCI::doStepMRT2(bool * flagacc)
{
    // propose a new position x
    this->computeNewX();

    // find the corresponding sampling function value
    this->computeNewSamplingFunction();

    //determine if the proposed x is accepted or not
    *flagacc = false;
    if ( this->computeAcceptance() > _rd(_rgen) ) *flagacc=true;

    //update some values according to the acceptance of the mrt2 step
    if ( *flagacc )
        {
            //accepted
            _acc++;
            //update the walker position x
            this->updateX();
            //update the sampling function values pdfx
            this->updateSamplingFunction();
            //if there are some call back functions, invoke them
            for (MCICallBackOnAcceptanceInterface * cback : _cback){
                cback->callBackFunction(_xold, _flagMC);
            }
        } else
        {
            //rejected
            _rej++;
        }
}


void MCI::applyPBC(double * v)
{
    for (int i=0; i<_ndim; ++i)
        {
            while ( *(v+i) < *(*(_irange+i)) )
                {
                    *(v+i) += *(*(_irange+i)+1) - *(*(_irange+i)) ;
                }
            while ( *(v+i) > *(*(_irange+i)+1) )
                {
                    *(v+i) -= ( *(*(_irange+i)+1) - *(*(_irange+i)) ) ;
                }
        }
}


void MCI::updateX()
{
    double * foo;
    foo = _xold;
    _xold = _xnew;
    _xnew = foo;
}


void MCI::newRandomX()
{
    //generate a new random x (within the irange)
    for (int i=0; i<_ndim; ++i)
        {
            *(_xold+i) = *(*(_irange+i)) + ( *(*(_irange+i)+1) - *(*(_irange+i)) ) * _rd(_rgen);
        }
    //move automatically accepted
    _acc++;
}


void MCI::resetAccRejCounters()
{
    _acc = 0;
    _rej = 0;
}


void MCI::computeNewX()
{
    for (int i=0; i<_ndim; ++i)
        {
            *(_xnew+i) = *(_xold+i) + (*(_mrt2step+i)) * (_rd(_rgen)-0.5);
        }
    applyPBC(_xnew);
}


void MCI::updateSamplingFunction()
{
    for (std::vector<MCISamplingFunctionInterface>::size_type i=0; i<_pdf.size(); ++i)
        {
            _pdf[i]->newToOld();
        }
}


double MCI::computeAcceptance()
{
    double acceptance=1.;
    for (std::vector<MCISamplingFunctionInterface>::size_type i=0; i<_pdf.size(); ++i)
        {
            acceptance*=_pdf[i]->getAcceptance();
        }
    return acceptance;
}


void MCI::computeOldSamplingFunction()
{
    for (std::vector<MCISamplingFunctionInterface>::size_type i=0; i<_pdf.size(); ++i)
        {
            _pdf[i]->computeNewSamplingFunction(_xold);
            _pdf[i]->newToOld();
        }
}


void MCI::computeNewSamplingFunction()
{
    for (std::vector<MCISamplingFunctionInterface>::size_type i=0; i<_pdf.size(); ++i)
        {
            _pdf[i]->computeNewSamplingFunction(_xnew);
        }
}


void MCI::computeObservables()
{
    for (std::vector<MCIObservableFunctionInterface>::size_type i=0; i<_obs.size(); ++i)
        {
            _obs[i]->computeObservables(_xold);
        }

}


void MCI::saveObservables()
{
    int idx=0;
    //save in _datax the observables contained in _obs
    for (std::vector<MCIObservableFunctionInterface>::size_type i=0; i<_obs.size(); ++i)
        {
            for (int j=0; j<_obs[i]->getNObs(); ++j)
                {
                    _datax[_ridx][idx]=_obs[i]->getObservable(j);
                    idx++;
                }
        }
}



//   --- Setters


void MCI::storeObservablesOnFile(const char * filepath, const int &freq)
{
    _pathobsfile.assign(filepath);
    _freqobsfile = freq;
    _flagobsfile=true;
}


void MCI::storeWalkerPositionsOnFile(const char * filepath, const int &freq)
{
    _pathwlkfile.assign(filepath);
    _freqwlkfile = freq;
    _flagwlkfile=true;
}


void MCI::clearCallBackOnAcceptance(){
    _cback.clear();
}


void MCI::addCallBackOnAcceptance(MCICallBackOnAcceptanceInterface * cback){
    _cback.push_back(cback);
}


void MCI::clearObservables()
{
    _obs.clear();
    _nobsdim=0;
}


void MCI::addObservable(MCIObservableFunctionInterface * obs)
{
    _obs.push_back(obs);
    _nobsdim+=obs->getNObs();
}


void MCI::clearSamplingFunctions()
{
    _pdf.clear();
    _flagpdf = false;
}


void MCI::addSamplingFunction(MCISamplingFunctionInterface * mcisf)
{
    _pdf.push_back(mcisf);
    _flagpdf = true;
}


void MCI::setTargetAcceptanceRate(const double * targetaccrate)
{
    _targetaccrate = *(targetaccrate);
}


void MCI::setMRT2Step(const double * mrt2step)
{
    for (int i=0; i<_ndim; ++i)
        {
            *(_mrt2step+i) = *(mrt2step+i);
        }
}


void MCI::setX(const double * x)
{
    for (int i=0; i<_ndim; ++i)
        {
            *(_xold+i) = *(x+i);
        }
}


void MCI::setIRange(const double * const * irange)
{
    // Set _irange and the initial walker position _x
    for (int i=0; i<_ndim; ++i)
        {
            *(*(_irange+i)) = *(*(irange+i));
            *(*(_irange+i)+1) = *(*(irange+i)+1);
            *(_xold+i) = ( *(*(_irange+i)+1) + *(*(_irange+i)) )*0.5;
        }
    // Set the integration volume
    _vol=1.;
    for (int i=0; i<_ndim; ++i)
        {
            _vol = _vol*( *(*(_irange+i)+1) - *(*(_irange+i)) );
        }
}


//   --- Constructor and Destructor

MCI::MCI(const int & ndim)
{
    int i;
    // _ndim
    _ndim = ndim;
    // _irange
    _irange = new double*[_ndim];
    for (i=0; i<_ndim; ++i)
        {
            *(_irange+i) = new double[2];
            *(*(_irange+i)) = -std::numeric_limits<double>::max()/(i+1);
            *(*(_irange+i)+1) = std::numeric_limits<double>::max()/(i+1);
        }
    // _vol
    _vol=0.;
    // _x
    _xold = new double[_ndim];
    for (i=0; i<_ndim; ++i)
        {
            _xold[i] = 0.;
        }
    _xnew = new double[_ndim];
    for (i=0; i<_ndim; ++i)
        {
            _xnew[i] = 0.;
        }
    _flagwlkfile=false;
    // _mrt2step
    _mrt2step = new double[_ndim];
    for (i=0; i<_ndim; ++i)
        {
            _mrt2step[i] = INITIAL_STEP;
        }
    // probability density function
    _flagpdf = false;
    // initialize info about observables
    _nobsdim=0;
    _flagobsfile=false;
    // initialize random generator
    _rgen = std::mt19937_64(_rdev());
    _rd = std::uniform_real_distribution<double>(0.,1.);
    // initialize all the other variables
    _targetaccrate=0.5;
    _flagMC=false;
}

MCI::~MCI()
{
    // _irange
    for (int i=0; i<_ndim; ++i){delete[] *(_irange+i);}
    delete[] _irange;
    // _xold and _xnew
    delete [] _xold;
    delete [] _xnew;
    // _mrt2step
    delete [] _mrt2step;
    // vectors
    _pdf.clear();
    _obs.clear();
}
