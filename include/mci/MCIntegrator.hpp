#ifndef MCI_MCINTEGRATOR_HPP
#define MCI_MCINTEGRATOR_HPP

#include "mci/MCICallBackOnAcceptanceInterface.hpp"
#include "mci/MCIObservableFunctionInterface.hpp"
#include "mci/MCISamplingFunctionInterface.hpp"
#include "mci/MCIAccumulatorInterface.hpp"
#include "mci/MCIObservableContainer.hpp"

#include <fstream>
#include <random>
#include <string>
#include <vector>

class MCI
{
protected:
    const double INITIAL_STEP=0.1;

    std::random_device _rdev;
    std::mt19937_64 _rgen;
    std::uniform_real_distribution<double> _rd;  //after initialization (done in the constructor) can be used with _rd(_rgen)

    int _ndim;  // number of dimensions
    double * _lbound; // integration lower bounds
    double * _ubound; // integration upper bounds
    double _vol;  // Integration volume

    double * _xold;  // walker position
    double * _xnew;  //walker proposed position
    double * _mrt2step;  // M(RT)^2 random step

    int _NfindMRT2steps; // how many MRT2 step adjustment iterations to do before integrating
    int _NdecorrelationSteps; // how many decorrelation steps to do before integrating

    double _targetaccrate;  // desired acceptance ratio

    // main object vectors
    std::vector<MCISamplingFunctionInterface *> _pdf; //vector of sampling functions
    bool _flagpdf;  // did the user provide a sampling function?

    std::vector< MCIObservableContainer > _obsc; // vector of observable containers (stored by value!)
    int _nobsdim; // total dimension of all observables

    std::vector<MCICallBackOnAcceptanceInterface *> _cback;  // Vector of acceptance callback functions

    // internal flags & counters
    int _acc, _rej;  // internal counters
    int _ridx;  // running index, which keeps track of the number of MC steps
    bool _flagMC; //flag that is true only when MCI is accumulating data for the integral

    // variables related to file IO
    std::ofstream _obsfile;  //ofstream for storing obs values while sampling
    std::string _pathobsfile;
    int _freqobsfile{};
    bool _flagobsfile;  // should write an output file with sampled obs values?

    std::ofstream _wlkfile;  //ofstream for storing obs values while sampling
    std::string _pathwlkfile;
    int _freqwlkfile{};
    bool _flagwlkfile;  // should write an output file with sampled obs values?


    // --- Internal methods
    void allocateObservables(int Nmc); // allocate data memory
    void accumulateObservables(bool flagacc); // process observable step (using _xold)
    void finalizeObservables(); // used after to apply all necessary data normalization
    void resetObservables(); // obtain clean state, but keep allocation
    void deallocateObservables(); // free data memory

    void computeNewSamplingFunction(); //compute the new sampling function with new coordinates
    void computeOldSamplingFunction(); //compute the new sampling function with the old coordinates
    // and it stores it in the old sampling
    void updateSamplingFunction(); //copy the new sampling function into the old one
    double computeAcceptance(); //compute the acceptance

    void resetAccRejCounters();

    void updateVolume();
    void applyPBC(double * v);
    void computeNewX();
    void updateX();
    bool doStepMRT2();  //use this if there is a pdf, returns whether step was accepted or not

    // these are used before sampling
    void findMRT2Step();
    void initialDecorrelation();

    // fill data with samples
    void sample(int npoints, bool flagobs);

    // store to file
    void storeObservables();
    void storeWalkerPositions();

public:
    explicit MCI(int ndim);  //Constructor, need the number of dimensions
    ~MCI();  //Destructor

    // --- Setters
    void setSeed(uint_fast64_t seed);

    // keep walkers within these bounds during integration (defaults to full range of double floats)
    void setIRange(double lbound, double ubound); // set the same range on all dimensions
    void setIRange(const double * lbound, const double * ubound);

    void setX(const double * x);
    void newRandomX();  // use if you want to take a new random _xold

    void setMRT2Step(const double * mrt2step);
    void setTargetAcceptanceRate(double targetaccrate); // acceptance rate target used in findMRT2Step
    // how many MRT2 step adjustment iterations to do
    void setNfindMRT2steps(int niterations /* -1 == auto, 0 == disabled */){_NfindMRT2steps=niterations;}
    // how many decorrelation steps to do
    void setNdecorrelationSteps(int nsteps /* -1 == auto, 0 == disabled */){_NdecorrelationSteps=nsteps;}


    void addObservable(MCIObservableFunctionInterface * obs /* MCI creates a container incl. accumulator for obs */,
                       bool flag_error = true /* need error on avg? */,
                       int nskip = 1 /* evaluate every n-th step */, int blocksize = 1 /* use fixed block size (if > 1)*/);
    void clearObservables(); // clear

    void addSamplingFunction(MCISamplingFunctionInterface * mcisf);
    void clearSamplingFunctions();

    void addCallBackOnAcceptance(MCICallBackOnAcceptanceInterface * cback);
    void clearCallBackOnAcceptance();

    void storeObservablesOnFile(const char * filepath, int freq);
    void storeWalkerPositionsOnFile(const char * filepath, int freq);

    // --- Getters
    int getNDim(){ return _ndim; }
    double getLBound(int i){ return _lbound[i]; }
    double getUBound(int i){ return _ubound[i]; }

    double getX(int i){ return _xold[i];}
    double getMRT2Step(int i){ return _mrt2step[i]; }
    int getNfindMRT2steps(){ return _NfindMRT2steps; }
    int getNdecorrelationSteps(){ return _NdecorrelationSteps; }

    MCIObservableFunctionInterface * getObservable(int i){ return _obsc[i].accu->getObservable(); }
    MCIAccumulatorInterface * getAccumulator(int i){ return _obsc[i].accu; }
    MCIObservableContainer & getObservableContainer(int i){ return _obsc[i]; }
    int getNObs(){ return _obsc.size(); }
    int getNObsDim(){ return _nobsdim; }

    MCISamplingFunctionInterface * getSamplingFunction(int i){ return _pdf[i]; }
    int getNSampF(){ return _pdf.size(); }

    MCICallBackOnAcceptanceInterface * getCallBackOnAcceptance(int i){ return _cback[i]; }
    int getNCallBacks(){ return _cback.size(); }

    double getTargetAcceptanceRate(){ return _targetaccrate; }
    double getAcceptanceRate(){ return (_acc>0) ? static_cast<double>(_acc)/(static_cast<double>(_acc)+_rej) : 0.; }

    // --- Integrate

    // Actual integrate implemention. With flags to skip the configured step adjustment/decorrelation.
    void integrate(int Nmc, double * average, double * error, bool doFindMRT2step = true, bool doDecorrelation = true);
};

#endif
