#ifndef MCI_MCINTEGRATOR_HPP
#define MCI_MCINTEGRATOR_HPP

#include "mci/MCICallBackOnAcceptanceInterface.hpp"
#include "mci/MCIObservableFunctionInterface.hpp"
#include "mci/MCISamplingFunctionInterface.hpp"
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <tuple>


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
    int _nblocks; // how many blocks to use for error estimation (0 for auto -> high RAM usage)

    double _targetaccrate;  // desired acceptance ratio

    std::vector<MCISamplingFunctionInterface *> _pdf; //vector of sampling functions
    bool _flagpdf;  // did the user provide a sampling function?

    std::vector<MCIObservableFunctionInterface *> _obs;  // Vector of observable functions
    int _nobsdim;

    std::vector<MCICallBackOnAcceptanceInterface *> _cback;  // Vector of observable functions

    int _acc, _rej;  // internal counters
    int _ridx;  // running index, which keeps track of the number of MC steps
    int _bidx; // index of the current block/datax element
    double * _datax;  // array that will contain all the measured observable (or block averages if used)
    bool _flagMC; //flag that is true only when MCI is accumulating data for the integral

    std::ofstream _obsfile;  //ofstream for storing obs values while sampling
    std::string _pathobsfile;
    int _freqobsfile{};
    bool _flagobsfile;  // should write an output file with sampled obs values?

    std::ofstream _wlkfile;  //ofstream for storing obs values while sampling
    std::string _pathwlkfile;
    int _freqwlkfile{};
    bool _flagwlkfile;  // should write an output file with sampled obs values?


    // --- Internal methods
    void computeObservables(); //compute observables with old coordinates
    void saveObservables(); //save the observables into _datax

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

    void findMRT2Step();
    void initialDecorrelation();

    void sample(const long &npoints, const bool &flagobs, const long &stepsPerBlock = 1);

    void storeObservables();
    void storeWalkerPositions();

public:
    explicit MCI(const int & ndim);  //Constructor, need the number of dimensions
    ~MCI();  //Destructor

    // --- Setters
    void setSeed(uint_fast64_t seed);

    // keep walkers within these bounds during integration (defaults to full range of double floats)
    void setIRange(const double &lbound, const double &ubound); // set the same range on all dimensions
    void setIRange(const double * lbound, const double * ubound);

    void setX(const double * x);
    void newRandomX();  // use if you want to take a new random _xold

    void setMRT2Step(const double * mrt2step);
    void setNfindMRT2steps(const int niterations /* -1 == auto, 0 == disabled */){_NfindMRT2steps=niterations;} // how many MRT2 step adjustment iterations to do before integrating
    void setNdecorrelationSteps(const int nsteps /* -1 == auto, 0 == disabled */){_NdecorrelationSteps=nsteps;} // how many decorrelation steps to do before integrating
    void setNBlocks(const int nblocks /* 0 == auto -> high RAM usage */){_nblocks=nblocks;} // how many blocks to use for error estimation
    void setTargetAcceptanceRate(double targetaccrate);

    void addObservable(MCIObservableFunctionInterface * obs);
    void clearObservables();

    void addSamplingFunction(MCISamplingFunctionInterface * mcisf);
    void clearSamplingFunctions();

    void addCallBackOnAcceptance(MCICallBackOnAcceptanceInterface * cback);
    void clearCallBackOnAcceptance();

    void storeObservablesOnFile(const char * filepath, const int &freq);
    void storeWalkerPositionsOnFile(const char * filepath, const int &freq);

    // --- Getters
    int getNDim(){return _ndim;}
    double getLBound(const int &i){return _lbound[i];}
    double getUBound(const int &i){return _ubound[i];}
    std::pair<double, double> getIRange(const int &i){return std::pair<double, double>(_lbound[i], _ubound[i]);}

    double getX(const int &i){return _xold[i];}
    double getMRT2Step(const int &i){return _mrt2step[i];}
    int getNfindMRT2steps(){return _NfindMRT2steps;}
    int getNdecorrelationSteps(){return _NdecorrelationSteps;}
    int getNBlocks(){return _nblocks;}

    MCIObservableFunctionInterface * getObservable(const int &i){return _obs[i];}
    int getNObs(){return _obs.size();}
    int getNObsDim(){return _nobsdim;}

    MCISamplingFunctionInterface * getSamplingFunction(const int &i){return _pdf[i];}
    int getNSampF(){return _pdf.size();}

    MCICallBackOnAcceptanceInterface * getCallBackOnAcceptance(const int &i){return _cback[i];}
    int getNCallBacks(){return _cback.size();}

    double getTargetAcceptanceRate(){return _targetaccrate;}
    double getAcceptanceRate(){return (_acc>0) ? static_cast<double>(_acc)/(_acc+_rej) : 0.;}

    // --- Integrate

    // Actual integrate implemention. With flags to skip the configured step adjustment/decorrelation.
    void integrate(const long &Nmc, double * average, double * error, bool doFindMRT2step = true, bool doDecorrelation = true);
};

#endif
