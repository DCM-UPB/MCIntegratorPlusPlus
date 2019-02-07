#ifndef MCINTEGRATOR
#define MCINTEGRATOR

#include "mci/MCISamplingFunctionInterface.hpp"
#include "mci/MCIObservableFunctionInterface.hpp"
#include "mci/MCICallBackOnAcceptanceInterface.hpp"
#include <vector>
#include <random>
#include <fstream>
#include <string>



class MCI
{
protected:
    const double INITIAL_STEP=0.1;

    std::random_device _rdev;
    std::mt19937_64 _rgen;
    std::uniform_real_distribution<double> _rd;  //after initialization (done in the constructor) can be used with _rd(_rgen)

    int _ndim;  // number of dimensions
    double ** _irange;  // integration ranges
    double _vol;  // Integration volume

    double * _xold;  // walker position
    double * _xnew;  //walker proposed position
    double * _mrt2step;  // M(RT)^2 random step

    double _targetaccrate;  // accepted and rejected moves
    int _acc, _rej;  // the MC integration will be done sampling from this pdf

    std::vector<MCISamplingFunctionInterface *> _pdf; //vector of sampling functions
    bool _flagpdf;  // did the user provide a sampling function?

    std::vector<MCIObservableFunctionInterface *> _obs;  // Vector of observable functions
    int _nobsdim;

    std::vector<MCICallBackOnAcceptanceInterface *> _cback;  // Vector of observable functions

    int _ridx;  // running index, which keeps track of the number of MC steps
    int _bidx; // index of the current block/datax element
    double ** _datax;  // array that will contain all the measured observable (or block averages if used)
    bool _flagMC; //flag that is true only when MCI is accumulating data for the integral

    std::ofstream _obsfile;  //ofstream for storing obs values while sampling
    std::string _pathobsfile;
    int _freqobsfile;
    bool _flagobsfile;  // should write an output file with sampled obs values?

    std::ofstream _wlkfile;  //ofstream for storing obs values while sampling
    std::string _pathwlkfile;
    int _freqwlkfile;
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

    void applyPBC(double * v);
    void computeNewX();
    void updateX();
    void doStepMRT2(bool * flagacc);  //use this if there is a pdf

    void findMRT2Step(const int &NfindMRT2stepIterations = -1);

    void sample(const long &npoints, const bool &flagobs, const long &stepsPerBlock = 1);

    void initialDecorrelation(const int &NdecorrelationSteps = -1);

    void storeObservables();
    void storeWalkerPositions();


public:
    MCI(const int & ndim);  //Constructor, need the number of dimensions
    ~MCI();  //Destructor

    // --- Setters
    void setSeed(const uint_fast64_t seed);

    void setIRange(const double * const * irange);

    void setX(const double * x);
    void newRandomX();  // use if you want to take a new random _xold
    void setMRT2Step(const double * mrt2step);

    void setTargetAcceptanceRate(const double * targetaccrate);

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
    double getIRange(const int &i, const int &j){return *(*(_irange+i)+j);}
    double getX(const int &i){return *(_xold+i);}
    double getMRT2Step(const int &i){return *(_mrt2step+i);}

    MCIObservableFunctionInterface * getObservable(const int &i){return _obs[i];}
    int getNObs(){return _obs.size();}
    int getNObsDim(){return _nobsdim;}

    MCISamplingFunctionInterface * getSamplingFunction(const int &i){return _pdf[i];}
    int getNSampF(){return _pdf.size();}

    MCICallBackOnAcceptanceInterface * getCallBackOnAcceptance(const int &i){return _cback[i];}
    int getNCallBacks(){return _cback.size();}

    double getTargetAcceptanceRate(){return _targetaccrate;}
    double getAcceptanceRate(){return (double(_acc)/(double(_acc+_rej)));}

    // --- Integrate

    // Wrapper for easiness/compatibility. Using automatic methods for findMRT2step and initial decorrelation (if enabled).
    void integrate(const long &Nmc, double * average, double * error, bool findMRT2step=true, bool initialdecorrelation = true, size_t nblocks = 0); // nblocks == 0 means using auto-blocking (no RAM benefit)

    // Actual integrate implemention. With integer controls: -1 -> auto, 0 -> disabled, >0 -> use fixed step/iteration counts, e.g. useful for parallel computation
    void integrate(const long &Nmc, double * average, double * error, int NfindMRT2stepIterations, int NdecorrelationSteps, size_t nblocks = 0);


};

#endif
