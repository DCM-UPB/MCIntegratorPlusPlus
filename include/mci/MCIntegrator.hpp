#ifndef MCI_MCINTEGRATOR_HPP
#define MCI_MCINTEGRATOR_HPP

#include "mci/AccumulatorInterface.hpp"
#include "mci/DomainInterface.hpp"
#include "mci/Factories.hpp"
#include "mci/ObservableContainer.hpp"
#include "mci/ObservableFunctionInterface.hpp"
#include "mci/SamplingFunctionContainer.hpp"
#include "mci/SamplingFunctionInterface.hpp"
#include "mci/TrialMoveInterface.hpp"
#include "mci/WalkerState.hpp"


#include <cstdint>
#include <fstream>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace mci
{
// Main class of the MCI library
// Holds all objects and provides the user interface for MC integration.
class MCI
{
private:
    const int _ndim;  // number of dimensions

    // Random
    std::random_device _rdev;
    std::mt19937_64 _rgen;
    std::uniform_real_distribution<double> _rd; // used to decide on acceptance (and for full random moves)

    // Main objects/vectors/containers
    WalkerState _wlkstate; // holds the current walker state (xold/xnew), including move information
    std::unique_ptr<DomainInterface> _domain; // holds the integration domain (init: unbound)
    std::unique_ptr<TrialMoveInterface> _trialMove; // holds the object to perform walker moves (init: uniform all-move)
    SamplingFunctionContainer _pdfcont; // sampling function container (init: empty)
    ObservableContainer _obscont; // observable container used during integration (init: empty)
    std::function<void(const MCI &)> _cback{}; // callback function (see setCallback() below)

    // Settings
    int _NfindMRT2Iterations; // how many MRT2 step adjustment iterations to do before integrating
    int64_t _NdecorrelationSteps; // how many decorrelation steps to do before integrating
    double _targetaccrate; // desired acceptance ratio

    // File-I/O parameters:
    // observables
    std::ofstream _obsfile; //ofstream for storing obs values while sampling
    std::string _pathobsfile;
    int _freqobsfile{};
    bool _flagobsfile; // should write an output file with sampled obs values?
    // walkers
    std::ofstream _wlkfile; //ofstream for storing obs values while sampling
    std::string _pathwlkfile;
    int _freqwlkfile{};
    bool _flagwlkfile; // should write an output file with sampled obs values?

    // internal counters
    // NOTE: All integers are int, except if they are directly counting MC steps (int64_t then)
    // or are required to be of a different integer type for other reasons (e.g. seed).
    int64_t _acc, _rej; // internal counters
    int64_t _ridx; // running index, which keeps track of the number of MC steps


    // --- Internal methods

    // these are used before sampling
    void findMRT2Step();
    void initialDecorrelation();

    // prepare new sampling run
    void initializeSampling(ObservableContainer * obsCont /*optional*/);

    // if there is a pdf, performs move and decides acc/rej
    void doStepMRT2(bool callbackPDF);
    // else we use this to sample randomly (mostly for testing/examples)
    void doStepRandom();

    // sample without taking data
    void sample(int64_t npoints);
    // fill data with samples and do things like file output, if flagMC (i.e. main sampling)
    void sample(int64_t npoints, ObservableContainer &container, bool flagMC);


    // store to file
    void storeObservables();
    void storeWalkerPositions();

public:
    explicit MCI(int ndim);  //Constructor, need the number of dimensions
    ~MCI() = default;  // Destructor (empty)

    // --- Setters

    void setSeed(uint_fast64_t seed); // seed internal random number generator

    // - manipulate initial position
    void setX(int i, double val);
    void setX(const double x[]);
    void moveX(); // use if you want to do a single move manually (uses configured trial move)
    void newRandomX(); // use to set a new random x within the configured domain
    void centerX() { _domain->getCenter(_wlkstate.xold); } // reset x back to domain center

    // - manipulate move stepsizes
    void setMRT2Step(double mrt2step); // set all identical
    void setMRT2Step(int i, double mrt2step); // set certain element
    void setMRT2Step(const double mrt2step[]); // set all elements

    // - manipulate automatic routines
    void setTargetAcceptanceRate(double targetaccrate); // acceptance rate target used in findMRT2Step
    // how many MRT2 step adjustment iterations to do
    void setNfindMRT2Iterations(int niterations /*N<0 -> auto with max abs(N) iterations, 0 -> off, N>0 -> fixed N iterations*/)
    {
        _NfindMRT2Iterations = niterations;
    }
    // how many decorrelation steps to do
    void setNdecorrelationSteps(int64_t nsteps /*N<0 -> auto with max abs(N) MC steps, 0 -> off, N>0 -> fixed N MC steps */)
    {
        _NdecorrelationSteps = nsteps;
    }


    // --- Adding objects to MCI
    // Note: Objects passed by raw-ref will be cloned by MCI

    // Domain Setters (previously set domain will be returned or destroyed if not taken)
    std::unique_ptr<DomainInterface> setDomain(std::unique_ptr<DomainInterface> domain); // move a domain to be owned by MCI
    std::unique_ptr<DomainInterface> setDomain(const DomainInterface &domain) { return this->setDomain(domain.clone()); } // pass a domain to be cloned by MCI
    std::unique_ptr<DomainInterface> resetDomain(); // reset the domain to unbound

    // keep walkers within these bounds during integration (using periodic boundaries)
    // NOTE: If you use these, any prior domain will be replaced with OrthoPeriodicDomain!!
    void setIRange(double lbound, double ubound); // set the same range on all dimensions
    void setIRange(const double lbounds[], const double ubounds[]);

    // Trial Move Setters (previously set moves will be returned or destroyed if not taken)
    std::unique_ptr<TrialMoveInterface> setTrialMove(std::unique_ptr<TrialMoveInterface> tmove); // move a move to be owned by MCI
    std::unique_ptr<TrialMoveInterface> setTrialMove(const TrialMoveInterface &tmove) { return this->setTrialMove(tmove.clone()); } // pass a move to be cloned by MCI
    std::unique_ptr<TrialMoveInterface> setTrialMove(MoveType move /*enum, see Factories.hpp*/); // set trial move to default version of chosen builtin move
    std::unique_ptr<TrialMoveInterface> setTrialMove(SRRDType srrd /*enum*/, // set builtin SRRD-class move, with distribution srrd
                                                     int veclen = 0, /*0 means all-move, > 0 means single-vector move*/
                                                     int ntypes = 1, int typeEnds[] = nullptr /* see TypedTrialMove.hpp */
                                                    );


    // Observables
    void addObservable(std::unique_ptr<ObservableFunctionInterface> obs /* MCI adds accumulator and estimator for this obs, with following options: */,
                       int blocksize, /* if > 1, use fixed block size, if 0, use no blocks and no error calculation */
                       int nskip, /* evaluate observable only every n-th step. NOTE: now one block is used for blocksize*nskip steps */
                       bool flag_equil, /* observable wants to be equilibrated when using automatic initial decorrelation (blocksize must be > 0) */
                       bool flag_correlated /* should block averages be treated as correlated samples? (blocksize must be > 0) */
                      );
    void addObservable(const ObservableFunctionInterface &obs, int blocksize, int nskip, bool flag_equil, bool flag_correlated)
    {
        this->addObservable(obs.clone(), blocksize, nskip, flag_equil, flag_correlated);
    }

    void addObservable(std::unique_ptr<ObservableFunctionInterface> obs, int blocksize = 1, int nskip = 1)
    {
        this->addObservable(std::move(obs), blocksize, nskip, blocksize > 0, blocksize == 1); // safe&easy defaults, appropriate for most cases
    }
    void addObservable(const ObservableFunctionInterface &obs, int blocksize = 1, int nskip = 1)
    {
        this->addObservable(obs.clone(), blocksize, nskip);
    }

    void addObservable(std::unique_ptr<ObservableFunctionInterface> obs, int blocksize, int nskip, bool flag_equil, EstimatorType estimType /*enum, see Factories-hpp*/);
    void addObservable(const ObservableFunctionInterface &obs, int blocksize, int nskip, bool flag_equil, EstimatorType estimType)
    {
        this->addObservable(obs.clone(), blocksize, nskip, flag_equil, estimType);
    }

    std::unique_ptr<ObservableFunctionInterface> popObservable(); // remove last observable (returns it for you to optionally take it back)
    void clearObservables() { _obscont.clear(); } // delete all observables


    // Sampling Functions
    void addSamplingFunction(std::unique_ptr<SamplingFunctionInterface> pdf);
    void addSamplingFunction(const SamplingFunctionInterface &pdf) { this->addSamplingFunction(pdf.clone()); }
    std::unique_ptr<SamplingFunctionInterface> popSamplingFunction() { return _pdfcont.pop_back(); } // remove last pdf (returns it for you to optionally take it back)
    void clearSamplingFunctions() { _pdfcont.clear(); } // delete all pdfs

    // Callback Function
    // Set a callback function which may read(!) const MCI after every move and do something with the data.
    // This should not be abused to somehow add MCI control logic via captured references to objects contained in MCI.
    void setCallback(const std::function<void(const MCI &)> &cback) { _cback = cback; }
    void clearCallback() { _cback = nullptr; } // set empty callback

    // enable file printout to given files, with frequency freq
    void storeObservablesOnFile(const std::string &filepath, int freq);
    void clearObservableFile();
    void storeWalkerPositionsOnFile(const std::string &filepath, int freq);
    void clearWalkerFile();

    // --- Getters

    int getNDim() const { return _ndim; }
    double getX(int i) const { return _wlkstate.xold[i]; }
    const double * getX() const { return _wlkstate.xold; }

    double getMRT2Step(int i) const;
    double getTargetAcceptanceRate() const { return _targetaccrate; }
    double getAcceptanceRate() const;

    int getNfindMRT2Iterations() const { return _NfindMRT2Iterations; }
    int64_t getNdecorrelationSteps() const { return _NdecorrelationSteps; }

    const DomainInterface &getDomain() const { return *_domain; }
    TrialMoveInterface &getTrialMove() const { return *_trialMove; }

    SamplingFunctionInterface &getSamplingFunction(int i) const { return _pdfcont.getSamplingFunction(i); }
    int getNPDF() const { return _pdfcont.getNPDF(); }

    ObservableFunctionInterface &getObservable(int i) const { return _obscont.getObservableFunction(i); }
    int getNObs() const { return _obscont.getNObs(); }
    int getNObsDim() const { return _obscont.getNObsDim(); }


    // --- Integrate

    // Actual integrate implemention. With flags to skip the configured step adjustment/decorrelation.
    void integrate(int64_t Nmc, double average[], double error[], bool doFindMRT2step = true, bool doDecorrelation = true);
};
}  // namespace mci

#endif
