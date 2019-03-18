#ifndef MCI_MCINTEGRATOR_HPP
#define MCI_MCINTEGRATOR_HPP

#include "mci/AccumulatorInterface.hpp"
#include "mci/CallBackOnMoveInterface.hpp"
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
    class MCI
    {
    protected:
        // Random
        std::random_device _rdev;
        std::mt19937_64 _rgen;
        std::uniform_real_distribution<double> _rd;  //after initialization (done in the constructor) can be used with _rd(_rgen)

        // Integration domain
        const int _ndim;  // number of dimensions
        double * _lbound; // integration lower bounds
        double * _ubound; // integration upper bounds
        double _vol;  // Integration volume

        // Main objects/vectors/containers
        WalkerState _wlkstate; // holds the current walker state (xold/xnew), including move information
        std::unique_ptr<TrialMoveInterface> _trialMove; // holds the object to perform walker moves
        SamplingFunctionContainer _pdfcont; // sampling function container

        ObservableContainer _obscont; // observable container used during integration
        std::vector< std::unique_ptr<CallBackOnMoveInterface> > _cbacks;  // Vector of callback-on-move functions

        // Settings
        int _NfindMRT2Iterations; // how many MRT2 step adjustment iterations to do before integrating
        int64_t _NdecorrelationSteps; // how many decorrelation steps to do before integrating
        double _targetaccrate; // desired acceptance ratio


        // File-I/O parameters

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
        // or are required to be of a different integer type for external reasons.
        int64_t _acc, _rej; // internal counters
        int64_t _ridx; // running index, which keeps track of the number of MC steps


        // --- Internal methods

        // these are used before sampling
        void findMRT2Step();
        void initialDecorrelation();

        // prepare new sampling run
        void initializeSampling(ObservableContainer * obsCont /*optional*/);

        // if there is a pdf, performs move and decides acc/rej
        void doStepMRT2();
        // else we use this to sample randomly (mostly for testing/examples)
        void doStepRandom();

        // sample without taking data
        void sample(int64_t npoints);
        // fill data with samples and do things like file output, if flagMC (i.e. main sampling)
        void sample(int64_t npoints, ObservableContainer &container, bool flagMC);

        // call callbacks
        void callBackOnMove();

        // store to file
        void storeObservables();
        void storeWalkerPositions();

        void checkIRange() const; // throws if internal ranges are invalid
        void updateVolume();
        double getMinBoxLen() const; // min(ubound[i]-lbound[i])
        void applyPBC(double v[]) const; // apply PBC to passed vector (len dim)
        void applyPBCUpdate(); // apply elementary PBC update to xnew of internal walker state

    public:
        explicit MCI(int ndim);  //Constructor, need the number of dimensions
        ~MCI();  //Destructor

        // --- Setters

        void setSeed(uint_fast64_t seed); // seed internal random number generator

        // keep walkers within these bounds during integration (defaults to full range of double floats)
        void setIRange(double lbound, double ubound); // set the same range on all dimensions
        void setIRange(const double lbound[], const double ubound[]);

        void setX(int i, const double val) { _wlkstate.xold[i] = val; }
        void setX(const double x[]);
        void newRandomX();  // use if you want to take a new random _xold

        void setMRT2Step(double mrt2step); // set all identical
        void setMRT2Step(int i, double mrt2step); // set certain element
        void setMRT2Step(const double mrt2step[]); // set all elements

        void setTargetAcceptanceRate(double targetaccrate); // acceptance rate target used in findMRT2Step
        // how many MRT2 step adjustment iterations to do
        void setNfindMRT2Iterations(int niterations /* -1 == auto, 0 == disabled */){_NfindMRT2Iterations=niterations;}
        // how many decorrelation steps to do
        void setNdecorrelationSteps(int64_t nsteps /* -1 == auto, 0 == disabled */){_NdecorrelationSteps=nsteps;}


        // --- Adding objects to MCI
        // Note: Objects passed by raw-ref will be cloned by MCI

        // Trial Moves
        void setTrialMove(const TrialMoveInterface &tmove); // pass an existing move to be cloned by MCI
        void setTrialMove(MoveType move /*enum, see Factories.hpp*/); // set trial move to default version of chosen builtin move
        void setTrialMove(SRRDType srrd /*enum*/, // set builtin SRRD-class move, with distribution srrd
                          int veclen = 0, /*0 means all-move, > 0 means single-vector move*/
                          int ntypes = 1, int typeEnds[] = nullptr /* see TypedTrialMove.hpp */
                          );

        // Observables
        void addObservable(const ObservableFunctionInterface &obs /* MCI adds accumulator and estimator for this obs, with following options: */,
                           int blocksize, /* if > 1, use fixed block size and assume uncorrelated samples, if 0, use no blocks and no error calculation */
                           int nskip, /* evaluate observable only every n-th step NOTE: now one block is used for blocksize*nskip steps */
                           bool flag_equil, /* observable wants to be equilibrated when using automatic initial decorrelation (blocksize must be > 0) */
                           bool flag_correlated /* should block averages be treated as correlated samples? (blocksize must be > 0) */
                           );
        void addObservable(const ObservableFunctionInterface &obs, int blocksize = 1, int nskip = 1) {
            addObservable(obs, blocksize, nskip, blocksize>0, blocksize==1); // safe&easy defaults, appropriate for most cases
        }
        void addObservable(const ObservableFunctionInterface &obs, int blocksize, int nskip, bool flag_equil, EstimatorType estimType /*enum, see Factories-hpp*/);
        void clearObservables(); // clear

        // Sampling Functions
        void addSamplingFunction(const SamplingFunctionInterface &mcisf);
        void clearSamplingFunctions();

        // Callbacks
        void addCallBack(const CallBackOnMoveInterface &cback);
        void clearCallBacks();

        // enable file printout to given files, with frequency freq
        void storeObservablesOnFile(const std::string &filepath, int freq);
        void storeWalkerPositionsOnFile(const std::string &filepath, int freq);


        // --- Getters

        int getNDim() const { return _ndim; }
        double getLBound(int i) const { return _lbound[i]; }
        double getUBound(int i) const { return _ubound[i]; }

        double getX(int i) const { return _wlkstate.xold[i];}
        const double * getX() const { return _wlkstate.xold; }

        double getMRT2Step(int i) const { return (i < _trialMove->getNStepSizes()) ? _trialMove->getStepSize(i) : 0.; } // this is easy to get wrong, so we make it safer
        double getTargetAcceptanceRate() const { return _targetaccrate; }
        double getAcceptanceRate() const { return (_acc>0) ? static_cast<double>(_acc)/(static_cast<double>(_acc)+_rej) : 0.; }

        int getNfindMRT2Iterations() const { return _NfindMRT2Iterations; }
        int getNdecorrelationSteps() const { return _NdecorrelationSteps; }

        const TrialMoveInterface & getTrialMove() const { return *_trialMove; }

        const SamplingFunctionInterface & getSamplingFunction(int i) const { return _pdfcont.getSamplingFunction(i); }
        int getNPDF() const { return _pdfcont.getNPDF(); }

        const ObservableFunctionInterface & getObservable(int i) const { return _obscont.getObservableFunction(i); }
        int getNObs() const { return _obscont.getNObs(); }
        int getNObsDim() const { return _obscont.getNObsDim(); }

        const CallBackOnMoveInterface & getCallBackOnMove(int i) const { return *_cbacks[i]; }
        int getNCallBacks() const { return _cbacks.size(); }

        // --- Integrate

        // Actual integrate implemention. With flags to skip the configured step adjustment/decorrelation.
        void integrate(int64_t Nmc, double average[], double error[], bool doFindMRT2step = true, bool doDecorrelation = true);
    };
}  // namespace mci

#endif
