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
        const double INITIAL_STEP=0.05;

        std::random_device _rdev;
        std::mt19937_64 _rgen;
        std::uniform_real_distribution<double> _rd;  //after initialization (done in the constructor) can be used with _rd(_rgen)

        int _ndim;  // number of dimensions
        double * _lbound; // integration lower bounds
        double * _ubound; // integration upper bounds
        double _vol;  // Integration volume

        double * _xold;  // walker position
        double * _xnew;  //walker proposed position

        int _NfindMRT2Iterations; // how many MRT2 step adjustment iterations to do before integrating
        int _NdecorrelationSteps; // how many decorrelation steps to do before integrating
        double _targetaccrate;  // desired acceptance ratio

        // main object vectors/containers
        std::unique_ptr<TrialMoveInterface> _trialMove;

        SamplingFunctionContainer _pdfcont; // sampling function container
        ObservableContainer _obscont; // observable container used during integration

        std::vector< std::unique_ptr<CallBackOnMoveInterface> > _cbacks;  // Vector of acceptance callback functions

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

        // prepare new sampling run
        void initializeSampling(ObservableContainer * obsCont);

        void updateVolume();
        double getMinBoxLen() const; // min(ubound[i]-lbound[i])
        void applyPBC(double v[]) const;
        //void computeNewX();
        void acceptX();
        void rejectX();
        //use this if there is a pdf, returns how many x were changed (0 if not accepted)
        int doStepMRT2(int changedIdx[] /*unused, because all-particle steps*/); // returns number of changed positions, i.e. 0 or ndim

        // these are used before sampling
        void findMRT2Step();
        void initialDecorrelation();

        // sample without taking data
        void sample(int npoints);
        // fill data with samples
        void sample(int npoints, ObservableContainer &container);

        // call callbacks
        void callBackOnMove(const double x[], bool accepted);

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
        void setIRange(const double lbound[], const double ubound[]);

        void setX(const double x[]);
        void newRandomX();  // use if you want to take a new random _xold

        void setMRT2Step(double mrt2step); // set all identical
        void setMRT2Step(int i, double mrt2step); // set certain element
        void setMRT2Step(const double mrt2step[]); // set all elements
        void setTargetAcceptanceRate(double targetaccrate); // acceptance rate target used in findMRT2Step
        // how many MRT2 step adjustment iterations to do
        void setNfindMRT2Iterations(int niterations /* -1 == auto, 0 == disabled */){_NfindMRT2Iterations=niterations;}
        // how many decorrelation steps to do
        void setNdecorrelationSteps(int nsteps /* -1 == auto, 0 == disabled */){_NdecorrelationSteps=nsteps;}

        // --- Adding objects to MCI (everything you pass will be cloned!)

        //void setTrialMove(const TrialMoveInterface &tmove);
        //void setTrialMove(const TrialMoveInterface &tmove);

        void addObservable(const ObservableFunctionInterface &obs /* MCI adds accumulator and estimator for this obs, with following options: */,
                           int blocksize, /* if > 1, use fixed block size and assume uncorrelated samples, if <= 0, use no blocks and no error calculation */
                           int nskip, /* evaluate observable only every n-th step NOTE: now a block consists of $blocksize non-skipped samples */
                           bool flag_equil, /* observable wants to be equilibrated when using automatic initial decorrelation (blocksize must be > 0) */
                           bool flag_correlated /* should block averages be treated as correlated samples? (blocksize must be > 0) */
                           );
        void addObservable(const ObservableFunctionInterface &obs, int blocksize = 1, int nskip = 1) {
            addObservable(obs, blocksize, nskip, blocksize>0, blocksize==1); // safe&easy defaults, appropriate for most cases
        }
        void addObservable(const ObservableFunctionInterface &obs, int blocksize, int nskip, bool flag_equil, EstimatorType estimType /*enumerator, see Factories*/);
        void clearObservables(); // clear

        void addSamplingFunction(const SamplingFunctionInterface &mcisf);
        void clearSamplingFunctions();

        void addCallBackOnMove(const CallBackOnMoveInterface &cback);
        void clearCallBackOnMove();

        void storeObservablesOnFile(const std::string &filepath, int freq);
        void storeWalkerPositionsOnFile(const std::string &filepath, int freq);

        // --- Getters
        int getNDim() const { return _ndim; }
        double getLBound(int i) const { return _lbound[i]; }
        double getUBound(int i) const { return _ubound[i]; }

        double getX(int i) const { return _xold[i];}
        const double * getX() const { return _xold; }
        double getMRT2Step(int i) const { return (i < _trialMove->getNStepSizes()) ? _trialMove->getStepSize(i) : 0.; } // this is easy to get wrong, so we make it safer
        int getNfindMRT2Iterations() const { return _NfindMRT2Iterations; }
        int getNdecorrelationSteps() const { return _NdecorrelationSteps; }

        const ObservableFunctionInterface & getObservable(int i) const { return _obscont.getObservableFunction(i); }
        int getNObs() const { return _obscont.getNObs(); }
        int getNObsDim() const { return _obscont.getNObsDim(); }

        const SamplingFunctionInterface & getSamplingFunction(int i) const { return _pdfcont.getSamplingFunction(i); }
        int getNPDF() const { return _pdfcont.getNPDF(); }

        const CallBackOnMoveInterface & getCallBackOnMove(int i) const { return *_cbacks[i]; }
        int getNCallBacks() const { return _cbacks.size(); }

        double getTargetAcceptanceRate() const { return _targetaccrate; }
        double getAcceptanceRate() const { return (_acc>0) ? static_cast<double>(_acc)/(static_cast<double>(_acc)+_rej) : 0.; }

        // --- Integrate

        // Actual integrate implemention. With flags to skip the configured step adjustment/decorrelation.
        void integrate(int Nmc, double average[], double error[], bool doFindMRT2step = true, bool doDecorrelation = true);
    };
}  // namespace mci

#endif
