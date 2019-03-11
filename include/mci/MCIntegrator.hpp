#ifndef MCI_MCINTEGRATOR_HPP
#define MCI_MCINTEGRATOR_HPP

#include "mci/AccumulatorInterface.hpp"
#include "mci/CallBackOnMoveInterface.hpp"
#include "mci/ObservableContainer.hpp"
#include "mci/ObservableFunctionInterface.hpp"
#include "mci/SamplingFunctionInterface.hpp"

#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <memory>

namespace mci
{
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

        int _NfindMRT2Iterations; // how many MRT2 step adjustment iterations to do before integrating
        int _NdecorrelationSteps; // how many decorrelation steps to do before integrating
        double _targetaccrate;  // desired acceptance ratio

        // main object vectors/containers
        std::vector< std::unique_ptr<SamplingFunctionInterface> > _pdfs; //vector of sampling functions
        bool _flagpdf;  // did the user provide a sampling function?

        ObservableContainer _obscont; // observable container

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

        void computeNewSamplingFunction(); //compute the new sampling function with new coordinates
        void computeOldSamplingFunction(); //compute the new sampling function with the old coordinates
        // and it stores it in the old sampling
        void updateSamplingFunction(); // swap old and new sampling function values
        double computeAcceptance() const; //compute the acceptance

        void resetAccRejCounters();

        void updateVolume();
        void applyPBC(double v[]) const;
        void computeNewX();
        void updateX();
        //use this if there is a pdf, returns whether step was accepted or not
        bool doStepMRT2(bool * flags_xchanged /* tells which x are changed */);

        // these are used before sampling
        void findMRT2Step();
        void initialDecorrelation();

        // sample without taking data
        void sample(int npoints);
        // fill data with samples
        void sample(int npoints, ObservableContainer &container);

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

        void setMRT2Step(const double mrt2step[]);
        void setTargetAcceptanceRate(double targetaccrate); // acceptance rate target used in findMRT2Step
        // how many MRT2 step adjustment iterations to do
        void setNfindMRT2Iterations(int niterations /* -1 == auto, 0 == disabled */){_NfindMRT2Iterations=niterations;}
        // how many decorrelation steps to do
        void setNdecorrelationSteps(int nsteps /* -1 == auto, 0 == disabled */){_NdecorrelationSteps=nsteps;}


        void addObservable(const ObservableFunctionInterface & obs /* MCI adds accumulator and estimator for this obs, with following options: */,
                           int blocksize, /* if > 1, use fixed block size and assume uncorrelated samples, if <= 0, use no blocks and no error calculation */
                           int nskip, /* evaluate observable only every n-th step NOTE: now a block consists of $blocksize non-skipped samples */
                           bool flag_equil, /* observable wants to be equilibrated when using automatic initial decorrelation (blocksize must be > 0) */
                           bool flag_correlated /* should block averages be treated as correlated samples? (blocksize must be > 0) */
                           );
        void addObservable(const ObservableFunctionInterface & obs, int blocksize = 1, int nskip = 1) {
            addObservable(obs, blocksize, nskip, blocksize>0, blocksize==1); // safe&easy defaults, appropriate for most cases
        }
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
        double getMRT2Step(int i) const { return _mrt2step[i]; }
        const double * getMRT2Step() const { return _mrt2step; }
        int getNfindMRT2Iterations() const { return _NfindMRT2Iterations; }
        int getNdecorrelationSteps() const { return _NdecorrelationSteps; }

        const ObservableFunctionInterface & getObservable(int i) const { return _obscont.getObservableFunction(i); }
        int getNObs() const { return _obscont.getNObs(); }
        int getNObsDim() const { return _obscont.getNObsDim(); }

        const SamplingFunctionInterface & getSamplingFunction(int i) const { return *_pdfs[i]; }
        int getNSampF() const { return _pdfs.size(); }

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
