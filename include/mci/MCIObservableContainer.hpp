#ifndef MCI_MCIOBSERVABLECONTAINER_HPP
#define MCI_MCIOBSERVABLECONTAINER_HPP

#include "mci/MCIAccumulatorInterface.hpp"
#include "mci/MCIObservableFunctionInterface.hpp"

#include <fstream>
#include <functional>
#include <vector>

class MCIObservableContainer
{
private:
    int _nobsdim {0}; // stores total dimension of contained observables

    // Accumulators
    std::vector< MCIAccumulatorInterface * > _accus;

    // Estimator functions used to obtain result of MC integration
    std::vector< std::function< void (double [] /*avg*/, double [] /*error*/) > > _estims; // corresponding accumulators are already bound

    // flags
    std::vector< unsigned char /*avoid special vector<bool>..*/ > _flags_equil; // equilibrate this observable when using automatic decorrelation?

public:
    explicit MCIObservableContainer() = default;
    ~MCIObservableContainer(){ this->clear(); }

    // simple getters
    int getNObs(){ return _accus.size(); }
    int getNObsDim(){ return _nobsdim; }

    MCIObservableFunctionInterface * getObservableFunction(int i){ return _accus[i]->getObservableFunction(); }
    bool getFlagEquil(int i){return (_flags_equil[i]>0);}

    // operational methods
    // add accumulator&estimator for an observable
    void addObservable(MCIAccumulatorInterface * accumulator,
                        const std::function< void (int /*nstored*/, int /*nobs*/, const double [] /*data*/, double [] /*avg*/, double [] /*error*/) > &estimator,
                        bool needsEquil);

    void allocate(int Nmc); // allocate data memory
    void accumulate(const double x[], bool flagacc); // process accumulation for position x (which is new if flagacc)
    void printObsValues(std::ofstream &file); // write last observables values to filestream
    void finalize(); // used after sampling to apply all necessary data normalization
    void estimate(double average[], double error[]); // eval estimators on finalized data and return average/error
    void reset(); // obtain clean state, but keep allocation
    void deallocate(); // free data memory
    void clear(); // clear everything
};


#endif
