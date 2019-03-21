#include "mci/MCIntegrator.hpp"
#include "mci/SamplingFunctionInterface.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

#include "../common/TestMCIFunctions.hpp"

using namespace std;
using namespace mci;

int main(){
    const int NMC_pre = 10000; // used for the initial two tests
    const int NMC = 25000; // used for all-moves, multipled by 1/changeRate for non-all-moves
    const double CORRECT_RESULT = 0.5;

    Gauss pdf(3);
    XSquared obs1d;
    X2 obs3d(3); // effectively an updateable XYZSquared

    MCI mci(3);
    mci.setSeed(5649871);
    mci.addSamplingFunction(pdf);
    mci.setTargetAcceptanceRate(0.7); // easier for findMRT2Step()

    // we will later clear&add these again with different settings
    mci.addObservable(obs1d);
    mci.addObservable(obs3d);

    // the integral should provide 0.5 as answer!

    double x[3];
    x[0] = 5.; x[1] = -5.; x[2] = 10.;

    double average[4];
    double error[4];

    // Starting with default all-particle moves

    // this integral should give a wrong answer
    mci.setX(x);
    mci.integrate(NMC_pre, average, error, false, false);
    for (int i=0; i<mci.getNObsDim(); ++i) {
        //std::cout << "i " << i << ", average[i] " << average[i] << ", error[i] " << error[i] << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
        assert( fabs(average[i]-CORRECT_RESULT) > 2.*error[i] );
    }
    //std::cout << std::endl;

    // this integral, instead, will provide the right answer
    mci.setX(x);
    mci.integrate(NMC_pre, average, error, true, true);
    for (int i=0; i<mci.getNObsDim(); ++i) {
        //std::cout << "i " << i << ", average[i] " << average[i] << ", error[i] " << error[i] << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
        assert( fabs(average[i]-CORRECT_RESULT) < 2.*error[i] );
    }
    //std::cout << std::endl;


    // Now do the same using default versions of all builtin moves:
    for (auto mtype : list_all_MoveType) { // from Factories.hpp
        //std::cout << "Setting MoveType " << static_cast<size_t>(mtype) << std::endl;
        mci.setTrialMove(mtype); // set different trial move
        // set obs nskip
        const int nskip = floor(1./mci.getTrialMove().getChangeRate());
        //std::cout << "Setting nskip " << nskip << std::endl;
        mci.clearObservables();
        mci.addObservable(obs1d, 1, nskip);
        mci.addObservable(obs3d, 1, nskip);
        // integrate
        mci.integrate(NMC*nskip, average, error, true, false);
        for (int i=0; i<mci.getNObsDim(); ++i) {
            //std::cout << "i " << i << ", average[i] " << average[i] << ", error[i] " << error[i] << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
            assert( fabs(average[i]-CORRECT_RESULT) < 3.*error[i] ); // for this test and those below to pass safely, factor 2 is a bit small
        }
        //std::cout << std::endl;
    }


    // Now using all/vec moves with uniform distribution and multiple settings
    for (int veclen=0; veclen<4; veclen=1+2*veclen) { // vector length parameter
        for (int ntypes=1; ntypes<3; ++ntypes) {
            if (veclen>1 && ntypes > 1) { continue; }
            // set typeEnds
            int typeEnds[ntypes];
            if (ntypes==1) { typeEnds[0] = 3; }
            else { typeEnds[0] = 1; typeEnds[1] = 3;}

            // set move
            mci.setTrialMove(SRRDType::Uniform, veclen, ntypes, typeEnds);

            // set obs nskip
            const int nskip = floor(1./mci.getTrialMove().getChangeRate());
            //std::cout << "Setting nskip " << nskip << std::endl;
            mci.clearObservables();
            mci.addObservable(obs1d, 1, nskip);
            mci.addObservable(obs3d, 1, nskip);

            // integrate
            mci.integrate(NMC*nskip, average, error, true, false);
            for (int i=0; i<mci.getNObsDim(); ++i) {
                //std::cout << "i " << i << ", average[i] " << average[i] << ", error[i] " << error[i] << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
                assert( fabs(average[i]-CORRECT_RESULT) < 3.*error[i] ); // for all these tests to pass safely, factor 2 is a bit small
            }
            //std::cout << std::endl;
        }
    }

    // Now try to (just) set all/vec moves with all possible distributions
    // No need to run integration again, if we trust the standard library
    for (auto srrd : list_all_SRRDType) { // from Factories.hpp
        for (int veclen=0; veclen<2; ++veclen) { // vector length 0 (all) / 1 (vec)
            // set move
            //std::cout << "Setting SRRDType " << static_cast<size_t>(srrd) << ", veclen " << veclen << std::endl;
            mci.setTrialMove(srrd, veclen);
        }
    }
    //std::cout << std::endl;


    // now try to use (as example) customized student-t move and a multi-step move
    auto customStudentDist = std::student_t_distribution<double>(2);
    auto defaultUniformDist = std::uniform_real_distribution<double>(-1., 1.); // because we can
    StudentAllMove customAllMove(mci.getNDim(), 0.05, &customStudentDist);
    UniformVecMove customVecMove(mci.getNDim(), 1, 0.1, &defaultUniformDist);
    MultiStepMove customMultiMove(mci.getNDim());
    customMultiMove.setTrialMove(customVecMove); // add custom sub-move
    customMultiMove.addSamplingFunction( ExpNDPDF(mci.getNDim()) ); // add exponential sub-pdf (not the same as gauss!)

    const double origTargetRate = mci.getTargetAcceptanceRate();
    for (int i=0; i<2; ++i) {
        mci.setTargetAcceptanceRate(origTargetRate);
        if (i==0) {
            mci.setTrialMove(customAllMove);
        } else {
            mci.setTrialMove(customMultiMove);
            mci.setTargetAcceptanceRate(0.85); // should be higher for multi-step move
        }
        const int nskip = floor(1./mci.getTrialMove().getChangeRate());

        mci.clearObservables();
        mci.addObservable(obs1d, 1, nskip);
        mci.addObservable(obs3d, 1, nskip);
        mci.integrate(NMC*nskip, average, error, true, false);
        for (int i=0; i<mci.getNObsDim(); ++i) {
            //std::cout << "i " << i << ", average[i] " << average[i] << ", error[i] " << error[i] << ", CORRECT_RESULT" << CORRECT_RESULT << std::endl;
            assert( fabs(average[i]-CORRECT_RESULT) < 3.*error[i] );
        }
        //std::cout << "acceptance " << mci.getAcceptanceRate() << std::endl;
        //std::cout << std::endl;
    }


    return 0;
}
