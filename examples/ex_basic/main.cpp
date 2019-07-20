#include <iostream>

#include "mci/MCIntegrator.hpp"

#include "../common/ExampleFunctions.hpp" // look here for the used example functions

int main()
{
    using namespace std;
    using namespace mci; // everything in MCI lib is namespaced

    // intro
    cout << "We want to compute the integral" << endl;
    cout << "    Integral[-1:3] dx (4x - x^2)" << endl << endl;
    cout << "Notice that we can re-write the integral as" << endl;
    cout << "    Integral[-1:3] dx (5*sign(x)*(4-x) * |x|/5)" << endl;
    cout << "where g(x)=|x|/5 is a candidate sampling function since it's positive on the domain [-1:3] and its integral on the domain is equal to 1." << endl << endl;

    cout << "We start by setting the MCI:" << endl;


    // declare a 1-dimensional integrator
    const int ndim = 1;
    MCI mci(ndim);

    cout << "ndim = " << mci.getNDim() << endl;


    // set the integration range to [-1:3]
    mci.setIRange(-1., 3.);

    cout << "irange = [ " << -1 << " ; " << 3 << " ]" << endl;


    // initial walker position
    double initpos[ndim];
    initpos[0] = -0.5;
    mci.setX(initpos); // alternatively use mci.setX(0, -0.5);

    cout << "initial walker position = " << mci.getX(0) << endl;


    // initial MRT2 step
    mci.setMRT2Step(0.25);

    cout << "MRT2 step = " << mci.getMRT2Step(0) << endl;


    // target acceptance rate
    mci.setTargetAcceptanceRate(0.7);

    cout << "Acceptance rate = " << mci.getTargetAcceptanceRate() << endl;
    cout << endl << endl;


    // first way of integrating
    cout << "We first compute the integral without setting a sampling function." << endl;
    cout << "    f(x) = (4x - x^2) " << endl;
    cout << "    g(x) = - " << endl << endl;


    // observable
    Parabola obs; // we create obs on the stack
    mci.addObservable(obs); // obs will be cloned by MCI
    // instead you may pass unique pointers, as shown later

    cout << "Number of observables set = " << mci.getNObs() << endl;
    cout << "Dimension of observables set = " << mci.getNObsDim() << endl;

    // sampling function
    cout << "Number of sampling function set = " << mci.getNPDF() << endl;


    // integrate
    const int Nmc = 1000000;
    double average[mci.getNObsDim()];
    double error[mci.getNObsDim()];
    mci.integrate(Nmc, average, error);

    cout << "The integral gives as result = " << average[0] << "   +-   " << error[0] << endl;
    cout << "--------------------------------------------------------" << endl << endl;


    // second way of integrating
    cout << "Now we compute the integral using a sampling function." << endl;
    cout << "    f(x) = 5*sign(x)*(4-x) " << endl;
    cout << "    g(x) = |x|/5 " << endl << endl;


    // observable
    mci.clearObservables();  // we first remove the old observable
    mci.addObservable(std::make_unique<NormalizedParabola>()); // you may also directly transfer a unique pointer into MCI

    cout << "Number of observables set = " << mci.getNObs() << endl;
    cout << "Dimension of observables set = " << mci.getNObsDim() << endl;


    // sampling function
    std::unique_ptr<SamplingFunctionInterface> pdf = std::make_unique<NormalizedLine>(); // if you create the unique pointer beforehand
    mci.addSamplingFunction(std::move(pdf)); // you have to explicitly move it in
    // But remember that the local pdf ptr is now invalid and may not be used!

    cout << "Number of sampling function set = " << mci.getNPDF() << endl;


    // integrate
    mci.integrate(Nmc, average, error);

    cout << "The integral gives as result = " << average[0] << "   +-   " << error[0] << endl;
    cout << "--------------------------------------------------------" << endl << endl;


    // final comments
    cout << "Using a sampling function in this case gives worse performance. In fact, the error bar is larger." << endl;
    cout << "This implies that the variance of the re-factored f(x) written for introducing a sampling function, is larger than the original f(x)." << endl;

    // if we wanted to safe the created sampling function from construction we could do:
    cout << endl << endl << "Before MCI gets deleted, let's retrieve some objects that we added..." << endl;
    pdf = mci.popSamplingFunction(); // pdf is now valid again and contains the original
    auto obs2 = mci.popObservable(); // we may also obtain the observable that we never owned directly
    // NOTE: explicitly, obs2 is of type std::unique_ptr<ObservableFunctionInterface>

    cout << "OK. Now MCI is empty:" << endl << endl;
    cout << "Number of observables set = " << mci.getNObs() << endl;
    cout << "Dimension of observables set = " << mci.getNObsDim() << endl;
    cout << "Number of sampling function set = " << mci.getNPDF() << endl;

    // end (we didn't use unguarded "new" allocation, so we don't need to delete anything)
    return 0;
}
