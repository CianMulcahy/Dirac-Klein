#include "wavefunction.h"
#include "hamiltonian.h"
#include "simulator.h"

int main() {
    int N = 1400;     
    double dx = 0.05;   

    double m = 1.0;      // mass
    double k0 = 2.2;     // initial momm
    double sigma = 1.0;  // width of Gaussian
    double x0 = 12.0;    // initial center 

    int stepIdx = (int)((N * dx) * 0.45 / dx); 
    double V0 = 8.0;

    double dt = 0.0015; 
    int steps = 22000;   
    int saveEvery = 8;   
    int printEvery = 400;

    Wavefunction wf(N, dx);
    wf.initializeGaussian(x0, k0, sigma);

    Hamiltonian H(N, dx);
    H.setStep(V0, stepIdx);

    Simulator sim(wf, H, dt, steps, m);

    sim.absorbWidth = 140;
    sim.absorbStrength = 12.0;

    sim.runSimulation(saveEvery, printEvery, "data/klein_");

    return 0;
}
