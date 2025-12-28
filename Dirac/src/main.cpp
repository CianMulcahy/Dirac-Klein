#include "wavefunction.h"
#include "hamiltonian.h"
#include "simulator.h"

int main() {
    int N = 800;
    double dx = 0.05;
    double dt = 0.002;
    int steps = 2000;
    double mass = 1.0;

    Wavefunction wf(N, dx);

    double x0 = 0.25 * N * dx;
    double k0 = 2.0;
    double sigma = 1.0;
    wf.initializeGaussian(x0, k0, sigma);

    Hamiltonian H(N, dx);

    int mid = N / 2;
    H.setBarrier(5.0, mid-10, mid+10);

    Simulator sim(wf, H, dt, steps, mass);
    sim.runSimulation(10, 100, "data/psi_"); // save every 10 steps, print every 100

    return 0;
}
