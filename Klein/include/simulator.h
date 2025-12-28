#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>
#include <string>
#include "wavefunction.h"
#include "hamiltonian.h"

class Simulator {
public:
    Wavefunction &wf;
    Hamiltonian &H;
    double dt;
    int steps;
    double mass;

    int absorbWidth;     // grid points per side
    double absorbStrength; // damping strength

    Simulator(Wavefunction &wf_, Hamiltonian &H_, double dt_, int steps_, double mass_);

    void timeStep();
    void runSimulation(int saveEvery = 10, int printEvery = 200, const std::string &outPrefix = "data/psi_");

    double computeNorm() const;
    double computeExpectationX() const;

    // x, Re(psi1), Im(psi1), Re(psi2), Im(psi2), |psi1|^2, |psi2|^2, total
    void writeCSV(const std::string &filename) const;

private:
    void computeRHS(const std::vector<Spinor> &psi, std::vector<Spinor> &rhs) const;
    void timeStepRK4();
    void applyAbsorbingBoundary(); 
};

#endif
