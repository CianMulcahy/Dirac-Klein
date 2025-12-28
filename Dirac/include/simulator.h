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

    Simulator(Wavefunction &wf_, Hamiltonian &H_, double dt_, int steps_, double mass_);

    void timeStep();

    void runSimulation(int saveEvery = 10, int printEvery = 50, const std::string &outPrefix = "data/psi_");

    double computeNorm() const;

    double computeExpectationX() const;

    void writeCSV(const std::string &filename) const;

private:

    void computeRHS(const std::vector<Spinor> &psi, std::vector<Spinor> &rhs) const;
    void timeStepRK4();
};

#endif // SIMULATOR_H
