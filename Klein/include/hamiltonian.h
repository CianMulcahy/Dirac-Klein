#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <vector>

class Hamiltonian {
public:
    std::vector<double> V;
    int N;
    double dx;

    Hamiltonian(int N_, double dx_);

    void setBarrier(double V0, int startIdx, int endIdx);

    void setStep(double V0, int stepIdx);

    void setZero();
};

#endif
