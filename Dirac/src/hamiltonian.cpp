#include "hamiltonian.h"

Hamiltonian::Hamiltonian(int N_, double dx_) : N(N_), dx(dx_) {
    V.resize(N, 0.0);
}

void Hamiltonian::setBarrier(double V0, int startIdx, int endIdx){
    for(int i=startIdx; i<=endIdx; ++i) V[i] = V0;
}