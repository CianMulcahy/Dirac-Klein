#include "hamiltonian.h"

Hamiltonian::Hamiltonian(int N_, double dx_) : N(N_), dx(dx_) {
    V.assign(N, 0.0);
}

void Hamiltonian::setBarrier(double V0, int startIdx, int endIdx) {
    for (int i = 0; i < N; ++i) V[i] = 0.0;
    if (startIdx < 0) startIdx = 0;
    if (endIdx >= N) endIdx = N - 1;
    for (int i = startIdx; i <= endIdx; ++i) V[i] = V0;
}

void Hamiltonian::setStep(double V0, int stepIdx) {
    for (int i = 0; i < N; ++i) V[i] = (i >= stepIdx ? V0 : 0.0);
}

void Hamiltonian::setZero() {
    for (int i = 0; i < N; ++i) V[i] = 0.0;
}
