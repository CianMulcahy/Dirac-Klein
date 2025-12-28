#include "wavefunction.h"
#include <cmath>
#include <complex>

Wavefunction::Wavefunction(int N_, double dx_) : N(N_), dx(dx_) {
    psi.resize(N);
}

void Wavefunction::initializeGaussian(double x0, double k0, double sigma) {
    for(int i=0; i<N; ++i){
        double x = i*dx;
        double gauss = std::exp(-(x-x0)*(x-x0)/(2*sigma*sigma));
        std::complex<double> phase = std::exp(cd(0, k0*x));
        psi[i].psi1 = gauss * phase;
        psi[i].psi2 = 0.0;
    }
}