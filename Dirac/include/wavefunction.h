#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include <vector>
#include <complex>

using cd = std::complex<double>;

struct Spinor {

    cd psi1;
    cd psi2;
    Spinor() : psi1(0), psi2(0) {}

};

class Wavefunction {
public:
    std::vector<Spinor> psi;
    int N;
    double dx;

    Wavefunction(int N_, double dx_);
    void initializeGaussian(double x0, double k0, double sigma);
};

#endif