#include "simulator.h"

#include <complex>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <sys/stat.h> // mkdir
#include <errno.h>

using cd = std::complex<double>;

Simulator::Simulator(Wavefunction &wf_, Hamiltonian &H_, double dt_, int steps_, double mass_)
    : wf(wf_), H(H_), dt(dt_), steps(steps_), mass(mass_) {

    if ((int)H.V.size() != wf.N) {
        std::cerr << "Warning: Hamiltonian V size does not match wavefunction N\n";
    }
}

static inline Spinor applySigmaZ(const Spinor &s) {
    Spinor r;
    r.psi1 = s.psi1;
    r.psi2 = -s.psi2;
    return r;
}

void Simulator::computeRHS(const std::vector<Spinor> &psi, std::vector<Spinor> &rhs) const {
    int N = psi.size();

    for (int i = 1; i < N - 1; ++i) {
        Spinor dpsi_dx;
        dpsi_dx.psi1 = (psi[i+1].psi1 - psi[i-1].psi1) / (2.0 * wf.dx);
        dpsi_dx.psi2 = (psi[i+1].psi2 - psi[i-1].psi2) / (2.0 * wf.dx);

        Spinor term1;
        term1.psi1 = -dpsi_dx.psi2;
        term1.psi2 = -dpsi_dx.psi1;

        // term2
        Spinor sz = applySigmaZ(psi[i]);
        Spinor term2;
        term2.psi1 = cd(0.0, -mass) * sz.psi1; // -i m * sz.psi1
        term2.psi2 = cd(0.0, -mass) * sz.psi2; // -i m * sz.psi2

        double Vi = (i < (int)H.V.size() ? H.V[i] : 0.0);
        Spinor term3;
        term3.psi1 = cd(0.0, -Vi) * psi[i].psi1;
        term3.psi2 = cd(0.0, -Vi) * psi[i].psi2;

        // sum
        rhs[i].psi1 = term1.psi1 + term2.psi1 + term3.psi1;
        rhs[i].psi2 = term1.psi2 + term2.psi2 + term3.psi2;
    }

    rhs[0].psi1 = rhs[0].psi2 = cd(0.0, 0.0);
    rhs[N-1].psi1 = rhs[N-1].psi2 = cd(0.0, 0.0);
}

void Simulator::timeStepRK4() {
    int N = wf.N;
    std::vector<Spinor> k1(N), k2(N), k3(N), k4(N), temp(N);

    computeRHS(wf.psi, k1);

    for (int i = 0; i < N; ++i) {
        temp[i].psi1 = wf.psi[i].psi1 + k1[i].psi1 * (dt * 0.5);
        temp[i].psi2 = wf.psi[i].psi2 + k1[i].psi2 * (dt * 0.5);
    }
    computeRHS(temp, k2);

    for (int i = 0; i < N; ++i) {
        temp[i].psi1 = wf.psi[i].psi1 + k2[i].psi1 * (dt * 0.5);
        temp[i].psi2 = wf.psi[i].psi2 + k2[i].psi2 * (dt * 0.5);
    }
    computeRHS(temp, k3);

    for (int i = 0; i < N; ++i) {
        temp[i].psi1 = wf.psi[i].psi1 + k3[i].psi1 * dt;
        temp[i].psi2 = wf.psi[i].psi2 + k3[i].psi2 * dt;
    }
    computeRHS(temp, k4);

    for (int i = 0; i < N; ++i) {
        wf.psi[i].psi1 += (dt / 6.0) * (k1[i].psi1 + cd(2.0, 0.0) * k2[i].psi1 + cd(2.0, 0.0) * k3[i].psi1 + k4[i].psi1);
        wf.psi[i].psi2 += (dt / 6.0) * (k1[i].psi2 + cd(2.0, 0.0) * k2[i].psi2 + cd(2.0, 0.0) * k3[i].psi2 + k4[i].psi2);
    }

    wf.psi[0].psi1 = wf.psi[0].psi2 = cd(0.0, 0.0);
    wf.psi[N-1].psi1 = wf.psi[N-1].psi2 = cd(0.0, 0.0);
}

void Simulator::timeStep() {
    timeStepRK4();
}

double Simulator::computeNorm() const {
    double sum = 0.0;
    for (int i = 0; i < wf.N; ++i) {
        sum += std::norm(wf.psi[i].psi1) + std::norm(wf.psi[i].psi2);
    }
    return sum * wf.dx;
}

double Simulator::computeExpectationX() const {
    double sum = 0.0;
    for (int i = 0; i < wf.N; ++i) {
        double x = i * wf.dx;
        sum += x * (std::norm(wf.psi[i].psi1) + std::norm(wf.psi[i].psi2));
    }
    return sum * wf.dx;
}

void Simulator::writeCSV(const std::string &filename) const {
    std::ofstream file(filename);
    file << std::setprecision(12);
    for (int i = 0; i < wf.N; ++i) {
        double x = i * wf.dx;
        double p1 = std::norm(wf.psi[i].psi1);
        double p2 = std::norm(wf.psi[i].psi2);
        file << x << "," << p1 << "," << p2 << "," << (p1 + p2) << "\n";
    }
    file.close();
}

void Simulator::runSimulation(int saveEvery, int printEvery, const std::string &outPrefix) {

    std::string dir = "data";
    if (outPrefix.find('/') != std::string::npos) {
        dir = outPrefix.substr(0, outPrefix.find('/'));
    }

    struct stat st = {0};
    if (stat(dir.c_str(), &st) == -1) {
        if (mkdir(dir.c_str(), 0755) != 0) {
            if (errno != EEXIST) {
                std::cerr << "Warning: could not create directory '" << dir << "'\n";
            }
        }
    }

    for (int n = 0; n < steps; ++n) {
        timeStepRK4();

        if ( (n % printEvery) == 0 ) {
            double norm = computeNorm();
            double xexp = computeExpectationX();
            std::cout << "step " << n << " / " << steps << "  norm=" << norm << "  <x>=" << xexp << std::endl;
        }

        if ( (n % saveEvery) == 0 ) {
            std::string fname = outPrefix + std::to_string(n) + ".csv";
            writeCSV(fname);
        }
    }
}
