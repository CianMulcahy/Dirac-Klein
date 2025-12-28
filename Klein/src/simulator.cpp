#include "simulator.h"

#include <complex>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <sys/stat.h>
#include <errno.h>

using cd = std::complex<double>;

Simulator::Simulator(Wavefunction &wf_, Hamiltonian &H_, double dt_, int steps_, double mass_)
    : wf(wf_), H(H_), dt(dt_), steps(steps_), mass(mass_) {

    absorbWidth = 120;       
    absorbStrength = 10.0;   
}

static inline Spinor applySigmaZ(const Spinor &s) {
    Spinor r;
    r.psi1 = s.psi1;
    r.psi2 = -s.psi2;
    return r;
}

void Simulator::computeRHS(const std::vector<Spinor> &psi, std::vector<Spinor> &rhs) const {
    int N = (int)psi.size();

    for (int i = 1; i < N - 1; ++i) {
        Spinor dpsi_dx;
        dpsi_dx.psi1 = (psi[i+1].psi1 - psi[i-1].psi1) / (2.0 * wf.dx);
        dpsi_dx.psi2 = (psi[i+1].psi2 - psi[i-1].psi2) / (2.0 * wf.dx);

        Spinor term1;
        term1.psi1 = -dpsi_dx.psi2;
        term1.psi2 = -dpsi_dx.psi1;

        Spinor sz = applySigmaZ(psi[i]);
        Spinor term2;
        term2.psi1 = cd(0.0, -mass) * sz.psi1;
        term2.psi2 = cd(0.0, -mass) * sz.psi2;

        double Vi = H.V[i];
        Spinor term3;
        term3.psi1 = cd(0.0, -Vi) * psi[i].psi1;
        term3.psi2 = cd(0.0, -Vi) * psi[i].psi2;

        rhs[i].psi1 = term1.psi1 + term2.psi1 + term3.psi1;
        rhs[i].psi2 = term1.psi2 + term2.psi2 + term3.psi2;
    }

    rhs[0].psi1 = rhs[0].psi2 = cd(0.0, 0.0);
    rhs[N-1].psi1 = rhs[N-1].psi2 = cd(0.0, 0.0);
}

void Simulator::applyAbsorbingBoundary() {

    int N = wf.N;
    int w = absorbWidth;
    if (w <= 0) return;
    if (2*w >= N) w = N/4;

    for (int i = 0; i < w; ++i) {
        double s = (double)(w - i) / (double)w; // 1 at edge -> 0 inward
        double damp = std::exp(-absorbStrength * s * s);
        wf.psi[i].psi1 *= damp;
        wf.psi[i].psi2 *= damp;
    }
    for (int i = N - w; i < N; ++i) {
        double s = (double)(i - (N - w)) / (double)w; // 0 inward -> 1 at edge
        double damp = std::exp(-absorbStrength * s * s);
        wf.psi[i].psi1 *= damp;
        wf.psi[i].psi2 *= damp;
    }
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
        wf.psi[i].psi1 += (dt / 6.0) * (k1[i].psi1 + cd(2.0,0.0)*k2[i].psi1 + cd(2.0,0.0)*k3[i].psi1 + k4[i].psi1);
        wf.psi[i].psi2 += (dt / 6.0) * (k1[i].psi2 + cd(2.0,0.0)*k2[i].psi2 + cd(2.0,0.0)*k3[i].psi2 + k4[i].psi2);
    }

    wf.psi[0].psi1 = wf.psi[0].psi2 = cd(0.0, 0.0);
    wf.psi[N-1].psi1 = wf.psi[N-1].psi2 = cd(0.0, 0.0);
    applyAbsorbingBoundary();
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
        cd a = wf.psi[i].psi1;
        cd b = wf.psi[i].psi2;
        double p1 = std::norm(a);
        double p2 = std::norm(b);
        double total = p1 + p2;

        file << x << ","
             << std::real(a) << "," << std::imag(a) << ","
             << std::real(b) << "," << std::imag(b) << ","
             << p1 << "," << p2 << "," << total << "\n";
    }
    file.close();
}

void Simulator::runSimulation(int saveEvery, int printEvery, const std::string &outPrefix) {
    // ensure output directory exists
    std::string dir = "data";
    if (outPrefix.find('/') != std::string::npos) {
        dir = outPrefix.substr(0, outPrefix.find('/'));
    }

    struct stat st = {0};
    if (stat(dir.c_str(), &st) == -1) {
        if (mkdir(dir.c_str(), 0755) != 0 && errno != EEXIST) {
            std::cerr << "Warning: could not create directory '" << dir << "'\n";
        }
    }

    for (int n = 0; n < steps; ++n) {
        timeStepRK4();

        if (n % printEvery == 0) {
            std::cout << "step " << n << "/" << steps
                      << "  norm=" << computeNorm()
                      << "  <x>=" << computeExpectationX()
                      << std::endl;
        }

        if (n % saveEvery == 0) {
            writeCSV(outPrefix + std::to_string(n) + ".csv");
        }
    }
}
