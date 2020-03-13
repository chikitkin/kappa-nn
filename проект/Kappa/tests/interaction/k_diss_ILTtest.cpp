/* 
 * File:   omegaTest.cpp
 */

#include <iostream>
#include "kappa.hpp"
#include <fstream>

using namespace std;
using namespace kappa;

void calculate(Molecule m, Atom a, Approximation appr) {
    std::cout << "Loading interaction parameters for " << m.name << "+" << a.name << endl;

    Interaction inter(m, a);

    ofstream outf;
    
    std::cout << "Calculating k_diss" << endl;

    std::vector<double> T_vals = { 500., 1000., 2000., 5000., 10000., 15000., 20000., 25000., 30000., 40000. };
    std::vector<int> i_vals;

    if (m.name == "O2") {
        i_vals = { 0, 5, 10, 20 };
    }
    else {
        i_vals = { 0, 5, 10, 20, 30 };
    }
    double res, tmp;
    outf.open(m.name + "_" + a.name + "_k_ILT_test.txt");

    outf << "T;";

    int counter;

    for (auto i : i_vals) {
        outf << "ILT,i=" << i << "(" << m.vibr_energy[0][i] / K_CONST_EV << " eV);";
        outf << "approx,i=" << i << "(" << m.vibr_energy[0][i] / K_CONST_EV << " eV);";
    }

    outf << endl;
    double c1, c2;
    
    for (auto T : T_vals) {
        outf << T << ";";

        for (auto i : i_vals) {
            outf << appr.k_diss(T, m, inter, i, model_k_diss_ilt) << ";"; 
            if (m.name == "O2") {
                c1 = 0.3867 * i * i * i - 2.7425 * i * i - 1901.9 * i + 61696;
                if (i <= 31) {
                    c2 = 1.63e-9 * i * i * i - 1.25e-7 * i * i + 3.24e-6 * i + 7.09e-5;
                }
                else if (i<=37) {
                    c2 = - 6.67e-6 * i * i + 4.65e-4 * i - 7.91e-3;
                }
                else {
                    c2 = 7.83e-7 * i * i * i * i - 1.31e-4 * i * i * i + 8.24e-3 * i * i - 0.23 * i + 2.4049;
                }
				outf << pow(T, -0.1) * 1.53e-10 * c2 * exp(-c1 / T) << ";";
            }
            else {
                if (i <= 8) {
                  c1 = 1.786e-18;
                }
                else if (i<=34) {
                  c1 = 1.71e-18;
                }
                else if (i<=52) {
                  c1 = 1.68e-18;
                }
                else {
                  c1 = 1.66e-18;
                }
                c2 = 4e-19 * i * i * i * i + 5.24e-19 * i * i * i - 7.41e-17 * i * i + 6.42e-15 * i + 7.3e-14;
                outf << 7.16e-2 * pow(T, -0.25) * c2 * exp((m.vibr_energy[0][i] - c1) / (K_CONST_K * T)) << ";";
            }
        }
        outf << endl;
    }
}

int main(int argc, char** argv) {
    std::cout << "Start Test for k_diss coefficients, loading particle data" << endl;

    Molecule N2("N2", "harmonic");
    Molecule O2("O2", "harmonic");
    Atom N("N");
    Atom O("O");
    Approximation ApproximationTest{};

    calculate(N2, N, ApproximationTest);
    calculate(O2, O, ApproximationTest);

    string a;
    cout << "Enter anything to quit: ";
    cin >> a;
    return 0;
}

