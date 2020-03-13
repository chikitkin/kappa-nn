#include <iostream>
#include "kappa.hpp"
#include <chrono>

using namespace std;
using namespace kappa;

int main(int argc, char** argv) {
    std::cout << "Start Test for k_exch coefficients, loading particle data" << endl;

    Molecule Molecule1("N2");
    Atom Atom1("O");

    std::cout << "Loading interaction parameters" << endl;

    Interaction inter(Molecule1, Atom1);
    Approximation ApproximationTest{};

    ofstream outf;
	
	double T = 0.0;
	double res;
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    for (int j=0; j<100; j++) {
        T = 500 + j * 100;
        res = 0.0;
        for (int i=0; i<Molecule1.num_vibr_levels[0]; i++) {
            res += ApproximationTest.k_exch(T, Molecule1, Atom1, inter, i, models_k_exch::model_k_exch_polak);
        }
        std::cout << res << std::endl;
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    std::cout << duration << std::endl;
    
    string a;
    cout << "Enter anything to quit: ";
    cin >> a;
    return 0;
}

