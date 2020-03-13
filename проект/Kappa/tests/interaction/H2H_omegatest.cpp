/* 
 * File:   omegaTest.cpp
 */

#include <iostream>
#include "kappa.hpp"
#include <fstream>

using namespace std;
using namespace kappa;


int main(int argc, char** argv) {
    std::cout << "Start Test for Omega integrals, loading particle data" << endl;
	Approximation ApproximationTest{};

	Molecule H2("H2");

	Atom H("H");

	Interaction H2H2(H2, H2);
	Interaction H2H(H2, H);
	Interaction HH(H, H);

	std::cout << ApproximationTest.omega_integral(2000., H2H2, 1, 1, models_omega::model_omega_esa) << " " << ApproximationTest.omega_integral(10000., H2H2, 1, 1, models_omega::model_omega_esa) << std::endl;
	std::cout << ApproximationTest.omega_integral(2000., H2H2, 1, 1, models_omega::model_omega_rs) << " " << ApproximationTest.omega_integral(10000., H2H2, 1, 1, models_omega::model_omega_rs) << std::endl;

	std::cout << ApproximationTest.omega_integral(2000., H2H, 1, 1, models_omega::model_omega_esa) << " " << ApproximationTest.omega_integral(10000., H2H, 1, 1, models_omega::model_omega_esa) << std::endl;
	std::cout << ApproximationTest.omega_integral(2000., H2H, 1, 1, models_omega::model_omega_rs) << " " << ApproximationTest.omega_integral(10000., H2H, 1, 1, models_omega::model_omega_rs) << std::endl;
	
	std::cout << ApproximationTest.omega_integral(2000., HH, 1, 1, models_omega::model_omega_esa) << " " << ApproximationTest.omega_integral(10000., HH, 1, 1, models_omega::model_omega_esa) << std::endl;
	std::cout << ApproximationTest.omega_integral(2000., HH, 1, 1, models_omega::model_omega_rs) << " " << ApproximationTest.omega_integral(10000., HH, 1, 1, models_omega::model_omega_rs) << std::endl;

    string a;
    cout << "Enter anything to quit: ";
    cin >> a;
    return 0;
}

