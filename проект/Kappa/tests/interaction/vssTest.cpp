#include <iostream>
#include "kappa.hpp"
#include <fstream>

using namespace std;
using namespace kappa;


int main(int argc, char** argv) {
    std::cout << "Start Test for Omega integrals, loading particle data" << endl;
	Approximation ApproximationTest{};

	Molecule N2("N2");
	Atom N("N");
	Atom H("H");

	Interaction N2N2(N2, N2);
	Interaction N2N(N2, N);
	Interaction HH(H, H);

	std::cout << "ESA: " << ApproximationTest.omega_integral(2000., N2N2, 1, 1, models_omega::model_omega_esa) << " " << ApproximationTest.omega_integral(10000., N2N2, 1, 1, models_omega::model_omega_esa) << std::endl;
	std::cout << "RS: " << ApproximationTest.omega_integral(2000., N2N2, 1, 1, models_omega::model_omega_rs) << " " << ApproximationTest.omega_integral(10000., N2N2, 1, 1, models_omega::model_omega_rs) << std::endl;
	std::cout << "VSS: " << ApproximationTest.omega_integral(2000., N2N2, 1, 1, models_omega::model_omega_vss) << " " << ApproximationTest.omega_integral(10000., N2N2, 1, 1, models_omega::model_omega_vss) << std::endl << std::endl;


	std::cout << "ESA: " << ApproximationTest.omega_integral(2000., N2N2, 2, 3, models_omega::model_omega_esa) << " " << ApproximationTest.omega_integral(10000., N2N2, 2, 3, models_omega::model_omega_esa) << std::endl;
	std::cout << "RS: " << ApproximationTest.omega_integral(2000., N2N2, 2, 3, models_omega::model_omega_rs) << " " << ApproximationTest.omega_integral(10000., N2N2, 2, 3, models_omega::model_omega_rs) << std::endl;
	std::cout << "VSS: " << ApproximationTest.omega_integral(2000., N2N2, 2, 3, models_omega::model_omega_vss) << " " << ApproximationTest.omega_integral(10000., N2N2, 2, 3, models_omega::model_omega_vss) << std::endl << std::endl;



	std::cout << "ESA: " << ApproximationTest.omega_integral(2000., N2N, 1, 1, models_omega::model_omega_esa) << " " << ApproximationTest.omega_integral(10000., N2N, 1, 1, models_omega::model_omega_esa) << std::endl;
	std::cout << "RS: " << ApproximationTest.omega_integral(2000., N2N, 1, 1, models_omega::model_omega_rs) << " " << ApproximationTest.omega_integral(10000., N2N, 1, 1, models_omega::model_omega_rs) << std::endl;
	std::cout << "VSS: " << ApproximationTest.omega_integral(2000., N2N, 1, 1, models_omega::model_omega_vss) << " " << ApproximationTest.omega_integral(10000., N2N, 1, 1, models_omega::model_omega_vss) << std::endl << std::endl;


	std::cout << "ESA: " << ApproximationTest.omega_integral(2000., N2N, 2, 3, models_omega::model_omega_esa) << " " << ApproximationTest.omega_integral(10000., N2N, 2, 3, models_omega::model_omega_esa) << std::endl;
	std::cout << "RS: " << ApproximationTest.omega_integral(2000., N2N, 2, 3, models_omega::model_omega_rs) << " " << ApproximationTest.omega_integral(10000., N2N, 2, 3, models_omega::model_omega_rs) << std::endl;
	std::cout << "VSS: " << ApproximationTest.omega_integral(2000., N2N, 2, 3, models_omega::model_omega_vss) << " " << ApproximationTest.omega_integral(10000., N2N, 2, 3, models_omega::model_omega_vss) << std::endl << std::endl;


	try {
		std::cout << ApproximationTest.omega_integral(2000., HH, 2, 3, models_omega::model_omega_vss);
	}
	catch (const DataNotFoundException &e) {
        std::cout << e.what() << endl;
    }

    string a;
    cout << "Enter anything to quit: ";
    cin >> a;
    return 0;
}

