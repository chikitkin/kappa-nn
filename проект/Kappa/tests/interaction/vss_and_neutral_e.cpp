#include <iostream>
#include "kappa.hpp"
#include <fstream>

using namespace std;
using namespace kappa;


int main(int argc, char** argv) {
    std::cout << "Start Test for Omega integrals, loading particle data" << endl;
	Approximation ApproximationTest{};

	Molecule N2("N2");
	Particle e("e-");
	Atom Np("N+");
	Atom N("N");

	Interaction N2Np(N2, Np);
	Interaction ee(e, e);
	Interaction Ne(N, e);
	// Interaction N2N(N2, N);
	// Interaction HH(H, H);
	std::vector<int> l_arr = {1, 1, 2, 3};
	std::vector<int> r_arr = {1, 4, 2, 3};
	int l, r;

	std::cout << "collision-reduced mass of two electrons: " << ee.collision_mass << ", e- mass:" << e.mass << std::endl;

	std::cout << Ne["_Om22_0"] << std::endl;

	for (int i=0;i<4; i++) {
		l = l_arr[i];
		r = r_arr[i];

		std::cout << N2.name << " + " << Np.name << std::endl;
		std::cout << "ESA: " << ApproximationTest.omega_integral(2000., N2Np, l, r, models_omega::model_omega_esa) << " " << ApproximationTest.omega_integral(10000., N2Np, l, r, models_omega::model_omega_esa) << std::endl;
		std::cout << "RS: " << ApproximationTest.omega_integral(2000., N2Np, l, r, models_omega::model_omega_rs) << " " << ApproximationTest.omega_integral(10000., N2Np, l, r, models_omega::model_omega_rs) << std::endl;
		std::cout << "VSS: " << ApproximationTest.omega_integral(2000., N2Np, l, r, models_omega::model_omega_vss) << " " << ApproximationTest.omega_integral(10000., N2Np, l, r, models_omega::model_omega_vss) << std::endl << std::endl;


		std::cout << N.name << " + " << e.name << std::endl;
		std::cout << "ESA: " << ApproximationTest.omega_integral(2000., Ne, l, r, models_omega::model_omega_esa) << " " << ApproximationTest.omega_integral(10000., Ne, l, r, models_omega::model_omega_esa) << std::endl;
		std::cout << "RS: " << ApproximationTest.omega_integral(2000., Ne, l, r, models_omega::model_omega_rs) << " " << ApproximationTest.omega_integral(10000., Ne, l, r, models_omega::model_omega_rs) << std::endl;

		std::cout << std::endl;
	}
	
    string a;
    cout << "Enter anything to quit: ";
    cin >> a;
    return 0;
}

