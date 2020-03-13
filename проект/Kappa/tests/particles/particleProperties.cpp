/* 
 * File:   particleProperties.cpp
 */

#include <iostream>

#include "kappa.hpp"

int main(int argc, char** argv) {
    // ---------------------------------------
	std::cout << "Particle properties test" << std::endl;
    kappa::Molecule O2 = kappa::Molecule("O2");
	std::cout << "Finished loading O2" << std::endl;
    kappa::Molecule N2 = kappa::Molecule("N2");
	std::cout << "Finished loading N2" << std::endl;
    
    kappa::Atom O = kappa::Atom("O");
	std::cout << "Finished loading O" << std::endl;
    kappa::Atom N = kappa::Atom("N");
	std::cout << "Finished loading N" << std::endl << std::endl;

	kappa::Approximation appr = kappa::Approximation{};

	double dE = kappa::convert_cm_to_Joule(1000);
	std::cout << "N, n_el_levels = " << appr.max_electron_level(N, dE) << std::endl;
	std::cout << "O, n_el_levels = " << appr.max_electron_level(O, dE) << std::endl;
    
    std::cout << "N2, n_el_levels = " << N2.num_electron_levels << std::endl;
    std::cout << "O2, n_el_levels = " << O2.num_electron_levels << std::endl;
    std::cout << "N, n_el_levels (FULL) = " << N.num_electron_levels << std::endl;
    std::cout << "O, n_el_levels (FULL) = " << O.num_electron_levels << std::endl;
    
    int i;
    for (i=0; i<N2.num_electron_levels; i++) {
		std::cout << "Electron level = " << i << std::endl;
        std::cout << "N2, n vibr levels = " << N2.num_vibr_levels[i] << std::endl;
        std::cout << "N2, n rot levels = " << N2.num_rot_levels[i][0] << std::endl << std::endl;
    }
    
    for (i=0; i<O2.num_electron_levels; i++) {
		std::cout << "Electron level = " << i << std::endl;
        std::cout << "O2, n vibr levels = " << O2.num_vibr_levels[i] << std::endl;
        std::cout << "O2, n rot levels = " << O2.num_rot_levels[i][0] << std::endl << std::endl;
    }

	std::cout << "N2, n vibr levels in ground state = " << N2.num_vibr_levels[0] << std::endl;
	std::cout << "O2, n vibr levels in ground state = " << O2.num_vibr_levels[0] << std::endl << std::endl;

	std::cout << "N2 frequency of ground state = " << N2.vibr_frequency[0] << std::endl;

	std::cout << "N2, vibr energy in ground state size = " << N2.vibr_energy[0].size() << std::endl;
	std::cout << "O2, vibr energy in ground state size = " << O2.vibr_energy[0].size() << std::endl << std::endl;

	std::cout << "N2 vibr energy of last level = " << N2.vibr_energy[0][N2.num_vibr_levels[0] - 1] << std::endl;
	std::cout << "O2 vibr energy of last level = " << O2.vibr_energy[0][O2.num_vibr_levels[0] - 1] << std::endl << std::endl;

	std::cout << "N2 Ediss = " << N2.diss_energy[0] << std::endl;
	std::cout << "O2 Ediss = " << O2.diss_energy[0] << std::endl << std::endl;

	std::cout << "N2 Ediss - E_el = " << N2.diss_energy - N2.electron_energy << std::endl << std::endl;

	std::cout << "O2 Ediss - E_el = " << O2.diss_energy - O2.electron_energy << std::endl << std::endl;
    
	std::string a;

	std::cin >> a;

    return 1;
}

