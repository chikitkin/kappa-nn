/* 
 * File:   mixture-sts-shear.cpp
 * Test for shear viscosity in binary mixtures
 */

#include <iostream>
#include "kappa.hpp"
#include <fstream>


int main(int argc, char** argv) {
    std::cout << "Start test: computation of shear viscosity" << std::endl;

    std::vector<kappa::Molecule> molecules;
    std::vector<kappa::Atom> atoms;

    std::cout << "Loading particles data" << std::endl;

    kappa::Molecule mol("N2", true, true, "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
	kappa::Atom at("N", "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
    molecules.push_back(mol);
    atoms.push_back(at);

	std::cout << "Finished loading particles data" << std::endl;

	kappa::Mixture mixture(molecules, atoms, "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml");

    std::vector<double> T_vals = { 500., 1000., 5000., 10000., 20000., 40000. };
    double x_atom = 0.5;

    std::vector<arma::vec> mol_ndens;
    mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), mol));
    arma::vec atom_ndens(1);

    std::ofstream outf;
    outf.open("c:/Users/st024385/Documents/Code/kappa-tests/shear_" + mol.name + "_" + at.name + "_xat50.txt");
    outf << "T; eta; H_N2N2; H_N2N; H_NN; omega_11_N2N2; omega_22_N2N2; omega_11_NN; omega_22_NN" << std::endl;
    double tot_ndens;

    for (auto T : T_vals) {
        tot_ndens =  101325.0 / (K_CONST_K * T);
        mol_ndens[0] = mixture.Boltzmann_distribution(T, (1 - x_atom) * tot_ndens, mol);
        atom_ndens[0] = x_atom * tot_ndens;
        outf << T << ";" << mixture.shear_viscosity(T, mol_ndens, atom_ndens) << std::endl;
    }

    outf.close();
    std::string a;
	std::cout << "Enter anything to quit: ";
	std::cin >> a;
    return 0;
}

