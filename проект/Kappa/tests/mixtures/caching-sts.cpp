/* 
 * Test for caching in mixture transport coefficients computation
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
	kappa::Mixture mixture_pure(molecules, "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml");

    std::vector<double> T_vals = { 300., 500., 1000., 2000., 10000., 40000. };
    double x_atom = 0.5;

    std::vector<arma::vec> mol_ndens;
    mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), mol));
    arma::vec atom_ndens(1);

    std::ofstream outf;
    std::ofstream outf2;
    outf.open("c:/Users/st024385/Documents/Code/kappa-tests/nocache_" + mol.name + "_" + at.name + "_xat50.txt");
    outf2.open("c:/Users/st024385/Documents/Code/kappa-tests/cache_" + mol.name + "_" + at.name + "_xat50.txt");
    outf << "T; lambda; eta" << std::endl;
    outf2 << "T; lambda; eta" << std::endl;
    double tot_ndens;
    double sh_v, th_c;

    for (auto T : T_vals) {
        tot_ndens =  101325.0 / (K_CONST_K * T);
        mol_ndens[0] = mixture.Boltzmann_distribution(T, (1 - x_atom) * tot_ndens, mol);
        atom_ndens[0] =  x_atom * tot_ndens;
        sh_v = mixture.thermal_conductivity(T, mol_ndens, atom_ndens);
        th_c = mixture.shear_viscosity(T, mol_ndens, atom_ndens);
        outf << T << ";" << sh_v << ";" << th_c << std::endl;
        mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);
        sh_v = mixture.get_thermal_conductivity();
        th_c = mixture.get_shear_viscosity();
        outf2 << T << ";" << sh_v << ";" << th_c << std::endl;
    }

    outf.close();
    outf2.close();

    std::string a;
	std::cout << "Enter anything to quit: ";
	std::cin >> a;
    return 0;
}

