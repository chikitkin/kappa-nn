/* 
 * Test for heat conductivity in binary mixtures; molecules are assumed to be non-rigid rotators
 */

#include <iostream>
#include "kappa.hpp"
#include <fstream>


int main(int argc, char** argv) {
    std::cout << "Start test: computation of shear viscosity" << std::endl;

    std::vector<kappa::Molecule> molecules;
    std::vector<kappa::Atom> atoms;

    std::cout << "Loading particles data" << std::endl;

    kappa::Molecule mol("N2", true, false, "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
	kappa::Atom at("N", "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
    molecules.push_back(mol);
    atoms.push_back(at);

	std::cout << "Finished loading particles data" << std::endl;

    kappa::Mixture mixture_pure(molecules, "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml");
    kappa::Mixture mixture(molecules, atoms, "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml");
    std::vector<double> T_vals = { 300., 500., 1000., 2000., 10000., 40000. };
    double x_atom = 0.2;

    std::vector<arma::vec> mol_ndens;
    mol_ndens.push_back(mixture_pure.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), mol));
    arma::vec atom_ndens(1);

    std::ofstream outf, outf2;
    outf.open("c:/Users/st024385/Documents/Code/kappa-tests/heatcond_nonrig_" + mol.name + "_" + at.name + "_pure.txt");
    outf << "T; lambda;" << std::endl;
    outf2.open("c:/Users/st024385/Documents/Code/kappa-tests/heatcond_nonrig_" + mol.name + "_" + at.name + "_xat20.txt");
    outf2 << "T; lambda;" << std::endl;
    double tot_ndens;

    for (auto T : T_vals) {
        tot_ndens =  101325.0 / (K_CONST_K * T);
        mol_ndens[0] = mixture_pure.Boltzmann_distribution(T, tot_ndens, mol);
        std::cout << T << std::endl;
		mixture_pure.compute_transport_coefficients(T, mol_ndens);
        outf << T << ";" << mixture_pure.get_thermal_conductivity() << std::endl;


        mol_ndens[0] = mixture.Boltzmann_distribution(T, (1 - x_atom) * tot_ndens, mol);
        atom_ndens[0] = x_atom * tot_ndens;
		mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);
        outf2 << T << ";" << mixture.get_thermal_conductivity() << std::endl;
    }


    outf.close();
    outf2.close();
	std::string a;
	std::cout << "Enter anything to quit: ";
	std::cin >> a;
    return 0;
}

