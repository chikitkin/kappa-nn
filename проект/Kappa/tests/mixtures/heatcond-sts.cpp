/* 
 * File:   heatcond-sts.cpp
 * Test for heat conductivity in binary mixtures and a pure gas
 */

#include <iostream>
#include "kappa.hpp"
#include <fstream>


int main(int argc, char** argv) {
    std::cout << "Start test: computation of heat conductivity" << std::endl;
	
    std::cout << "Loading particles data" << std::endl;

    kappa::Molecule mol("O2", true, true, "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
    kappa::Molecule mol_nonrig("O2", true, false, "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
	kappa::Atom at("O", "c:/Users/st024385/Documents/DB/Particles/particles.yaml");

	std::cout << "Finished loading particles data" << std::endl;
	std::cout << "Molecule vibrational levels " << mol.num_vibr_levels[0] << std::endl;

	kappa::Mixture mixture(mol, at, "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml");
	kappa::Mixture mixture_pure(mol, "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml");
    kappa::Mixture mixture_pure_at(at, "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml");
    kappa::Mixture mixture_nonrig(mol_nonrig, at, "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml");
    kappa::Mixture mixture_nonrig_pure(mol_nonrig, "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml");

	std::vector<double> T_vals;

	int i;
	for (i = 0; i < 159; i++) {
		T_vals.push_back(500 + i * 250);
	}
	int x_atom_perc = 50;

    double x_atom = x_atom_perc / 100.;
    std::vector<arma::vec> mol_ndens;
    mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), mol));
    arma::vec atom_ndens(1);

    std::ofstream outf, outf2, outf3, outf4, outf5;
    outf.open("c:/Users/st024385/Documents/Code/kappa-tests/thermal_conductivity/thcond_" + mol.name + "_" + at.name + "_xat" + std::to_string(x_atom_perc) + ".txt");
	outf2.open("c:/Users/st024385/Documents/Code/kappa-tests/thermal_conductivity/thcond_pure_" + mol.name + ".txt");
    outf3.open("c:/Users/st024385/Documents/Code/kappa-tests/thermal_conductivity/thcond_pure_" + at.name + ".txt");
    outf4.open("c:/Users/st024385/Documents/Code/kappa-tests/thermal_conductivity/thcond_sts_" + mol.name + "_" + at.name + "_xat" + std::to_string(x_atom_perc) + ".txt");
    outf5.open("c:/Users/st024385/Documents/Code/kappa-tests/thermal_conductivity/thcond_pure_sts_" + mol.name + ".txt");
    outf << "T; lambda;" << std::endl;
	outf2 << "T; lambda;" << std::endl;
    outf3 << "T; lambda;" << std::endl;
    outf4 << "T; lambda;" << std::endl;
    outf5 << "T; lambda;" << std::endl;
    double tot_ndens;

    for (auto T : T_vals) {
        tot_ndens =  101325.0 / (K_CONST_K * T);
        mol_ndens[0] = mixture.Boltzmann_distribution(T, (1 - x_atom) * tot_ndens, mol);
        atom_ndens[0] = x_atom * tot_ndens;
        std::cout << T << std::endl;

		mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);
		mixture_nonrig.compute_transport_coefficients(T, mol_ndens, atom_ndens);

        outf << T << ";" << mixture.get_thermal_conductivity() << std::endl;
        outf4 << T << ";" << mixture_nonrig.get_thermal_conductivity() << std::endl;
		
		mol_ndens[0] = mixture.Boltzmann_distribution(T, tot_ndens, mol);

		mixture_pure.compute_transport_coefficients(T, mol_ndens);
		mixture_nonrig_pure.compute_transport_coefficients(T, mol_ndens);
		mixture_pure_at.compute_transport_coefficients(T, atom_ndens);
		outf2 << T << ";" << mixture_pure.get_thermal_conductivity() << std::endl;
        outf5 << T << ";" << mixture_nonrig_pure.get_thermal_conductivity() << std::endl;
        
		atom_ndens[0] = tot_ndens;
		outf3 << T << ";" << mixture_pure_at.get_thermal_conductivity() << std::endl;
    }

    outf.close();
	outf2.close();
	outf3.close();
    outf4.close();
    outf5.close();

    std::string a;
	std::cout << "Enter anything to quit: ";
	std::cin >> a;
    return 0;
}

