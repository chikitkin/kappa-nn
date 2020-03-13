#include <iostream>
#include "kappa.hpp"
#include <fstream>


int main(int argc, char** argv) {
    std::cout << "Start test: computation of bulk viscosity" << std::endl;

    std::vector<kappa::Molecule> molecules;
    std::vector<kappa::Atom> atoms;

    std::cout << "Loading particles data" << std::endl;

    kappa::Molecule mol("N2", true, false, "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
	kappa::Atom at("N", "c:/Users/st024385/Documents/DB/Particles/particles.yaml");

	std::cout << "Finished loading particles data" << std::endl;

	kappa::Mixture mixture(mol, at, "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml");
    kappa::Mixture mixture_pure(mol, "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml");

    std::vector<double> T_vals;
    int i;
    for (i = 0; i < 159; i++) {
        T_vals.push_back(500 + i * 250);
    }
    int x_atom_perc = 20;
    double x_atom = x_atom_perc / 100.0;

    std::vector<arma::vec> mol_ndens;
    mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), mol));
    arma::vec atom_ndens(1);

    std::ofstream outf, outf_pure, outf_xat99, outf_xat01;
    outf.open("c:/Users/st024385/Documents/Code/kappa-tests/bulk_viscosity/bulk_sts_" + mol.name + "_" + at.name + "_xat" + std::to_string(x_atom_perc) + ".txt");
    outf_pure.open("c:/Users/st024385/Documents/Code/kappa-tests/bulk_viscosity/bulk_sts_" + mol.name + "_pure.txt");
    outf_xat99.open("c:/Users/st024385/Documents/Code/kappa-tests/bulk_viscosity/bulk_sts_" + mol.name + "_" + at.name + "_xat99.txt");
    outf_xat01.open("c:/Users/st024385/Documents/Code/kappa-tests/bulk_viscosity/bulk_sts_" + mol.name + "_" + at.name + "_xat01.txt");
    outf << "T; zeta; zeta / eta" << std::endl;
    outf_pure << "T; zeta; zeta / eta" << std::endl;
    outf_xat99 << "T; zeta; zeta / eta" << std::endl;
    outf_xat01 << "T; zeta; zeta / eta" << std::endl;
    double tot_ndens;

    for (auto T : T_vals) {
        tot_ndens =  101325.0 / (K_CONST_K * T);
        mol_ndens[0] = mixture.Boltzmann_distribution(T, (1 - x_atom) * tot_ndens, mol);
        atom_ndens[0] = x_atom * tot_ndens;

        mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);
        outf << T << ";" << mixture.get_bulk_viscosity() << ";" << mixture.get_bulk_viscosity() / mixture.get_shear_viscosity() << std::endl;

        mol_ndens[0] = mixture.Boltzmann_distribution(T, tot_ndens, mol);
        mixture_pure.compute_transport_coefficients(T, mol_ndens);
        outf_pure << T << ";" << mixture_pure.get_bulk_viscosity() << ";" << mixture_pure.get_bulk_viscosity() / mixture_pure.get_shear_viscosity() << std::endl;

        mol_ndens[0] = mixture.Boltzmann_distribution(T, 0.01 * tot_ndens, mol);
        atom_ndens[0] = 0.99 * tot_ndens;
        mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);
        outf_xat99 << T << ";" << mixture.get_bulk_viscosity() << ";" << mixture.get_bulk_viscosity() / mixture.get_shear_viscosity() << std::endl;

        mol_ndens[0] = mixture.Boltzmann_distribution(T, 0.99 * tot_ndens, mol);
        atom_ndens[0] = 0.01 * tot_ndens;
        mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);
        outf_xat01 << T << ";" << mixture.get_bulk_viscosity() << ";" << mixture.get_bulk_viscosity() / mixture.get_shear_viscosity() << std::endl;
    }

    outf.close();
    outf_pure.close();
    outf_xat01.close();
    outf_xat99.close();
    std::string a;
	std::cout << "Enter anything to quit: ";
	std::cin >> a;
    return 0;
}

