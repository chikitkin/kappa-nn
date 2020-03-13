#include <iostream>
#include "kappa.hpp"
#include <fstream>
#include <chrono>


int main(int argc, char** argv) {
    std::cout << "Start test: computation of shear viscosity" << std::endl;

    std::vector<kappa::Molecule> molecules;
    std::vector<kappa::Atom> atoms;

    std::cout << "Loading particles data" << std::endl;

    kappa::Molecule N2("N2", true, true, "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
    kappa::Molecule O2("O2", true, true, "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
    kappa::Molecule NO("NO", true, true, "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
	kappa::Atom N("N", "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
    kappa::Atom O("O", "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
    molecules.push_back(N2); 
    molecules.push_back(O2);
    molecules.push_back(NO);
    atoms.push_back(N);
    atoms.push_back(O);

	std::cout << "Finished loading particles data" << std::endl;

	kappa::Mixture mixture(molecules, atoms, "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml");

    std::vector<double> T_vals = { 500., 1000., 5000., 10000., 20000., 40000. };
    double x_N2 = 0.2, x_O2 = 0.2, x_NO = 0.2, x_N = 0.2, x_O = 0.2;

    std::vector<arma::vec> mol_ndens;
	arma::vec atom_ndens(2);

	mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), N2));
    mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), O2));
    mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), NO));
    
    double tot_ndens;

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    for (auto T : T_vals) {
        tot_ndens =  101325.0 / (K_CONST_K * T);
        mol_ndens[0] = mixture.Boltzmann_distribution(T, x_N2 * tot_ndens, N2);
        mol_ndens[1] = mixture.Boltzmann_distribution(T, x_O2 * tot_ndens, O2);
        mol_ndens[2] = mixture.Boltzmann_distribution(T, x_NO * tot_ndens, NO);
        atom_ndens[0] = x_N * tot_ndens;
        atom_ndens[1] = x_O * tot_ndens;
        mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);
        std::cout << T << " " << mixture.get_shear_viscosity() << " " << mixture.get_bulk_viscosity() << " " << mixture.get_thermal_conductivity() << std::endl;
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    std::cout << duration << std::endl;
    // outf.close();
    std::string a;
	std::cout << "Enter anything to quit: ";
	std::cin >> a;
    return 0;
}

