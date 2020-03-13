/* 
 */

#include <iostream>
#include "kappa.hpp"
#include <fstream>


int main(int argc, char** argv) {
    std::cout << "Start Test state-to-state mixture" << std::endl;

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

    for (auto at: atoms) {
        for (auto mo: molecules) {
            std::cout << at.name << "+" << mo.name << ", interaction=" << mixture.interaction(at, mo).particle1_name << "+" << mixture.interaction(at, mo).particle2_name << std::endl;
        }
    }

    for (auto at: atoms) {
        for (auto at2: atoms) {
            std::cout << at.name << "+" << at2.name << ", interaction=" << mixture.interaction(at, at2).particle1_name << " " << mixture.interaction(at, at2).particle2_name << std::endl;
        }
    }

    for (auto mo: molecules) {
        for (auto mo2: molecules) {
            std::cout << mo.name << "+" << mo2.name << ", interaction=" << mixture.interaction(mo, mo2).particle1_name << " " << mixture.interaction(mo, mo2).particle2_name << std::endl;
        }
    }

  //   std::vector<arma::vec> mol_ndens;
  //   mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), N2));


  //   for (auto T : T_vals) {
  //       mol_ndens[0] = mixture.Boltzmann_distribution(T, 101325.0 / (K_CONST_K * T), N2);
		// std::cout << T << " " << mol_ndens[0][0] / (101325.0 / (K_CONST_K * T));
  //   }

    std::string a;
	std::cout << "Enter anything to quit: ";
	std::cin >> a;
    return 0;
}

