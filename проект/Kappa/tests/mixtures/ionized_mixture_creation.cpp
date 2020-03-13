#include <iostream>
#include "kappa.hpp"
#include <fstream>


int main(int argc, char** argv) {
    std::cout << "Loading particles data" << std::endl;

    kappa::Molecule mol("N2", true, true, "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
    kappa::Molecule mol_ion("N2+", true, true, "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
	kappa::Atom at("N", "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
    kappa::Atom at_ion("N+", "c:/Users/st024385/Documents/DB/Particles/particles.yaml");

	std::cout << "Finished loading particles data" << std::endl;

    kappa::Mixture mixture1(mol_ion, "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml", "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
    std::cout << "Finished loading N2+,e- Mixture" << std::endl;

    kappa::Mixture mixture2(at_ion, "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml", "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
    std::cout << "Finished loading N+,e- Mixture" << std::endl;

	kappa::Mixture mixture3(mol_ion, at, "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml", "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
    std::cout << "Finished loading N2+,N,e- Mixture" << std::endl;

    kappa::Mixture mixture4("N2+,N2", "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml", "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
    std::cout << "Finished loading N2+,N2,e- Mixture" << std::endl;

    kappa::Mixture mixture5("N2,N,N+", "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml", "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
    std::cout << "Finished loading N2,N,N+,e- Mixture" << std::endl;

    std::cout << mixture1.get_names() << std::endl;
    std::cout << mixture2.get_names() << std::endl;
    std::cout << mixture3.get_names() << std::endl;
    std::cout << mixture4.get_names() << std::endl;
    std::cout << mixture5.get_names() << std::endl << std::endl;

    kappa::Interaction inter(mol, mol, "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml");


    inter = mixture5.interaction(mol, mol); // N2 + N2
    std::cout << inter.particle1_name << " " << inter.particle2_name << std::endl;
    inter = mixture5.interaction(mol, at); // N2 + N
    std::cout << inter.particle1_name << " " << inter.particle2_name << std::endl;
    inter = mixture5.interaction(mol, at_ion); // N2 + N+
    std::cout << inter.particle1_name << " " << inter.particle2_name << std::endl << std::endl;
    inter = mixture5.interaction(at, at_ion); // N + N+
    std::cout << inter.particle1_name << " " << inter.particle2_name << std::endl << std::endl;

    std::string a;
	std::cout << "Enter anything to quit: ";
	std::cin >> a;
    return 0;
}

