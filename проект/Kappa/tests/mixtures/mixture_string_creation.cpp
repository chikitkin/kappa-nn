
#include "kappa.hpp"
#include <iostream>

int main(int argc, char** argv) {
	kappa::Mixture N2N("N2, N", "c:/Users/st024385/Documents/DB/Particles/particles.yaml", "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml");

	std::cout << N2N.get_names();

	std::string a;
	std::cout << "Enter anything to quit: ";
	std::cin >> a;
}

