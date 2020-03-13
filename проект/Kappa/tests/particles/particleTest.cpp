/* 
 * File:   particleTest.cpp
 * Author: aspera
 *
 * Created on 4 марта 2016 г., 12:36
 */

#include <cstdlib>
#include <iostream>

#include "kappa.hpp"




using namespace std;
using namespace kappa;

int main(int argc, char** argv) {
    cout << "Start Test for Particle classes" << endl;
    // ---------------------------------------
    

    try {
        Particle Ar_badfilename("Ar","articles.yaml");  // test that trying to read a non-existent database file throws the correct error
    }
    catch (const UnopenedFileException& e) {
        std::cout << e.what() << endl;
    }

	try {
		Particle Bad("B");  // test that trying to load particles that are not in the database throw the correct error
	}
	catch (const DataNotFoundException& e) {
		std::cout << e.what() << endl;
	}

    Particle Ar("Ar", "particles.yaml");
    Particle C("C");
    Particle e("e-");
    Particle N("N");
    Particle Nplus("N+");
    Particle O("O");
    Particle Oplus("O+");
    Particle Ominus("O-");
    
    // --------------------------------------- 
    Atom Ar_atom("Ar","particles.yaml");
    Atom C_atom("C");
    Atom N_atom("N");
    Atom Nplus_atom("N+");
    Atom O_atom("O");
    Atom Oplus_atom("O+");
    Atom Ominus_atom("O-");
    
    // --------------------------------------- 
    Molecule C2("C2",true,false,"particles.yaml");
    Molecule CO("CO");
    Molecule N2("N2");
    Molecule N2plus("N2+");
    Molecule NO("NO");
    Molecule O2("O2");
    Molecule O2plus("O2+");

    std::cout << "Number of electron levels in CO: " << CO.num_electron_levels << endl;  // should be 28
    // std::cout << "Vibrational spectrum of CO: " << CO.vibr_spectrum << endl;  // should be "anharmonic"
    std::cout << "Rigid rotator model used for CO molecule: " << CO.rigid_rotator << endl; // should be true (1)
    std::cout << "N2 number of vibrational levels in ground state: " << N2.num_vibr_levels[0] << endl;  // should be 48


    try {
        Molecule CObad2("CO",true,false,"articles.yaml");  // test that trying to read a non-existent database file throws the correct error
    }
    catch (const UnopenedFileException& e) {
        std::cout << e.what() << endl;
    }

    Molecule C2_h("C2", false);
    Molecule CO_h("CO", false);
    Molecule N2_h("N2", false);
    Molecule N2plus_h("N2+", false);
    Molecule NO_h("NO", false);
    Molecule O2_h("O2", false);
    Molecule O2plus_h("O2+", false);
    std::cout << "N2 number of vibrational levels in ground state: " << N2_h.num_vibr_levels[0] << endl;  // should be 33
    std::cout << "O2 number of vibrational levels in ground state: " << O2_h.num_vibr_levels[0] << endl;
    
	std::string a;
	std::cin >> a;

	Atom wrong_at("N2");

	std::cout << "loaded N2 as atom";

	// Molecule wrong_mol("N");
	// std::cout << "loaded N as molecule";

    return 0;
}

