/* 
 * File:   ZintTest.cpp
 */

#include <cstdlib>
#include <iostream>

#include "kappa.hpp"




using namespace std;
using namespace kappa;

int main(int argc, char** argv) {
    cout << "Start Test for Zint" << endl;
    // ---------------------------------------
    
    double T;
	double DE = convert_cm_to_Joule(1000); // value of Delta_E to cut off electron levels of atoms

    Molecule MoleculeTest("N2",true,false,"particles.yaml");  // molecule name, type of vibrational spectrum, whether the rigid rotator model is used
    Atom AtomTest("N", "particles.yaml");  // atom name
	ApproximationOneT ApproximationOneTTest{};

    T = 1000;
    std::cout << MoleculeTest.name << ".Z_int(" << T << ") = " << ApproximationOneTTest.Z_int(T, MoleculeTest) << endl;
    std::cout << MoleculeTest.name << ".E_int(" << T << ")/RT = " << ApproximationOneTTest.avg_energy(T, MoleculeTest) / (K_CONST_K * T) << endl;
	std::cout << MoleculeTest.name << ".c_int(" << T << ")/R = " << ApproximationOneTTest.c_int(T, MoleculeTest) * MoleculeTest.mass / K_CONST_K << endl << endl;

	T = 10000;
	std::cout << MoleculeTest.name << ".Z_int(" << T << ") = " << ApproximationOneTTest.Z_int(T, MoleculeTest) << endl;
	std::cout << MoleculeTest.name << ".E_int(" << T << ")/RT = " << ApproximationOneTTest.avg_energy(T, MoleculeTest) / (K_CONST_K * T) << endl;
	std::cout << MoleculeTest.name << ".c_int(" << T << ")/R = " << ApproximationOneTTest.c_int(T, MoleculeTest) * MoleculeTest.mass / K_CONST_K << endl << endl;

	T = 20000;
	std::cout << MoleculeTest.name << ".Z_int(" << T << ") = " << ApproximationOneTTest.Z_int(T, MoleculeTest) << endl;
	std::cout << MoleculeTest.name << ".E_int(" << T << ")/RT = " << ApproximationOneTTest.avg_energy(T, MoleculeTest) / (K_CONST_K * T) << endl;
	std::cout << MoleculeTest.name << ".c_int(" << T << ")/R = " << ApproximationOneTTest.c_int(T, MoleculeTest) * MoleculeTest.mass / K_CONST_K << endl << endl;

    T = 30000;
	std::cout << MoleculeTest.name << ".Z_int(" << T << ") = " << ApproximationOneTTest.Z_int(T, MoleculeTest) << endl;
	std::cout << MoleculeTest.name << ".E_int(" << T << ")/RT = " << ApproximationOneTTest.avg_energy(T, MoleculeTest) / (K_CONST_K * T) << endl;
	std::cout << MoleculeTest.name << ".c_int(" << T << ")/R = " << ApproximationOneTTest.c_int(T, MoleculeTest) * MoleculeTest.mass / K_CONST_K << endl << endl;

    T = 1000;
	std::cout << AtomTest.name << ".Z_int(" << T << ") = " << ApproximationOneTTest.Z_int(T, AtomTest, DE) << endl << endl;
    std::cout << AtomTest.name << ".E_int(" << T << ")/RT = " << ApproximationOneTTest.avg_energy(T, AtomTest, DE) / (K_CONST_K * T) << endl;
    std::cout << AtomTest.name << ".c_int(" << T << ")/R = " << ApproximationOneTTest.c_int(T, AtomTest, DE) * AtomTest.mass / K_CONST_K << endl << endl;

    T = 10000;
	std::cout << AtomTest.name << ".Z_int(" << T << ") = " << ApproximationOneTTest.Z_int(T, AtomTest, DE) << endl << endl;
	std::cout << AtomTest.name << ".E_int(" << T << ")/RT = " << ApproximationOneTTest.avg_energy(T, AtomTest, DE) / (K_CONST_K * T) << endl;
	std::cout << AtomTest.name << ".c_int(" << T << ")/R = " << ApproximationOneTTest.c_int(T, AtomTest, DE) * AtomTest.mass / K_CONST_K << endl << endl;

    T = 20000;
	std::cout << AtomTest.name << ".Z_int(" << T << ") = " << ApproximationOneTTest.Z_int(T, AtomTest, DE) << endl << endl;
	std::cout << AtomTest.name << ".E_int(" << T << ")/RT = " << ApproximationOneTTest.avg_energy(T, AtomTest, DE) / (K_CONST_K * T) << endl;
	std::cout << AtomTest.name << ".c_int(" << T << ")/R = " << ApproximationOneTTest.c_int(T, AtomTest, DE) * AtomTest.mass / K_CONST_K << endl << endl;

    T = 30000;
	std::cout << AtomTest.name << ".Z_int(" << T << ") = " << ApproximationOneTTest.Z_int(T, AtomTest, DE) << endl << endl;
	std::cout << AtomTest.name << ".E_int(" << T << ")/RT = " << ApproximationOneTTest.avg_energy(T, AtomTest, DE) / (K_CONST_K * T) << endl;
	std::cout << AtomTest.name << ".c_int(" << T << ")/R = " << ApproximationOneTTest.c_int(T, AtomTest, DE) * AtomTest.mass / K_CONST_K << endl << endl;

	std::cout << "Enter anything to quit: ";
    std::string a;
    std::cin >> a;
    return 1; 
}

