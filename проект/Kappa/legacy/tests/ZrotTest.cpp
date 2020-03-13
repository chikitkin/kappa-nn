/* 
 * File:   ZrotTest.cpp
 */

#include <cstdlib>
#include <iostream>
#include <fstream>

#include "kappa.hpp"




using namespace std;
using namespace kappa;

int main(int argc, char** argv) {
    cout << "Start Test for Zrot" << endl;
    // ---------------------------------------
    
    int i = 0; // vibrational level
    double Zr, Zr_simple;
		
	std::string molecule_name = "N2";
    Molecule MoleculeTest(molecule_name,"anharmonic",true);
	ApproximationSTS ApproximationTest{};

	bool write_to_file = true;
	ofstream output;
	if (write_to_file) {
		output.open(MoleculeTest.name + ".txt");
	}

	std::vector<double> T_vals = {1000., 10000., 15000., 20000., 25000., 30000., 40000.};

	for (auto T : T_vals) {
		Zr = ApproximationTest.Z_rot(MoleculeTest, T, 0, i, false);
		Zr_simple = ApproximationTest.Z_rot(MoleculeTest, T, 0, i, true);
		std::cout << MoleculeTest.name << ".Z_rot_simplified(" << T << ") = " << Zr_simple << endl;
		std::cout << MoleculeTest.name << ".Z_rot(" << T << ") = " << Zr << endl;
		std::cout << MoleculeTest.name << ".c_rot_simplified(" << T << ") = " << ApproximationTest.c_rot(MoleculeTest, T, 0, i, true) << endl;
	    std::cout << MoleculeTest.name << ".c_rot_Z_simplified(" << T << ") = " << ApproximationTest.c_rot(MoleculeTest, T, 0, i, false) * Zr / Zr_simple << endl;
	    std::cout << MoleculeTest.name << ".c_rot(" << T << ") = " << ApproximationTest.c_rot(MoleculeTest, T, 0, i, false) << endl << endl;
		if (write_to_file) {
			output << MoleculeTest.name << ".Z_rot_simplified(" << T << ") = " << Zr_simple << endl;
			output << MoleculeTest.name << ".Z_rot(" << T << ") = " << Zr << endl;
		    output << MoleculeTest.name << ".c_rot_simplified(" << T << ") = " << ApproximationTest.c_rot(MoleculeTest, T, 0, i, true) << endl;
		    output << MoleculeTest.name << ".c_rot_Z_simplified(" << T << ") = " << ApproximationTest.c_rot(MoleculeTest, T, 0, i, false) * Zr / Zr_simple << endl;
		    output << MoleculeTest.name << ".c_rot(" << T << ") = " << ApproximationTest.c_rot(MoleculeTest, T, 0, i, false) << endl << endl;
		}
	}

	if (write_to_file) {
		output.close();
	}

    std::cout << "Enter anything to quit: ";
    std::string a;
    std::cin >> a;
    return 1; 
}

