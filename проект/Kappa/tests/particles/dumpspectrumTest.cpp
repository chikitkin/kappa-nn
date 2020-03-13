/* 
 * Dump vibrational spectrum
 */

#include <iostream>
#include "kappa.hpp"
#include <fstream>

using namespace std;
using namespace kappa;


int main(int argc, char** argv) {
    std::cout << "Start Test for dump_spectrum" << endl;

    Molecule N2("N2", true, true, "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
	Molecule NO("NO", true, true, "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
	Molecule O2("O2", true, true, "c:/Users/st024385/Documents/DB/Particles/particles.yaml");

	int i;

	ofstream f;
	f.open("c:/Users/st024385/Documents/Code/kappa-tests/N2_spectrum.csv");
	for (i=0; i<N2.num_vibr_levels[0]; i++) {
		f << N2.vibr_energy[0][i] / K_CONST_EV;
		f << ",";
	}
	f.close();

	f.open("c:/Users/st024385/Documents/Code/kappa-tests/NO_spectrum.csv");
	for (i = 0; i<NO.num_vibr_levels[0]; i++) {
		f << NO.vibr_energy[0][i] / K_CONST_EV;
		f << ",";
	}
	f.close();


	f.open("c:/Users/st024385/Documents/Code/kappa-tests/O2_spectrum.csv");
	for (i = 0; i<O2.num_vibr_levels[0]; i++) {
		f << O2.vibr_energy[0][i] / K_CONST_EV;
		f << ",";
	}
	f.close();

    string a;
    cout << "Enter anything to quit: ";
    cin >> a;
    return 0;
}

