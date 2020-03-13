#include <iostream>
#include <fstream>
#include <iomanip> 

#include "kappa.hpp"
using namespace std;
using namespace kappa;

int main1(int argc, char** argv) {

	std::cout << "Start test: computation of diffusion coefficients" << std::endl;

	std::cout << "Loading particles data" << std::endl;

	// N2 molecule (non-rigid model)
	Molecule mol("N2", true, false);

	// N atom
	Atom at("N");

	cout << "Finished loading particles data" << std::endl;

	cout << "Molecule's name " << mol.name << std::endl;
	cout << "Atom's name " << at.name << std::endl;
	cout << "Molecule vibrational levels " << mol.num_vibr_levels[0] << std::endl;

	// N2/N binary mixture creation
	Mixture mixture(mol, at);

	// some check print
	cout << "particles: " << mixture.get_n_particles() << std::endl;
	cout << "names: " << mixture.get_names() << std::endl;

	// set a range for temperature
	//vector<double> T_vals = { 50000.0 };
	vector<double> T_vals = {2500.0, 5000.0, 20000.0, 50000.0};
	// vector<double> T_vals = {15000.0};
	// vector<double> T_vals = {30000.0};
	// vector<double> T_vals = {1.1986210569919167e+04};

	double p = 101325.;

	// for (i = 0; i < 80; i++) { // assume a max num. of vibr. levels a priori
	//   T_vals.push_back(500 + i * 500);
	// }

	// set an arbitrary atom mass fraction 
	double x_atom_perc = 20.;
	double x_atom = x_atom_perc / 100.;

	// arma vector for atom number density
	arma::vec atom_ndens(1);

	// vector of arma vector for molecular number density
	vector<arma::vec> mol_ndens;

	double tot_ndens;
	tot_ndens = p / (K_CONST_K * T_vals[0]);
	cout << " tot_ndens " << tot_ndens << std::endl;

	mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], (1. - x_atom) * tot_ndens, mol));

	// atom number density
	atom_ndens = x_atom * tot_ndens;

	for (int at = 0; at < atom_ndens.size(); at++) {
		cout << " atom_ndens " << atom_ndens.at(at) << std::endl;
	}
	ofstream diff("diffusion.txt");
	for (int i = 0; i < size(T_vals); i++)
	{
		diff << "T= " << T_vals[i] << endl;
		mixture.compute_transport_coefficients(T_vals[i], mol_ndens, atom_ndens);
		cout << "compute finish" << endl;

		diff << mixture.get_diffusion() << endl << endl;
	}

	string x;
	cout << "Enter anything to quit: ";
	cin >> x;
	return 0;
}