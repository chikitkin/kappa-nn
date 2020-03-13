#include <iostream>
#include <fstream>
#include <iomanip> 


#include "kappa.hpp"
using namespace std;
using namespace kappa;

int main(int argc, char** argv) {

	std::cout << "Start test: computation of diffusion coefficients" << std::endl;
	std::cout << "Loading particles data" << std::endl;

	Molecule mol_O2("O2", true, false);
	Molecule mol_N2("N2", true, false);
	Molecule mol_NO("NO", true, false);
	Atom at_O("O");
	Atom at_N("N");

	//Mixture mixture( mol, at);
	Mixture air({ mol_O2, mol_N2, mol_NO }, { at_O,at_N });
	Approximation approx{};
	cout << "Finished loading particles data" << std::endl;


	// set a range for temperature
	vector<double> T_vals = { 2500.0 }; //, 5000.0, 7500.0, 10000.0, 15000.0, 20000.0, 30000.0};
	double p = 101325.;
	double x_N2 = 0.3;
	double x_O2 = 0.3;
	double x_NO = 0.3;
	double x_N = 0.05;
	double x_O = 0.05;

	// vector of arma vector for molecular number density
	double tot_ndens;
	ofstream diffus("diffusion.txt");
	for (int i = 0; i < size(T_vals); i++)
	{
		cout << "T= " << T_vals[i] << endl;
		diffus << "T= " << T_vals[i] << endl;
		tot_ndens = p / (K_CONST_K * T_vals[i]);

		arma::vec atom_ndens(2);		
		std::vector<arma::vec> mol_ndens;
		mol_ndens.push_back(air.Boltzmann_distribution(T_vals[i], x_N2 * tot_ndens, mol_O2));
		mol_ndens.push_back(air.Boltzmann_distribution(T_vals[i], x_O2 * tot_ndens, mol_N2));
		mol_ndens.push_back(air.Boltzmann_distribution(T_vals[i], x_NO * tot_ndens, mol_NO));
		// atoms number density

		atom_ndens[0] =  x_N * tot_ndens;
		atom_ndens[1] = x_O * tot_ndens;		

		air.compute_transport_coefficients(T_vals[i], mol_ndens, atom_ndens, kappa::models_omega::model_omega_esa, 1e-09);
		cout << "compute transport coefficient  finish" << endl;
	
		diffus << air.get_diffusion() << endl << endl;

	}

	string x;
	cout << "Enter anything to quit: ";
	cin >> x;
	return 0;
}