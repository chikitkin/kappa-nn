#include <iostream>
#include <fstream>
#include <iomanip> 


#include "kappa.hpp"
using namespace std;
using namespace kappa;

int main_old(int argc, char** argv) {

	std::cout << "Start test: computation of diffusion coefficients" << std::endl;
	std::cout << "Loading particles data" << std::endl;

	Molecule mol("O2", true, false); // "N2", "NO"
	Atom at("O"); //"N"

	cout << "Finished loading particles data" << std::endl;

	cout << "Molecule's names: " << mol.name <<endl;
	cout << "Atom's names: " << at.name << endl;

	Mixture mixture( mol, at);
	Approximation approx{};

	//// set a range for temperature
	//vector<double> T_vals = { 2500.0 , 5000.0, 7500.0, 10000.0, 15000.0, 20000.0, 30000.0 };
	//double p = 101325.;

	//// vector of arma vector for molecular number density
	arma::vec mol_ndens;
	double tot_ndens;
	//ofstream diffus("diffusion.txt");

	//for (int i = 0; i < size(T_vals); i++)
	//{
	//	vector <arma::vec> mol_ndens;
	//	arma::vec atom_ndens(1);
	//	cout << "T= " << T_vals[i] << endl;
	//	diffus << "T= " << T_vals[i] << endl;
	//	tot_ndens = p / (K_CONST_K * T_vals[i]);

	//	atom_ndens[0] =  0.2*tot_ndens;

	//	// Coefficients are stored in the following order : DA2iA2i | DA2iA2k | DA2A | DAA
	//	mol_ndens = { mixture.Boltzmann_distribution(T_vals[i], 0.8*tot_ndens, mol) };
	//	mixture.compute_transport_coefficients(T_vals[i], mol_ndens, atom_ndens, kappa::models_omega::model_omega_esa, 1e-16);
	//	cout << "compute finish" << endl;
	//
	//	diffus << mixture.get_binary_diffusion() << endl << endl;

	//}
	
	// Читаем файл с точками p, T, n_dens_atom в которых нужно посчитать коээффициенты переноса 
	string input_file = "input.txt";
	string output_file = "output.txt";

	ifstream input(input_file);
	ofstream output(output_file);

	output.precision(6);
	output << scientific;
	double p_, T_, n_atom_;
	int number_of_points;
	input >> number_of_points;
	for (int j = 0; j < number_of_points; ++j) {

		input >> p_ >> T_ >> n_atom_;
		tot_ndens = p_ / (K_CONST_K * T_);
		
		vector <arma::vec> mol_ndens;
		arma::vec atom_ndens(1);

		atom_ndens[0] = n_atom_ * tot_ndens;
		mol_ndens = { mixture.Boltzmann_distribution(T_, (1 - n_atom_) * tot_ndens, mol) };
		mixture.compute_transport_coefficients(T_, mol_ndens, atom_ndens, kappa::models_omega::model_omega_esa, 1e-9);
		auto coeffs = mixture.get_binary_diffusion();
		for (int i = 0; i < coeffs.size(); ++i) {
			output << coeffs[i] << " ";
		}
		output << endl;
	}

	return 0;
}