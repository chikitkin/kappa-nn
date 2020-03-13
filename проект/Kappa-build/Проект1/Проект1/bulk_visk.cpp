#include <iostream>
#include "kappa.hpp"
#include <fstream>
#include <cmath>
#include <ctime>
#include <numeric>

using namespace std;
using namespace kappa;

int main()
{
	cout << "start test" << endl;
	Approximation approx{};
	Molecule mol1("N2", true, false);
	Molecule mol2("N2", true, false);
	Atom at1("N");
	//Mixture mix(const std::vector<kappa::Molecule> &mol1,
		//const std::vector<kappa::Atom> &at1, const std::string &interactions_filename = "interaction.yaml", 
		//const std::string &particles_filename = "particles.yaml");
	
	//Interaction inter(mol2, mol2);
	//ofstream Omega("Omega.txt");
	//models_omega model = kappa::models_omega::model_omega_esa;

	//Omega << "NO-NO" << endl;
	double p = 101325;
	//Omega << inter.collision_mass << endl;
	//Omega << K_CONST_K << endl;
	//Omega << mol1.mass << endl;
	//Omega << "T       |     Omega^(1,1) |     Omega^(2,2)    |    Omega^(1,2)     |    Omega^(2,3)   |    Omega^(1,3)" << endl;
	//Omega << approx.omega_integral(500.0, inter, 1, 1, model,false)<<endl;

	//std::string get_names();
	Mixture mix2((mol1,mol1) );
	double n = p / (1000 * K_CONST_K);
	arma::vec nci = approx.Boltzmann_distribution(1000, n, mol1);
	
	cout << nci<< endl << endl;
	int num = nci.size();
	cout << num << endl << endl;

	arma:: vec nc(1) ;
	nc[0] = 0;
	nc[1]=0;
	cout << nc[0] << endl;

	for (int i = 0; i < num ; i++)
	{
		nc[0] = nc[0] + nci[i];
	
	}
	nc[1] = nc[0];
	cout << nc << endl << endl;
	
	cout << mix2.get_n_particles()<<endl;
	cout << mix2.get_names()<<endl<<endl;

	mix2.compute_transport_coefficients(1000, nc);
	cout << "thermal conductivity" << endl;
	cout<<mix2.get_thermal_conductivity()<<endl;
	cout << "x" << endl;
	int x;
	cin >> x;
	return 0;
}