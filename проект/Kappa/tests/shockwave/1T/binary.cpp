#include <iostream>
#include "kappa.hpp"
#include <fstream>
#include <string>
#define _CRT_SECURE_NO_WARNINGS

#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;

typedef std::vector< double > state_type;
typedef runge_kutta_fehlberg78< state_type > error_stepper_type;

class binary_flow {
	kappa::Molecule molecule;
	kappa::Atom atom;
	kappa::Approximation approx;
	kappa::Interaction interaction_mm, interaction_ma;
	arma::mat44 coeff_matrix;
	arma::vec4 rhs;

	double avg_vibr_energy(double T, const arma::vec &vibr_energy) {
    	return arma::dot(vibr_energy, arma::exp(-vibr_energy / (K_CONST_K * T))) / approx.Z_vibr_eq(T, molecule);
	}

	double avg_vibr_energy_sq(double T, const arma::vec &vibr_energy) {
    	return arma::dot(vibr_energy % vibr_energy, arma::exp(-vibr_energy / (K_CONST_K * T))) / approx.Z_vibr_eq(T, molecule);
	}

	double E_tr(double T, double n) {
		return 1.5 * n * K_CONST_K * T;
	}

	double E_rot(double T, double nmol) {
		return approx.avg_rot_energy(T, molecule) * nmol;
	}

	double E_vibr(double T, double nmol) {
		return avg_vibr_energy(T, molecule.vibr_energy[0]) * nmol;
	}

	double E_f(double T, double natom) {
		return atom.formation_energy * natom;
	}

	double c_tr(double n) {
		return 1.5 * K_CONST_K * n;
	}

	double c_rot(double T) {
		return approx.c_rot(T, molecule);
	}

	double c_vibr(double T) {
		double avg_ve = avg_vibr_energy(T, molecule.vibr_energy[0]);
		return (avg_vibr_energy_sq(T, molecule.vibr_energy[0]) - avg_ve * avg_ve) / (K_CONST_K * T * T * molecule.mass);
	}

	double R_diss(double T, double nmol, double natom) {
		double k_diss_mol_mol = approx.k_diss(T, molecule, interaction_mm, 0, kappa::models_k_diss::model_k_diss_arrh_scanlon);
		double k_diss_mol_atom = approx.k_diss(T, molecule, interaction_ma, 0, kappa::models_k_diss::model_k_diss_arrh_scanlon);
		double K_rec_diss = pow(molecule.mass / (atom.mass * atom.mass), 1.5) * K_CONST_H * K_CONST_H * K_CONST_H
		                  * pow(2 * K_CONST_PI * K_CONST_K * T, -1.5)
		                  * approx.Z_vibr_eq(T, molecule) * approx.Z_rot(T, molecule) * exp(molecule.diss_energy[0] / (K_CONST_K * T));
		return nmol * (natom * natom * K_rec_diss - nmol) * k_diss_mol_mol + natom * (natom * natom * K_rec_diss - nmol) * k_diss_mol_atom; 
	}

public:

	binary_flow(std::string molecule_name, std::string atom_name, bool anharmonic_spectrum=true) : molecule(molecule_name, anharmonic_spectrum), atom(atom_name), interaction_mm(molecule, molecule),
	                                                                                               interaction_ma(molecule, atom), approx() {
		rhs.zeros();
	};

	double lambda, v_init, T_init, nmol_init;;

	void operator() (const state_type &x, state_type &dxdt, const double t)
	{
		// dimensionless variables:
		// v = x[0], T=x[1], n_molecule=x[2], n_atom=x[3]
		
		coeff_matrix.zeros();
		double n = x[2] + x[3]; // dimensionless total numeric density
		double rho = x[2] + atom.mass * x[3] / molecule.mass; // dimensionless density

		// some dimensional variables
		double p = n * nmol_init * K_CONST_K * x[1]; // pressure
		double T = x[1] * T_init; // temperature
		double nmol = x[2] * nmol_init;

		// calculate energies
		double E_rot_val = E_rot(T, nmol);
		double E_vibr_val = E_vibr(T, nmol);
		double E = E_tr(T, n * nmol_init) + E_rot_val + E_vibr_val + E_f(T, x[3] * nmol_init);

		// calculate specific heats
		double dEdnatom = 1.5 * K_CONST_K * T + atom.formation_energy;
		double dEdnmol = 1.5 * K_CONST_K * T  + E_rot_val / nmol + E_vibr_val / nmol;
		double dEdT = c_tr(n) + c_rot(T) * molecule.mass * nmol + c_vibr(T) * molecule.mass * nmol;

		// calculate relaxation terms
		double R = R_diss(T, nmol, x[3] * nmol_init) * (lambda / (nmol_init * v_init));

		coeff_matrix.at(0,0) = rho * x[0] * molecule.mass * v_init * v_init / (K_CONST_K * T_init);
		coeff_matrix.at(0,1) = n;
		coeff_matrix.at(0,2) = x[1];
		coeff_matrix.at(0,3) = x[1];

		coeff_matrix.at(1,0) = (E + p) / (nmol_init * K_CONST_K * T_init);
		coeff_matrix.at(1,1) = x[0] * dEdT / (nmol_init * K_CONST_K);
		coeff_matrix.at(1,2) = x[0] * dEdnmol / (K_CONST_K * T_init);
		coeff_matrix.at(1,3) = x[0] * dEdnatom / (K_CONST_K * T_init);

		coeff_matrix.at(2,0) = x[2];
		coeff_matrix.at(2,2) = x[0];

		coeff_matrix.at(3,0) = x[3];
		coeff_matrix.at(3,3) = x[0];

		rhs[2] = R;
		rhs[3] = -2 * R;
		
		// get the derivatives vector
		arma::vec derivatives = arma::solve(coeff_matrix, rhs);
		
		dxdt[0] = derivatives[0];
		dxdt[1] = derivatives[1];
		dxdt[2] = derivatives[2];
		dxdt[3] = derivatives[3];
	}

	state_type rg_init(double n1, double v1, double T1) {
		state_type res(4); // return dimensionless vector of v, T, n_molecules, n_atom
		// double gamma = 1.4;

		// double Rbar = K_CONST_R / (molecule.mass * K_CONST_NA);

		// double rho0 = p0 * molecule.mass / (K_CONST_K * T0);

		// double p1 = p0 * ((2 * gamma) * M0 * M0 - (gamma - 1)) / (gamma + 1);
		// double rho1 = rho0 * ((gamma + 1) * M0 * M0 / ((gamma - 1) * M0 * M0 + 2));
		// double T1 = (p1 / p0) * (rho0 / rho1) * T0;

		// double M1 = sqrt((1 + 0.5 * (gamma - 1) * M0 * M0) / (gamma * M0 * M0 - 0.5 * (gamma - 1)));

		v_init = v1;
		T_init = T1;
		nmol_init = n1;
		lambda = 1 / (sqrt(2.) * nmol_init * molecule.diameter);

		res[0] = 1;
		res[1] = 1;
		res[2] = 1;
		res[3] = 0.00000;
		return res;
	}

	std::string get_name() {
		return molecule.name + "_" + atom.name;
	}
};

struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};

int main(void) {
	int flow_type;

	std::cout << "Enter 1 for N2/N flow, any other number for O2/O flow: ";
	std::cin >> flow_type;

	binary_flow flow("N2", "N");

	if (flow_type != 1) {
		binary_flow fl2("O2", "O");
		flow = fl2;
	}

	double n1, v1, T1;

	std::cout << "Enter numeric density behind shockwave: ";
	std::cin >> n1;
	std::cout << "Enter velocity behind shockwave: ";
	std::cin >> v1;
	std::cout << "Enter temperature behind shockwave: ";
	std::cin >> T1;

	state_type init_cond = flow.rg_init(n1, v1, T1);
	std::ofstream outf(flow.get_name() + "_sw.txt");
	// std::cout << init_cond[0] << " " << init_cond[1] << " " << init_cond[2] << std::endl;

	std::vector<state_type> macro_vars;
	std::vector<double> x;
	double max_x = 0.5;
	double x_step = 0.001;

	std::cout << "Enter maximum x (in meters): ";
	std::cin >> max_x;
	std::cout << "Enter step size (in meters):";
	std::cin >> x_step;
	std::cout << std::endl;

	size_t steps = integrate_adaptive(make_controlled(1.0e-10, 1.0e-6, error_stepper_type()), flow, init_cond, 0.0, max_x / flow.lambda, x_step / flow.lambda, push_back_state_and_time(macro_vars, x));
	outf << "x;v;T;n_mol;n_atom" << std::endl;
	for(size_t i=0; i<=steps; i++ )
	{
		outf << x[i] * flow.lambda << ";" << macro_vars[i][0] * flow.v_init << ";" << macro_vars[i][1] * flow.T_init << ";" << macro_vars[i][2] * flow.nmol_init << ";" << macro_vars[i][3] * flow.nmol_init << std::endl;
	    std::cout << "x=" << x[i] * flow.lambda << ", v=" << macro_vars[i][0] * flow.v_init;
	    std::cout << ", T=" << macro_vars[i][1] * flow.T_init << ", n_mol=" << macro_vars[i][2]  * flow.nmol_init << ", n_atom=" << macro_vars[i][3] * flow.nmol_init << std::endl;
	}



	std::string a;
	std::cin >> a;
	return 1;
}