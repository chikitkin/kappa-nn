#include <iostream>
#include "kappa.hpp"
#include <fstream>
#include <string>
#define _CRT_SECURE_NO_WARNINGS

#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;

typedef std::vector< double > state_type;
typedef runge_kutta_fehlberg78< state_type > error_stepper_type;

//const
const double k = 1.3806488E-23;
const double h = 6.62606957E-34;
const double c = 2.9979E8;
const double Na = 6.0221E23;
const double pi = 3.14;


const double e100 = 2.757135054E-20;
const double e010 = 1.32493446E-20;
const double e001 = 4.66607366E-20;

const int om_CO = 216981;
const double I_CO_2 = 7.5E-46;
const double Be_CO = 1.93128E2;

const double m_CO_2 = 7.3064E-26;
const double m_CO = 4.6495E-26;
const double m_O = 2.6569E-26;

const double e_CO_2 =-3.955E5/Na;
const double e_CO = -5.209E5 / Na;
const double e_O = 2.458E5/ Na;

//const Arrenius
const double A_f1 = 6.9E15 / Na;
const double A_f2 = 1.38E16 / Na;
const double b_f = -1.5;
const int E_f1 = 63275;




class binary_flow {
	arma::mat55 coeff_matrix;
	arma::vec5 rhs;

	//stat_sum
	
	double Zvibr_CO2(double T) {
		int i1, i2, i3;
		double S = 0;
		for (i1 = 0; i1 < 65; i1++)
		{
			for (i2 = 0; i2 < 36; i2++)
			{
				for (i3 = 0; i3 < 20; i3++)
				{
					S += (i2 + 1)*exp(-i1*e100 / k / T)*exp(-i2*e010/k/T)*exp(-i3*e001/k/T);
				}
			}
		}
		return S;
	}

	double Zvibr_CO(double T)
	{
		return 1 / (1 - exp(-h*c*om_CO / (k*T)));
	}

	double Zrot_CO_2(double T)
	{
		return 8 * pow(pi, 2)*I_CO_2*k*T / (pow(h,2));
	}

	double Zrot_CO(double T)
	{
		return k*T / (Be_CO*h*c);
	}


	//energy

	double E_tr(double T, double n) {
		return 1.5 * n *k* T;
	}

	double E_rot(double T, double n_CO2, double n_CO) {
		return (n_CO2+n_CO)*k*T;
	}

	double E_vibr_CO2(double T, double n_CO2) {
		int j1, j2, j3;
		double D = 0;
		for (j1 = 0; j1 < 65; j1++)
		{
			for (j2 = 0; j2 < 36; j2++)
			{
				for (j3 = 0; j3 < 20; j3++)
				{
					D += (j2 + 1)*(j1*e100 + j2*e010 + j3*e001)*exp(-(j1*e100 + j2*e010 + j3*e001) / k / T);
				}
			}
		}
		return D*n_CO2/Zvibr_CO2(T);
	}

	double E_vibr_CO(double T, double n_CO)
	{
		return n_CO*h*c*om_CO/(exp(h*c*om_CO/k/T-1));
	}

	double E_f( double n_CO2,double n_CO, double n_O ) {
		return n_CO2*e_CO_2+n_CO*e_CO+n_O*e_O;
	}

	//help_f
	double dUCOdT(double T)
	{
		return k*pow((h*c*om_CO/(k*T)),2) * exp(h*c*om_CO/(k*T)) / (pow((exp(h*c*om_CO / (k*T)) - 1),2));
	}

	double dUdnCO2(double T)
	{
		double P = 0;
		int k1, k2, k3;
		for (k1 = 0; k1 < 65; k1++)
		{
			for (k2 = 0; k2 < 36; k2++)
			{
				for (k3 = 0; k3 < 20; k3++) 
				{
					
					P +=  ((k2 + 1)*(k1*e100 + k2*e010 + k3*e001)*exp(-(2 * k1 + k2)*e010 / (k*T)))*exp(-k3*e001 / (k*T));
				}

			}
		}
		return 5 * k*T / 2 + e_CO_2 + P / Zvibr_CO2(T);
	}

	double dZdT(double T)
	{
		double S = 0;
		int i1, i2, i3;
		for (i1 = 0; i1 < 65; i1++)
		{
			for (i2 = 0; i2 < 36; i2++)
			{
				for (i3 = 0; i3 < 20; i3++)
				{
					S += (i1*e100*exp(-(i1*e100) / (k*T))) / (k*pow(T,2))*((i2 + 1)*exp(-(i2*e010) / (k*T)))*(exp(-(i3*e001) / (k*T))) 
						+ (exp(-(i1*e100) / (k*T)))*(exp(-(i3*e001) / (k*T)))*((i2 + 1)*i2*e010*exp(-i2*e010 / (k*T)) / (k*pow(T,2))) 
						+ (exp(-(i1*e100) / (k*T)))*((i2 + 1)*exp(-(i2*e010) / (k*T)))*(i3*e001*exp(-i3*e001 / (k*T)) / (k*pow(T,2)));
				}
			}
		}
	}


	double dUdT(double T)
	{
		double S = 0;
		int i1, i2, i3;
		for (i1 = 0; i1 < 65; i1++)
		{
			for (i2 = 0; i2 < 36; i2++)
			{
				for (i3 = 0; i3 < 20; i3++)
				{
					S += (i2 + 1)*(i1*e100 + i2*e010 + i3*e001)*exp(-((2 * i1 + i2)*e010) / (k*T))*exp(-i3*e001 / (k*T));
				}
			}
		}
		double F = 0;
		for (i1 = 0; i1 < 65; i1++)
		{
			for (i2 = 0; i2 < 36; i2++)
			{
				for (i3 = 0; i3 < 20; i3++)
				{
					F+= (i2 + 1)*(i1*e100 + i2*e010 + i3*e001)*exp(-((2 * i1 + i2)*e010) / (k*T))*exp(-i3*e001 / (k*T))*Zvibr_CO2(T)*((2 * i1 + i2)*e010 + i3*e001) / (k*pow(T,2));
				}
			}
		}
		return (F - S*dZdT(T)) / (pow(Zvibr_CO2(T),2));
	}

	//hearts

	double dE_dT(double T, double n_CO2, double n_CO, double n_O)
	{
		return 3 *(n_CO2+n_CO+n_O)*k / 2 + (n_CO2+n_CO)*k + n_CO*dUCOdT(T) + n_CO2*dUdnCO2(T);
	}
	double dE_dn_CO2(double T)
	{
		return 3 * k*T / 2 + k*T + dUdnCO2(T);
			
	}

	double dE_dn_CO (double T)
	{ return 3 * k*T / 2 + k*T + e_CO + h*c*om_CO / (exp(h*c*om_CO / (k*T)) - 1); }

	double dE_dn_O(double T)
	{
		return 3 * k*T / 2 + e_O;
	}

	double R(double T, double n_CO2, double n_CO, double n_O)
	{
		double F;
		F = exp(-(e_CO_2 - e_CO - e_O) / (k*T))*pow(m_CO_2 / (m_O*m_CO),1.5) *pow(h, 3) * pow(2 * pi*k*T,-1.5) *Zrot_CO_2(T)*Zvibr_CO2(T) / (Zrot_CO(T)*Zvibr_CO(T));
		return (F*n_CO*n_O - n_CO2)*(A_f1*pow(T,b_f)*exp(-E_f1 / T)*n_CO2 + (A_f1*pow(T,b_f)*exp(-E_f1 / T)*n_CO) + (A_f2*pow(T,b_f)*exp(-E_f1 / T)*n_O));
	}

public:
	double v_init, T_init, n_init, lambda;

	binary_flow(double v1, double T1, double n1) {

		v_init = v1;
		T_init = T1;
		n_init = n1;
		lambda = 1 / (sqrt(2.) * n_init * CO2_diameter);
	};

	void operator() (const state_type &x, state_type &dxdt, const double t)
	{
		// dimensionless variables:
		// v = x[0], T=x[1], n_CO2=x[2], n_CO=x[3], n_O=x[4]

		coeff_matrix.zeros();
		double n = x[2] + x[3]+x[4]; // dimensionless total numeric density
		double rho = m_CO_2*x[2]+m_CO*x[3]+m_O*x[4]; // dimensionless density

															  // some dimensional variables
		double p = n * k *x[1]; // pressure
		double T = x[1] ; // temperature
		double n_CO2 = x[2];
		double n_CO = x[3];
		double n_O = x[4];
		double v = x[0];

		// calculate energies
		
		double E = E_tr(T, n) + E_rot(T, n_CO2, n_CO) + E_f(n_CO2,n_CO,n_O) + E_vibr_CO(T, n_CO) + E_vibr_CO2(T, n_CO2);
		
		double R_diss = R(T,n_CO2,n_CO,n_O);

		coeff_matrix.at(0, 0) = rho*v;
		coeff_matrix.at(0, 1) = n*k;
		coeff_matrix.at(0, 2) = k*T;
		coeff_matrix.at(0, 3) = k*T;
		coeff_matrix.at(0, 4) = k*T;

		coeff_matrix.at(1, 0) = (E + p) ;
		coeff_matrix.at(1, 1) = v*dE_dT(T,n_CO2,n_CO,n_O);
		coeff_matrix.at(1, 2) = v*dE_dn_CO2(T);
		coeff_matrix.at(1, 3) = v*dE_dn_CO(T);
		coeff_matrix.at(1, 4) = v*dE_dn_O(T);

		coeff_matrix.at(2, 0) = n_CO2;
		coeff_matrix.at(2, 2) = v;

		coeff_matrix.at(3, 0) = n_CO;
		coeff_matrix.at(3, 3) = v;

		coeff_matrix.at(4, 0) = n_O;
		coeff_matrix.at(4, 4) = v;

		rhs[2] = R_diss;
		rhs[3] = -R_diss;
		rhs[4] = -R_diss;

		// get the derivatives vector
		arma::vec derivatives = arma::solve(coeff_matrix, rhs);

		dxdt[0] = derivatives[0];
		dxdt[1] = derivatives[1];
		dxdt[2] = derivatives[2];
		dxdt[3] = derivatives[3];
		dxdt[4]= derivatives[4];
	}

	
};

struct push_back_state_and_time
{
	std::vector< state_type >& m_states;
	std::vector< double >& m_times;

	push_back_state_and_time(std::vector< state_type > &states, std::vector< double > &times)
		: m_states(states), m_times(times) { }

	void operator()(const state_type &x, double t)
	{
		m_states.push_back(x);
		m_times.push_back(t);
	}
};

int main(void) {
	double T1 = 3047;
	double n1 = 1.7E23;
	double v1 = 309.11;

	binary_flow flow = binary_flow(v1, T1, n1);
	
	// std::cout << "Enter Mach Number: ";
	// std::cin >> M0;

	// state_type init_cond = flow.rg_init(p0, T0, M0);

	std::vector< double > init_cond = {v1, T1, n1, 0, 0};

	std::cout << init_cond[0] << " " << init_cond[1] << " " << init_cond[2] << std::endl;

	std::vector<state_type> macro_vars;
	std::vector<double> x;
	double max_x = 0.5;
	double x_step = 0.001;

	std::cout << "Enter maximum x: ";
	std::cin >> max_x;
	std::cout << "Enter step size:";
	std::cin >> x_step;
	std::cout << std::endl;

	size_t steps = integrate_adaptive(make_controlled(1.0e-10, 1.0e-6, error_stepper_type()), flow, init_cond, 0.0, max_x, x_step, push_back_state_and_time(macro_vars, x));

	for (size_t i = 0; i <= steps; i++)
	{
		std::cout << "x=" << x[i] << ", v=" << macro_vars[i][0];
		std::cout << ", T=" << macro_vars[i][1]  << ", n_CO2=" << macro_vars[i][2]  << ", n_CO=" << macro_vars[i][3] << ", n_O=" << macro_vars[i][4] << std::endl;
	}

	std::string a;
	std::cin >> a;
	return 1;
}