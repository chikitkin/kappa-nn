//
//approximation_multit.cpp

#include "../include/approximation_multit.hpp"

kappa::ApproximationMultiT::ApproximationMultiT() :kappa::Approximation() { } 

	double kappa::ApproximationMultiT::Z_rot(double T, const kappa::Molecule &molecule) {
		return p_Z_rot(T, molecule.rot_energy[0][0], molecule.num_rot_levels[0][0], molecule.rot_symmetry);
	}

	double kappa::ApproximationMultiT::avg_vibr_energy(double T, double T1, const kappa::Molecule & molecule) {
		int i_star = max_i(T, T1, molecule);
		if (i_star >= molecule.num_vibr_levels[0]-1) {
			return p_avg_vibr_energy(T, T1, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
		}
		else {
			return p_avg_vibr_energy(T, T1, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);
			}
	}
	
	double kappa::ApproximationMultiT::avg_vibr_energy_sq(double T, double T1, const kappa::Molecule & molecule) {
		int i_star = max_i(T, T1, molecule);
		if (i_star >= molecule.num_vibr_levels[0] - 1) {
			return p_avg_vibr_energy_sq(T, T1, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
		}
		else {
			return p_avg_vibr_energy_sq(T, T1, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);
		}
	}

	double kappa::ApproximationMultiT::avg_vibr_i(double T, double T1, const kappa::Molecule & molecule) {
		int i_star = max_i(T, T1, molecule);
		if (i_star >= molecule.num_vibr_levels[0] - 1) {
			return p_avg_vibr_i(T, T1, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
		}
		else {
			return p_avg_vibr_i(T, T1, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);
		}
	}

	double kappa::ApproximationMultiT::avg_vibr_i_sq(double T, double T1, const kappa::Molecule & molecule) {
		int i_star = max_i(T, T1, molecule);
		if (i_star >= molecule.num_vibr_levels[0] - 1) {
			return p_avg_vibr_i_sq(T, T1, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
		}
		else {
			return p_avg_vibr_i_sq(T, T1, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);
		}
	}

	double kappa::ApproximationMultiT::avg_vibr_i_energy(double T, double T1, const kappa::Molecule & molecule) {
		int i_star = max_i(T, T1, molecule);
		if (i_star >= molecule.num_vibr_levels[0] - 1) {
			return p_avg_vibr_i_energy(T, T1, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
		}
		else {
			return p_avg_vibr_i_energy(T, T1, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);
		}
	}

	double kappa::ApproximationMultiT::c_rot(double T, double T1, const kappa::Molecule & molecule) {
		return p_c_rot(T, molecule.mass, molecule.rot_energy[0][0], molecule.num_rot_levels[0][0], molecule.rot_symmetry);
	}

	double kappa::ApproximationMultiT::c_vibr_T(double T, double T1, const kappa::Molecule & molecule) {
		int i_star = max_i(T, T1, molecule);
		if (i_star >= molecule.num_vibr_levels[0] - 1) {
			return p_c_vibr_T(T, T1, molecule.mass, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
		}
		else {
			return p_c_vibr_T(T, T1, molecule.mass, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);
		}
	}

	double kappa::ApproximationMultiT::c_vibr_T1(double T, double T1, const kappa::Molecule & molecule) {
		int i_star = max_i(T, T1, molecule);
		if (i_star >= molecule.num_vibr_levels[0] - 1) {
			return p_c_vibr_T1(T, T1, molecule.mass, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
		}
		else {
			return p_c_vibr_T1(T, T1, molecule.mass, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);
		}
	}

	double kappa::ApproximationMultiT::c_W_T(double T, double T1, const kappa::Molecule & molecule) {
		int i_star = max_i(T, T1, molecule);
		if (i_star >= molecule.num_vibr_levels[0] - 1) {
			return p_c_W_T(T, T1, molecule.mass, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
		}
		else {
			return p_c_W_T(T, T1, molecule.mass, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);
		}
	}

	double kappa::ApproximationMultiT::c_W_T1(double T, double T1, const kappa::Molecule & molecule) {
		int i_star = max_i(T, T1, molecule);
		if (i_star >= molecule.num_vibr_levels[0] - 1) {
			return p_c_W_T1(T, T1, molecule.mass, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
		}
		else {
			return p_c_W_T1(T, T1, molecule.mass, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);
		}
	}

	kappa::ApproximationMultiT::k_Dissociation kappa::ApproximationMultiT::k_diss(double T, double T1, kappa::Molecule const & molecule,
		kappa::Atom const &atom1, kappa::Atom const &atom2, kappa::Interaction const & interaction,	kappa::models_k_diss model){
		k_Dissociation V; 
		double tmp = 0.0, k_bf = 0.0;
		int max_level = max_i(T, T1, molecule);
		for (int i = 0; i < max_level; i++) {
			tmp = k_diss(T, molecule, interaction, i, model) * exp(((-(molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0])
				+ (molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0])*i) / (K_CONST_K * T))
				+ (-(molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * i / (K_CONST_K * T1)));
			k_bf = k_bf_diss(T, molecule, atom1, atom2, i, 0);
			V.k_diss += tmp;
			V.k_diss_i += i*tmp;
			V.k_diss_ve += (molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0])*tmp;
			V.k_rec += tmp*k_bf;
			V.k_rec_i += i*tmp*k_bf;
			V.k_rec_ve += (molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0])*tmp*k_bf;

		}
		V.k_diss = V.k_diss / Z_vibr(T, T1, molecule);
		V.k_diss_i = V.k_diss_i / Z_vibr(T, T1, molecule);
		V.k_diss_ve = V.k_diss_ve / Z_vibr(T, T1, molecule);
		V.k_rec = V.k_rec / Z_vibr(T, T1, molecule);
		V.k_rec_i = V.k_rec_i / Z_vibr(T, T1, molecule);
		V.k_rec_ve = V.k_rec_ve / Z_vibr(T, T1, molecule);
		
		return V;
	}

	kappa::ApproximationMultiT::k_VT_transition kappa::ApproximationMultiT::k_VT(double T, double T1, kappa::Molecule const & molecule, kappa::Interaction const & interaction, int max_delta_i, kappa::models_k_vt model) {
		k_VT_transition V;
		double tmp = 0.0;
		int max_level = max_i(T, T1, molecule);
		for (int i = 0; i <= max_level; i++)
			for (int delta_i = std::max(-i, -max_delta_i); delta_i <= std::min(max_level - i, max_delta_i); delta_i++) {
				if (delta_i != 0) {
					tmp = exp(((-(molecule.vibr_energy[0][delta_i + i] - molecule.vibr_energy[0][0])
						+ (molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0])*delta_i) / (K_CONST_K * T))
						+ (-(molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * delta_i / (K_CONST_K * T1)))
						*k_VT(T, molecule, interaction, i + delta_i, -delta_i, model)
						- exp(((-(molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0])
							+ (molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0])*i) / (K_CONST_K * T))
							+ (-(molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * i / (K_CONST_K * T1)))
						*k_VT(T, molecule, interaction, i, delta_i, model);
					V.k_VT_i += i*tmp;
					V.k_VT_ve += (molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0])*tmp;

				}
			}
		return V;
	}

	kappa::ApproximationMultiT::k_Exchange kappa::ApproximationMultiT::k_exch(double T, double T1_r, double T1_p, kappa::Molecule const & molecule, kappa::Atom const & atom, kappa::Molecule const & molecule_prod, kappa::Atom const & atom_prod, kappa::Interaction const & interaction, kappa::models_k_exch model) {
		k_Exchange V;
		double k_exch_dir = 0, k_exch_rev = 0, k_exch_r_gain = 0, k_exch_r_loos = 0, k_exch_p_gain = 0, k_exch_p_loos = 0,
			k_exch_r_gain_ve = 0, k_exch_r_loos_ve = 0, k_exch_p_gain_ve = 0, k_exch_p_loos_ve = 0, tmp_dir = 0, tmp_rev = 0;
		int max_level_r = max_i(T, T1_r, molecule);
		int max_level_p = max_i(T, T1_p, molecule_prod);
		for (int i = 0; i < max_level_r; i++)
			for (int k = 0; k < max_level_p; k++) {
				tmp_dir = k_exch(T, molecule, atom, molecule_prod, interaction, i, k, model);
				tmp_rev = k_exch(T, molecule, atom, molecule_prod, interaction, i, k, model)
					*k_bf_exch(T, molecule, atom, molecule_prod, atom_prod, interaction, i, k, 0, 0);
				k_exch_dir += exp(((-(molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0]) + (molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0])*i) / (K_CONST_K * T))
					+ (-(molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * i / (K_CONST_K * T1_r)))
					*tmp_dir;
				k_exch_rev += exp(((-(molecule_prod.vibr_energy[0][k] - molecule_prod.vibr_energy[0][0]) + (molecule_prod.vibr_energy[0][1] - molecule_prod.vibr_energy[0][0])*k) / (K_CONST_K * T))
					+ (-(molecule_prod.vibr_energy[0][1] - molecule_prod.vibr_energy[0][0]) * k / (K_CONST_K * T1_p)))
					*tmp_rev;
				k_exch_r_gain += i*exp(((-(molecule_prod.vibr_energy[0][k] - molecule_prod.vibr_energy[0][0]) + (molecule_prod.vibr_energy[0][1] - molecule_prod.vibr_energy[0][0])*k) / (K_CONST_K * T))
					+ (-(molecule_prod.vibr_energy[0][1] - molecule_prod.vibr_energy[0][0]) * k / (K_CONST_K * T1_p)))
					*tmp_rev;
				k_exch_r_loos += i*exp(((-(molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0]) + (molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0])*i) / (K_CONST_K * T))
					+ (-(molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * i / (K_CONST_K * T1_r)))
					*tmp_dir;
				k_exch_p_gain += k*exp(((-(molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0]) + (molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0])*i) / (K_CONST_K * T))
					+ (-(molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * i / (K_CONST_K * T1_r)))
					*tmp_dir;
				k_exch_p_loos += k*exp(((-(molecule_prod.vibr_energy[0][k] - molecule_prod.vibr_energy[0][0]) + (molecule_prod.vibr_energy[0][1] - molecule_prod.vibr_energy[0][0])*k) / (K_CONST_K * T))
					+ (-(molecule_prod.vibr_energy[0][1] - molecule_prod.vibr_energy[0][0]) * k / (K_CONST_K * T1_p)))
					*tmp_rev;
				k_exch_r_gain_ve += (molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0])
					*exp(((-(molecule_prod.vibr_energy[0][k] - molecule_prod.vibr_energy[0][0]) + (molecule_prod.vibr_energy[0][1] - molecule_prod.vibr_energy[0][0])*k) / (K_CONST_K * T))
						+ (-(molecule_prod.vibr_energy[0][1] - molecule_prod.vibr_energy[0][0]) * k / (K_CONST_K * T1_p)))
					*tmp_rev;
				k_exch_r_loos_ve += (molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0])
					*exp(((-(molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0]) + (molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0])*i) / (K_CONST_K * T))
						+ (-(molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * i / (K_CONST_K * T1_r)))
					*tmp_dir;
				k_exch_p_gain_ve += (molecule_prod.vibr_energy[0][k] - molecule_prod.vibr_energy[0][0])
					*exp(((-(molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0]) + (molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0])*i) / (K_CONST_K * T))
						+ (-(molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * i / (K_CONST_K * T1_r)))
					*tmp_dir;
				k_exch_p_loos_ve += (molecule_prod.vibr_energy[0][k] - molecule_prod.vibr_energy[0][0])
					*exp(((-(molecule_prod.vibr_energy[0][k] - molecule_prod.vibr_energy[0][0]) + (molecule_prod.vibr_energy[0][1] - molecule_prod.vibr_energy[0][0])*k) / (K_CONST_K * T))
						+ (-(molecule_prod.vibr_energy[0][1] - molecule_prod.vibr_energy[0][0]) * k / (K_CONST_K * T1_p)))
					*tmp_rev;
			}

		V.k_exch_dir = k_exch_dir / Z_vibr(T, T1_r, molecule);
		V.k_exch_rev = k_exch_rev / Z_vibr(T, T1_p, molecule_prod);
		V.k_exch_r_gain = k_exch_r_gain / Z_vibr(T, T1_p, molecule_prod);
		V.k_exch_r_loos = k_exch_r_loos / Z_vibr(T, T1_r, molecule);
		V.k_exch_p_gain = k_exch_p_gain / Z_vibr(T, T1_r, molecule);
		V.k_exch_p_loos = k_exch_p_loos / Z_vibr(T, T1_p, molecule_prod);
		V.k_exch_r_gain_ve = k_exch_r_gain_ve / Z_vibr(T, T1_p, molecule_prod);
		V.k_exch_p_loos_ve = k_exch_r_loos_ve / Z_vibr(T, T1_r, molecule);
		V.k_exch_r_gain_ve = k_exch_p_gain_ve / Z_vibr(T, T1_r, molecule);
		V.k_exch_p_loos_ve = k_exch_p_loos_ve / Z_vibr(T, T1_p, molecule_prod);

		return V;
	}

	double kappa::ApproximationMultiT::Z_vibr(double T, double T1, const kappa::Molecule &molecule) {
		if ((molecule.anharmonic_spectrum == false) || (T >= T1)) {
			return p_Z_vibr(T, T1, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
		}
		else {
			int i_star = max_i(T, T1, molecule);
			if (i_star < 1) {
				std::string error_string = "No Treanor distribution possible for such values of T, T1";
				throw kappa::ModelParameterException(error_string.c_str());
			}
			else if (i_star >= molecule.num_vibr_levels[0] - 1) {
				return p_Z_vibr(T, T1, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
			}
			else {
				return p_Z_vibr(T, T1, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);
			}
		}
	}

	int kappa::ApproximationMultiT::max_i(double T, double T1, const kappa::Molecule &molecule) {
		if ((molecule.anharmonic_spectrum == false) || (T >= T1)) {
			return molecule.num_vibr_levels[0];
		}
		else {
			int i_star = p_max_i(T, T1, molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0], K_CONST_C * molecule.vibr_frequency[0], molecule.vibr_we_xe[0] / molecule.vibr_frequency[0]);
			if (i_star < 1) {
				std::string error_string = "No Treanor distribution possible for such values of T, T1";
				throw kappa::ModelParameterException(error_string.c_str());
			}
			else if (i_star >= molecule.num_vibr_levels[0] - 1) {
				return molecule.num_vibr_levels[0] - 1;
			}
			else {
				return p_max_i(T, T1, molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0], K_CONST_C * molecule.vibr_frequency[0], molecule.vibr_we_xe[0] / molecule.vibr_frequency[0]);
			}
		}
	}

	//privat methods:

	double kappa::ApproximationMultiT::p_Z_vibr(double T, double T1, const arma::vec &vibr_energy, int num_vibr_levels) {
		return arma::sum(arma::exp(((-vibr_energy+ (vibr_energy[1])*arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels) )/ (K_CONST_K * T))
			+ (-vibr_energy[1]* arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels) / (K_CONST_K * T1)) ));
	}

	double kappa::ApproximationMultiT::p_avg_vibr_energy(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels) {
		return (arma::dot(vibr_energy, arma::exp(((-vibr_energy
			+ (vibr_energy[1])*arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T))
			+ (-vibr_energy[1] * arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T1))))
			/ p_Z_vibr(T, T1, vibr_energy, num_vibr_levels);
	}
	
	double kappa::ApproximationMultiT::p_avg_vibr_energy_sq(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels) {
		return (arma::dot(vibr_energy%vibr_energy, arma::exp(((-vibr_energy
			+ vibr_energy[1]*arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T))
			+ (-vibr_energy[1] * arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T1))))
			/ p_Z_vibr(T, T1, vibr_energy, num_vibr_levels);
	}

	double kappa::ApproximationMultiT::p_avg_vibr_i(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels) {
		return (arma::dot(arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels), arma::exp(((-vibr_energy
			+ (vibr_energy[1])*arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T))
			+ (-vibr_energy[1] * arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T1))))
			/ p_Z_vibr(T, T1, vibr_energy, num_vibr_levels);
	}

	double kappa::ApproximationMultiT::p_avg_vibr_i_sq(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels) {
		return (arma::dot(arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)%arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels),
			arma::exp(((-vibr_energy+ (vibr_energy[1])*arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T))
				+ (-vibr_energy[1] * arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T1))))
			/ p_Z_vibr(T, T1, vibr_energy, num_vibr_levels);
	}

	double kappa::ApproximationMultiT::p_avg_vibr_i_energy(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels) {
		return (arma::dot(arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)% vibr_energy, arma::exp(((-vibr_energy
			+ (vibr_energy[1])*arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T))
			+ (-vibr_energy[1] * arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T1))))
			/ p_Z_vibr(T, T1, vibr_energy, num_vibr_levels);
	}

	double kappa::ApproximationMultiT::p_c_vibr_T(double T, double T1, double mass, const arma::vec & vibr_energy, int num_vibr_levels) {
		double tmp = p_avg_vibr_energy(T, T1, vibr_energy, num_vibr_levels);
		return (p_avg_vibr_energy_sq(T, T1, vibr_energy, num_vibr_levels) - tmp*tmp			
			+ vibr_energy[1]* p_avg_vibr_i(T, T1, vibr_energy, num_vibr_levels)*p_avg_vibr_energy(T, T1, vibr_energy, num_vibr_levels) 
			- vibr_energy[1]* p_avg_vibr_i_energy(T, T1, vibr_energy, num_vibr_levels) )/(K_CONST_K*K_CONST_K*T*T);
	}

	double kappa::ApproximationMultiT::p_c_vibr_T1(double T, double T1, double mass, const arma::vec & vibr_energy, int num_vibr_levels) {
		return ( vibr_energy[1]*(p_avg_vibr_i_energy(T, T1, vibr_energy, num_vibr_levels)
			- p_avg_vibr_i(T, T1, vibr_energy, num_vibr_levels)*p_avg_vibr_energy(T, T1, vibr_energy, num_vibr_levels)) )/(K_CONST_K*K_CONST_K*T1*T1);
	}

	double kappa::ApproximationMultiT::p_c_W_T(double T, double T1, double mass, const arma::vec & vibr_energy, int num_vibr_levels) {
		double tmp = p_avg_vibr_i(T, T1, vibr_energy, num_vibr_levels);
		return (vibr_energy[1] * (p_avg_vibr_i_energy(T, T1, vibr_energy, num_vibr_levels) 
			- vibr_energy[1]* p_avg_vibr_i_sq(T, T1, vibr_energy, num_vibr_levels)
			- p_avg_vibr_i(T, T1, vibr_energy, num_vibr_levels)*p_avg_vibr_energy(T, T1, vibr_energy, num_vibr_levels) 
			+ vibr_energy[1]*tmp*tmp)) / (K_CONST_K*K_CONST_K*T*T);
	}

	double kappa::ApproximationMultiT::p_c_W_T1(double T, double T1, double mass, const arma::vec & vibr_energy, int num_vibr_levels) {
		double tmp = p_avg_vibr_i(T, T1, vibr_energy, num_vibr_levels);
		return ( ( (vibr_energy[1]) * (vibr_energy[1]))*(p_avg_vibr_i_sq(T, T1, vibr_energy, num_vibr_levels)- tmp*tmp ) )/ (K_CONST_K*K_CONST_K*T1*T1);
	} 

	int kappa::ApproximationMultiT::p_max_i(double T, double T1, double vibr_energy1, double vibr_frequency, double alpha) {
		int p_i = -1;
		p_i = (vibr_energy1*T) / (2 * alpha*vibr_frequency*T1*K_CONST_H) + 0.5;
		return p_i;
	}
