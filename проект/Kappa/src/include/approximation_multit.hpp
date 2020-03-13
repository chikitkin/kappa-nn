//
//  approximation_multit.hpp
//

#ifndef kappa_approximation_multit_hpp
#define kappa_approximation_multit_hpp

#include "approximation.hpp"
#include <algorithm>

namespace kappa {
    class ApproximationMultiT : public Approximation {
		using Approximation::k_diss;
		using Approximation::k_exch;
		using Approximation::k_VT;
		
    protected:
		double p_Z_vibr(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels);
		double p_avg_vibr_energy(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels);
		double p_avg_vibr_energy_sq(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels);
		double p_avg_vibr_i(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels);
		double p_avg_vibr_i_sq(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels);
		double p_avg_vibr_i_energy(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels);
		double p_c_vibr_T(double T, double T1, double mass,  const arma::vec & vibr_energy, int num_vibr_levels);
		double p_c_vibr_T1(double T, double T1, double mass, const arma::vec & vibr_energy, int num_vibr_levels);
		double p_c_W_T(double T, double T1, double mass, const arma::vec & vibr_energy, int num_vibr_levels);
		double p_c_W_T1(double T, double T1, double mass, const arma::vec & vibr_energy, int num_vibr_levels);
		int p_max_i(double T, double T1, double vibr_energy1, double vibr_frequency, double alpha);

    public:
        ApproximationMultiT(); 

		class k_Exchange {
		public:
			//k_Exchange();
			double k_exch_dir, k_exch_rev, k_exch_r_gain, k_exch_r_loos, k_exch_p_gain, k_exch_p_loos,
				k_exch_r_gain_ve, k_exch_r_loos_ve, k_exch_p_gain_ve, k_exch_p_loos_ve;
		};

		class k_Dissociation {
		public:
			//k_Dissociation();

			double k_diss, k_diss_i, k_diss_ve, k_rec, k_rec_i, k_rec_ve;
		};

		class k_VT_transition {
		public:
			//k_VT_transition();

			double k_VT_i, k_VT_ve;
		};

		
		double Z_rot(double T, const kappa::Molecule &molecule);
        double Z_vibr(double T, double T1, const kappa::Molecule &molecule);
		int max_i(double T, double T1, const kappa::Molecule & molecule);
		double avg_vibr_energy(double T, double T1, const kappa::Molecule &molecule);
		double avg_vibr_energy_sq(double T, double T1, const kappa::Molecule &molecule);
		double avg_vibr_i(double T, double T1, const kappa::Molecule &molecule); 
		double avg_vibr_i_sq(double T, double T1, const kappa::Molecule &molecule);
		double avg_vibr_i_energy(double T, double T1, const kappa::Molecule &molecule);
		double c_rot(double T, double T1, const kappa::Molecule &molecule);
		double c_vibr_T(double T, double T1, const kappa::Molecule &molecule);
		double c_vibr_T1(double T, double T1, const kappa::Molecule &molecule);
		double c_W_T(double T, double T1, const kappa::Molecule &molecule);
		double c_W_T1(double T, double T1, const kappa::Molecule &molecule);
		k_Dissociation k_diss(double T, double T1, kappa::Molecule const &molecule, kappa::Atom const &atom1, kappa::Atom const &atom2, kappa::Interaction const &interaction,
			kappa::models_k_diss model = kappa::models_k_diss::model_k_diss_vss_thresh_cmass_vibr);
		k_VT_transition k_VT(double T, double T1, kappa::Molecule const &molecule, kappa::Interaction const &interaction, int max_delta_i,
			kappa::models_k_vt model = kappa::models_k_vt::model_k_vt_vss_fho);
		k_Exchange k_exch(double T, double T1_r, double T1_p, kappa::Molecule const &molecule, kappa::Atom const &atom, 
			kappa::Molecule const &molecule_prod, kappa::Atom const &atom_prod, kappa::Interaction const &interaction, kappa::models_k_exch model);
		

		};

};

#endif /* kappa_approximation_sts_hpp */