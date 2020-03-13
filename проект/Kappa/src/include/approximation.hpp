//
//  approximation.hpp
//

/*
<topdoc>
    <label>approximation-label</label>
    <type>classdoc</type>
    <doc lang="eng">
        The Approximation class contains
    </doc>
    <doc lang="rus">
        Класс Approximation содержит методы
    </doc>
</topdoc>
*/


#ifndef kappa_approximation_hpp
#define kappa_approximation_hpp

#include <armadillo>
#include <vector>
#include "interaction.hpp"
#include "numeric.hpp"
#include "models.h"


namespace kappa {

    class Approximation {
        // base approximation class
    protected:

        const static double p4e_n2n_vt_bij1[5][5];
        const static double p4e_n2n_vt_bij2[5][5];
        const static double p4e_n2n_vt_cij2[5][5];
        const static double p4e_n2n_vt_bij3[5][5];
        const static double p4e_n2n_vt_cij3[5][5];
        const static double p4e_n2n_vt_bij4[5][5];
        const static double p4e_n2n_vt_cij4[5][5];
        const static double p4e_n2n_vt_bij5[5][5];
        const static double p4e_n2n_vt_cij5[5][5];
        const static double p4e_n2n_vt_bij6[5][5];
        const static double p4e_n2n_vt_cij6[5][5];


        const static double p4e_o2o_vt_Cijk1[3][5][4];
        const static double p4e_o2o_vt_Cijk2[3][5][4];
        const static double p4e_o2o_vt_Cijk3[3][5][4];
        const static double p4e_o2o_vt_Cijk4[3][5][4];

        const static double p4e_n2n_diss_bjk[4][5];
        
        const static double p4e_o2o_diss_bjk[3][8];

        const static double cc_ar[5];
        const static double cc_cl[5];

        static int p_max_electron_level(const arma::vec &electron_energy, double ionization_potential, double Delta_E);

        static double convert_vibr_ladder_N2(double vibr_energy);
        static double convert_vibr_ladder_O2(double vibr_energy);

        static double p_Z_rot(double T, const arma::vec &rot_energy, int num_rot_levels, int rot_symmetry);
        // calculate the rotational partition function given a temperature T, values of the energies of the rotational levels,
        // the number of rotational levels and whether the molecule is symmetric

        static double p_Z_vibr_eq(double T, const arma::vec &vibr_energy);
        // calculate the equilibrium vibrational partition function given a temperature T and values of the energies of the vibrational levels

        static double p_Z_electron(double T, const arma::vec &electron_energy, const arma::Col<unsigned long> &statistical_weight, int num_electron_levels);
        // calculate the electron partition function given a temperature T, values of the energies of the electron levels,
        // the statistical weights and the number of electron levels

        static double Z_diss(double T, double U, const arma::vec &electron_energy, const arma::Col<unsigned long> &statistical_weight, int n_electron_levels,
                             const std::vector<arma::vec> &vibr_energy, int i, int e);

        static double Z_diss(double T, const arma::vec &electron_energy, const arma::Col<unsigned long> &statistical_weight, int n_electron_levels,
                             const std::vector<arma::vec> &vibr_energy, const std::vector<int> &num_vibr_levels, int i, int e);

        static double Z_diss(double T, double U, const arma::vec &vibr_energy, int i);

        static double Z_diss(double T, const arma::vec &vibr_energy, int num_vibr_levels, int i);

        static double C_aliat(double T, const arma::vec &electron_energy, const arma::Col<unsigned long> &statistical_weight, int n_electron_levels,
                              const std::vector<arma::vec> &vibr_energy, const std::vector<int> &num_vibr_levels,
                              double vibr_energy_product, double activation_energy, int i, int e, double U);

        static double C_aliat(double T, const arma::vec &electron_energy, const arma::Col<unsigned long> &statistical_weight, int n_electron_levels,
                              const std::vector<arma::vec> &vibr_energy, const std::vector<int> &num_vibr_levels,
                              double vibr_energy_product, double activation_energy, int i, int e);

        static double p_avg_rot_energy(double T, const arma::vec &rot_energy, int num_rot_levels, int rot_symmetry);
        // calculate the rotational energy averaged over the rotational spectrum given a temperature T, values of the energies of the rotational levels,
        // the number of rotational levels and whether the molecule is symmetric

        static double p_avg_rot_energy_sq(double T, const arma::vec &rot_energy, int num_rot_levels, int rot_symmetry);
        // calculate the squared rotational energy averaged over the rotational spectrum given a temperature T, values of the energies of the rotational levels,
        // the number of rotational levels and whether the molecule is symmetric

        static double p_c_tr(double T, double mass);
        // calculate the specific translational heat given a temperature T

        static double p_c_rot(double T, double mass, const arma::vec &rot_energy, int num_rot_levels, int rot_symmetry);
        // calculate the specific heat capacity of rotational degrees of freedom given a temperature T, the mass of the molecule, values of the energies of the rotational levels,
        // the number of rotational levels and whether the molecule is symmetric

        static double vel_avg_vv(double rel_vel, double coll_mass, double Delta_E_vibr);
        // calculate the mean of the relative velocity before and after a VV exchange using the energy conservation law, Delta_E_vibr is the difference
        // of the total vibrational energies of the molecules before and after the collision

        static double vel_avg_vt(double rel_vel, double coll_mass, double Delta_E_vibr);
        // calculate the mean of the relative velocity before and after a VT transition using the energy conservation law, Delta_E_vibr is the difference
        // of the vibrational energy of the molecule before and after the collision
        static double min_vel_VT(double coll_mass, double vibr_energy_before, double vibr_energy_after);
        static double min_vel_diss(double coll_mass, double diss_energy, double vibr_energy);
        // calculate the minimum velocity for which a dissociation probability is non-zero, accounting for vibrational energy
        static double min_vel_diss_ILT_N2N(double coll_mass, double vibr_energy, double i);
        static double min_vel_diss_ILT_O2O(double coll_mass, double i);

        static double probability_VV_FHO(double rel_vel, double coll_mass,
                                         double Delta_E_vibr,
                                         int i, int k, int delta_i,
                                         double omega1, double omega2, double omega, double alpha_FHO);
        static double probability_VT_FHO(double rel_vel, double coll_mass, double osc_mass,
                                         double Delta_E_vibr,
                                         int i, int delta_i,
                                         double omega, double alpha_FHO);
        
       
        static double p_probability_diss(double rel_vel, double coll_mass, double diss_energy, double vibr_energy, bool center_of_mass=true);
        // calculate the probability of dissociation, accounting for vibrational energy

        static double crosssection_elastic_RS(double diameter);
        // calculate the elastic crossection sigma_{tot} using the Rigid Sphere interaction potential

        static double crosssection_elastic_VSS(double rel_vel, double coll_mass, double vss_c, double vss_omega);

        static double crosssection_elastic_VSS(double rel_vel, double vss_c_cs, double vss_omega);

        static double crosssection_VV_FHO_RS();
        static double crosssection_VV_FHO_VSS();
        static double crosssection_VT_FHO_RS(double rel_vel, double coll_mass, double diameter, double osc_mass,
                                             double Delta_E_vibr, int i, int delta_i,
                                             double omega, double alpha_FHO);
        static double crosssection_VT_FHO_VSS(double rel_vel, double coll_mass, double vss_c_cs, double vss_omega, double osc_mass,
                                              double Delta_E_vibr, int i, int delta_i,
                                              double omega, double alpha_FHO);

        static double crosssection_diss_RS(double rel_vel, double coll_mass, double diameter, double diss_energy, double vibr_energy, bool center_of_mass=true);
        // dissociation cross-section (RS model), accounting for vibrational energy


        static double crosssection_diss_VSS(double rel_vel, double coll_mass, double vss_c_cs, double vss_omega, double diss_energy, double vibr_energy, bool center_of_mass=true);
        // dissociation cross-section (VSS model), accounting for vibrational energy

        static double crosssection_diss_ILT_N2N(double rel_vel, double coll_mass, double vibr_energy, double i);
        // dissociation cross-section for the N2+N collision obtained using the Inverse Laplace Transform applied to QCT data

        static double crosssection_diss_ILT_O2O(double rel_vel, double coll_mass, double vibr_energy, double i);
        // dissociation cross-section for the O2+O collision obtained using the Inverse Laplace Transform applied to QCT data

        static double p_Z_coll(double T, double n, double coll_mass, double diameter);
        static double integral_rigid_RS(double T, int degree, double coll_mass, double diameter);
        static double integral_rigid_VSS(double T, int degree, double coll_mass, double diameter);
        static double integral_VV_FHO_RS();
        static double integral_VV_FHO_VSS();
        static double integral_VT_FHO_RS(double T, int degree, double coll_mass, double diameter, double osc_mass,
                                         double Delta_E_vibr, int i, int delta_i, double omega, double alpha_FHO);
        static double integral_VT_FHO_VSS(double T, int degree, double coll_mass, double vss_c_cs, double vss_omega, double osc_mass,
                                          double Delta_E_vibr, int i, int delta_i,
                                          double omega, double alpha_FHO);
        static double integral_diss_RS(double T, int degree, double coll_mass, double diameter, double diss_energy, double vibr_energy, bool center_of_mass=true);
        static double integral_diss_VSS(double T, int degree, double coll_mass, double vss_c_cs, double vss_omega, double diss_energy, double vibr_energy, bool center_of_mass=true);
        static double integral_diss_ILT_N2N(double T, int degree, double coll_mass, double vibr_energy, double i);
        static double integral_diss_ILT_O2O(double T, int degree, double coll_mass, double vibr_energy, double i);

        static double P_SSH_VT_10(double T, double coll_mass, double omega_e, double epsilon, double diameter, double r_e);
        static double P_SSH_VV_01(double T, double omega_e, double epsilon, double osc_mass, double diameter, double r_e);
        static double k_VV_FHO_RS();
        static double k_VV_FHO_VSS();
        static double k_VV_SSH();
        static double k_VT_FHO_RS(double T, double coll_mass, double diameter, double osc_mass,
                                  double Delta_E_vibr, int i, int delta_i, double omega, double alpha_FHO);
        static double k_VT_FHO_VSS(double T, double coll_mass, double vss_c_cs, double vss_omega, double osc_mass,
                                   double Delta_E_vibr, int i, int delta_i,
                                   double omega, double alpha_FHO);
        static double k_VT_SSH(double T, int i, double coll_mass, double diameter, double omega_e, double epsilon, double r_e);
        static double k_VT_SSH(double T, int i, double coll_mass, double diameter, double omega_e, double epsilon, double r_e, double Delta_E_vibr, double vibr_energy_1); // anharmonic
        static double k_VV_Billing_N2N2(double T, int i, int k);
        static double k_VT_Billing_N2N2(double T, int i);
        static double k_VT_Billing_N2N(double T, int i, int delta_i);
        static double k_VT_N2N_p4e(double T, int i, int delta_i);
        static double k_VT_O2O_p4e(double T, int i, int delta_i);

        static double k_exch_WRFP(double T, double vibr_energy, double activation_energy, double alpha_exch, double beta_exch, double A_exch, double n_exch);
        
        static double k_diss_RS(double T, double coll_mass, double diameter, double diss_energy, double electron_vibr_energy, bool center_of_mass=true);
        static double k_diss_VSS(double T, double coll_mass, double vss_c_cs, double vss_omega, double diss_energy, double electron_vibr_energy, bool center_of_mass=true);
        static double k_diss_ILT_N2N(double T, double coll_mass, double vibr_energy, double i);
        static double k_diss_ILT_O2O(double T, double coll_mass, double vibr_energy, double i);
        static double diss_rate_N2N_p4e(double T, double i); // phys4entry
        static double diss_rate_O2O_p4e(double T, double i); // phys4entry

        static double k_Arrhenius(double T, double arrhenius_A, double arrhenius_n, double energy);
        // Calculates the equilibrium Arrhenius coefficient defined as
        // arrhenius_A * T^arrhenius_n * exp(-energy / kT)
        //
        // takes as parameters the temperature T, the Arrhenius model parameters arrhenius_A and arrhenius_n, and the dissociation/activation energy

        static double k_diss_Treanor_Marrone();

        static double p_rot_collision_number_parker(double T, double xi_inf, double epsilon);
        static double p_vibr_relaxation_time_MW(double T, double p, double char_vibr_temp, double coll_mass);
        static double p_vibr_relaxation_time_Park_corr(double T, double n, double coll_mass, double crosssection);

        static double omega_integral_RS(double T, int l, int r, double diameter, double coll_mass);
        // calculates the \Omega^{(l,r)} integral using the Rigid Spheres model
        //
        // takes as parameters the temperature T, l and r powers of the integrand, the collision diameter and collision-reduced mass
        
        static double omega_integral_VSS(double T, int l, int r, double coll_mass, double vss_c_cs, double vss_omega, double vss_alpha);
        // calculates the \Omega^{(l,r)} integral using the Rigid Spheres model
        //
        // takes as parameters the temperature T, l and r powers of the integrand, the vss_c and vss_omega VSS potential parameters
        //
        // if (l,r) != (1,1) or (l,r) != (2,2), returns -1

        static double omega_integral_LennardJones(double T, int l, int r, double epsilon);
        // calculates the reduced \Omega^{(l,r)}* integral using the Lennard-Jones potential
        //
        // takes as parameters the temperature T, l and r powers of the integrand, and the epsilon interaction parameter
        //
        // Approximation formulas according to [Kustova, Nagnibeda]
        //
        // if (l,r) != (1,1) or (l,r) != (2,2), returns -1

		static double omega_integral_Born_Mayer(double T, int l, int r, double beta, double phi_zero, double diameter, double epsilon);
		// calculates the reduced \Omega^{(l,r)}* integral using the Born-Mayer potential
		//
		// takes as parameters the temperature T, l and r powers of the integrand, the beta, phi_zero interaction parameters,
		// the collision diameter and the epsilon interaction parameter
		//
		// Approximation formulas according to [Kustova, Nagnibeda]
		//
		// if (l,r) != (1,1) or (l,r) != (2,2), returns -1

        static double omega_integral_ESA_nn(double T, int l, int r, double beta, double epsilon_zero, double r_e);
		// calculates the reduced \Omega^{(l,r)}* integral using the potential given in ESA STR 256 for neutral-neutral interactions
        //
        // takes as parameters the temperature T, l and r powers of the integrand, the beta, epsilon_zero, and r_e interaction parameters
        //
        // if (l,r) not in the set ((1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (2, 2), (2, 3), (2, 4), (3, 3), (4, 4)), returns -1

        static double omega_integral_ESA_cn(double T, int l, int r, double beta, double epsilon_zero, double r_e);
        // neutral-ion
        // -//-
        //
        //
        //

        static double omega_integral_ESA_cc(double T, int l, int r, int charge1, int charge2, double debye_length);
        // charged-charged

        static double omega_integral_ESA_ne(double T, double g1, double g2, double g3, double g4, double g5, double g6, double g7, double g8, double g9, double g10, double diameter);
        // neutral - electron

        static double omega_integral_ESA_hydrogen(double T, double diameter, double a1, double a2, double a3, double a4, double a5, double a6, double a7);

        static double omega_integral_ESA_H2H2(double T, int l, int r, double diameter);
        static double omega_integral_ESA_H2H(double T, int l, int r, double diameter);
        static double omega_integral_ESA_HH(double T, int l, int r, double diameter);
        // for H+H, H2+H2, H2+H


        static double omega_integral_ESA_cn_corr(double T, double d1, double d2, double d3, double beta, double r_e);

        static double p_k_bf_VT(double T, double vibr_energy_before, double vibr_energy_after, const arma::vec &rot_energy_before, int num_rot_levels_before, 
                              const arma::vec &rot_energy_after, int num_rot_levels_after, int rot_symmetry, bool is_rigid_rotator);

        static double p_k_bf_VV(double T, double vibr_energy1_before, double vibr_energy1_after, double vibr_energy2_before, double vibr_energy2_after,
                              const arma::vec &rot_energy1_before, int num_rot_levels1_before, const arma::vec &rot_energy1_after, int num_rot_levels1_after, int rot_symmetry1, 
                              const arma::vec &rot_energy2_before, int num_rot_levels2_before, const arma::vec &rot_energy2_after, int num_rot_levels2_after, int rot_symmetry2);

        static double p_k_bf_exch(double T, double molecule_before_mass, double atom_before_mass, double molecule_after_mass, double atom_after_mass,
                                double diss_energy_before, double diss_energy_after, double vibr_energy_before, double vibr_energy_after,
                                const arma::vec &rot_energy_before, int num_rot_levels_before, int rot_symmetry_before,
                                const arma::vec &rot_energy_after, int num_rot_levels_after, int rot_symmetry_after);

        static double p_k_bf_diss(double T, double molecule_mass, double atom1_mass, double atom2_mass, double diss_energy, double vibr_energy,
                                const arma::vec &rot_energy, int num_rot_levels, int rot_symmetry);

    public:
        // <constructor_docs>
        Approximation();
        /*
        <cdoc>
            <doc lang="eng">
                Create an Approximation instance
            </doc>
            <doc lang="rus">
                Создает объект типа Approximation
            </doc>
        </cdoc>
        */
        // </constructor_docs>


        // <method_docs>
        double debye_length(double T, const arma::vec &concentrations, const arma::vec &charges);
        /*
        <fdoc>
            <doc lang="rus">
                Производит расчет Дебаевской длины

                **Параметры**:
                    
                    * ``T`` - температура смеси

                    * ``const arma::vec &concentrations`` - числовые плотности частиц

                    * ``const arma::vec &charges`` - (размерные) заряды частиц
            </doc>
        </fdoc>
        */

        int max_electron_level(const kappa::Atom &atom, double Delta_E);
        /*
        <fdoc>
            <i1><math>\varepsilon^{el}_{e} \le \varepsilon_{ion} - \Delta E</math>, <math>\varepsilon^{el}_{e+1} > \varepsilon_{ion} - \Delta E</math></i1>
            <doc lang="eng">
                Calculate the maximum electron level <imath>e</imath> to take into account for an atom, which satisfies the relations:
                <insert>i1</insert>

                **Parameters**:
                
                    * ``const kappa::Atom &atom`` - the atom being considered

                    * ``double Delta_E`` - the value of <imath>\Delta E</imath> (in Joules)
            </doc>
            <doc lang="rus">
                Рассчитывает номер максимального допустимого электронного уровня <imath>e</imath> для атома, который удовлетворяет соотношениям:
                <insert>i1</insert>

                **Параметры**:
                    
                    * ``kappa::Atom const &atom`` - рассматриваемый атом

                    * ``double Delta_E`` - значение величины <imath>\Delta E</imath> (в Джоулях)
            </doc>
        </fdoc>
        */

        double Z_rot(double T, const kappa::Molecule &molecule, int i=0, int e=0);
        /*
        <fdoc>
            <doc lang="eng">
                Calculate the rotational partition function

                **Parameters**:
                
                    * ``double T`` - the temperature

                    * ``const kappa::Molecule &molecule`` - the molecule for which the calculations are perfomed

                    * ``int i`` - the vibrational level of the molecule (default value is 0)

                    * ``int e`` - the electron level of the molecule (default value is 0)
            </doc>
            <doc lang="rus">
                Рассчитывает вращательную статистическую сумму

                **Параметры**:
                
                    * ``double T`` - температура

                    * ``const kappa::Molecule &molecule`` - молекула, для которой производится расчет

                    * ``int i`` - Колебательный уровень молекулы (значение по умолчанию 0)

                    * ``int e`` - Электронный уровень молекулы (значение по умолчанию 0)
            </doc>
        </fdoc>
        */

        double Z_vibr_eq(double T, const kappa::Molecule &molecule, int e=0);
        /*
        <fdoc>
            <doc lang="eng">
                Calculate the equilibrium vibrational partition function

                **Parameters**:
                    
                    * ``double T`` - the temperature

                    * ``const kappa::Molecule &molecule`` - the molecule for which the calculations are perfomed

                    * ``int e`` - the electron level of the molecule (default value is 0)
            </doc>
            <doc lang="rus">
                Рассчитывает равновесную колебательную статистическую сумму

                **Параметры**:
                
                
                    * ``double T`` - температура

                    * ``const kappa::Molecule &molecule`` - молекула, для которой производится расчет

                    * ``int e`` - Электронный уровень молекулы (значение по умолчанию 0)
            </doc>
        </fdoc>
        */

        double Z_electron(double T, const kappa::Molecule &molecule, int num_electron_levels);
        /*
        <fdoc>
            <doc lang="eng">
                Calculate the equilibrium electron partition function

                **Parameters**:
                    
                    * ``double T`` - the temperature

                    * ``const kappa::Molecule &molecule`` - the molecule for which the calculations are perfomed

                    * ``int num_electron_levels`` - the amount of electron states to account for; if set to -1, all possible states will be included
            </doc>
            <doc lang="rus">
                Рассчитывает равновесную электронную статистическую сумму

                **Параметры**:
                   
                    * ``double T`` - температура

                    * ``const kappa::Molecule &molecule`` - молекула, для которой производится расчет

                    * ``int num_electron_levels`` - число учитываемых электронных уровней; если равно -1, то учитываются все электронные состояния
            </doc>
        </fdoc>
        */

        double avg_rot_energy(double T, const kappa::Molecule &molecule, int i=0, int e=0);
        /*
        <fdoc>
            <doc lang="eng">
                Calculate the averaged rotational energy

                **Parameters**:
                            
                    * ``double T`` - the temperature

                    * ``const kappa::Molecule &molecule`` - the molecule for which the calculations are perfomed

                    * ``int i`` - the vibrational level of the molecule (default value is 0)

                    * ``int e`` - the electron level of the molecule (default value is 0)
            </doc>
            <doc lang="rus">
                Рассчитывает осредненную вращательную энергию

                **Параметры**:
                
                    * ``double T`` - температура

                    * ``const kappa::Molecule &molecule`` - молекула, для которой производится расчет

                    * ``int i`` - Колебательный уровень молекулы (значение по умолчанию 0)

                    * ``int e`` - Электронный уровень молекулы (значение по умолчанию 0)
            </doc>
        </fdoc>
        */

        double avg_rot_energy_sq(double T, const kappa::Molecule &molecule, int i=0, int e=0);
        /*
        <fdoc>
            <doc lang="eng">
                Calculate the averaged squared rotational energy

                **Parameters**:
                
                    * ``double T`` - the temperature

                    * ``const kappa::Molecule &molecule`` - the molecule for which the calculations are perfomed

                    * ``int i`` - the vibrational level of the molecule (default value is 0)

                    * ``int e`` - the electron level of the molecule (default value is 0)
            </doc>
            <doc lang="rus">
                Рассчитывает осредненный квадрат вращательной энергии

                **Параметры**:
                                
                    * ``double T`` - температура

                    * ``const kappa::Molecule &molecule`` - молекула, для которой производится расчет

                    * ``int i`` - Колебательный уровень молекулы (значение по умолчанию 0)

                    * ``int e`` - Электронный уровень молекулы (значение по умолчанию 0)
            </doc>
        </fdoc>
        */

        double c_tr(double T, const kappa::Particle &particle);
        /*
        <fdoc>
            <doc lang="eng">
                Calculate the specific heat of translational degrees of freedom

                **Parameters**:
                
                    * ``double T`` - the temperature

                    * ``const kappa::Particle &particle`` - the particle for which the calculations are perfomed

            </doc>
            <doc lang="rus">
                Рассчитывает удельную теплоемкость поступательных степеней свободы

                **Параметры**:
                
                    * ``double T`` - температура

                    * ``const kappa::Particle &particle`` - частица, для которой производится расчет
            </doc>
        </fdoc>
        */

        double c_rot(double T, const kappa::Molecule &molecule, int i=0, int e=0);
        /*
        <fdoc>
            <doc lang="eng">
                Calculate the specific heat of rotational degrees of freedom

                **Parameters**:
                
                    * ``double T`` - the temperature

                    * ``const kappa::Molecule &molecule`` - the molecule for which the calculations are perfomed

                    * ``int i`` - the vibrational level of the molecule (default value is 0)

                    * ``int e`` - the electron level of the molecule (default value is 0)
            </doc>
            <doc lang="rus">
                Рассчитывает удельную теплоемкость вращательных степеней свободы

                **Параметры**:
                
                    * ``double T`` - температура

                    * ``const kappa::Molecule &molecule`` - молекула, для которой производится расчет

                    * ``int i`` - Колебательный уровень молекулы (значение по умолчанию 0)

                    * ``int e`` - Электронный уровень молекулы (значение по умолчанию 0)
            </doc>
        </fdoc>
        */

        double probability_VV(double rel_vel, const kappa::Molecule &molecule1, const kappa::Molecule &molecule2, const kappa::Interaction &interaction,
                              int i, int k, int delta_i, int e1, int e2, kappa::models_prob_vv model=kappa::models_prob_vv::model_prob_vv_fho);
        /*
        <fdoc>
            <i1><math>M_{1}(e1, i) + M_{2}(e2, k) \to M_{1}(e1, i+\Delta i) + M_{2}(e2, k-\Delta i)</math></i1>
            <doc lang="eng">
                Calculate the probability of a VV transition <insert>i1</insert>

                **Parameters**:
                
                    * ``double rel_vel`` - the relative velocity of the particles

                    * ``const kappa::Molecule &molecule1`` - the molecule <imath>M_{1}</imath>

                    * ``const kappa::Molecule &molecule2`` - the molecule <imath>M_{2}</imath>

                    * ``const kappa::Interaction &interaction`` - the interaction data for the colliding particles

                    * ``int i`` - the vibrational level of the molecule <imath>M_{1}</imath>

                    * ``int k`` - the vibrational level of the molecule <imath>M_{2}</imath>

                    * ``int delta_i`` - the change in the vibrational level of the molecule <imath>M_{1}</imath> (the quantity <imath>\Delta i</imath>)

                    * ``int e1`` - the electron level of the molecule <imath>M_{1}</imath>

                    * ``int e2`` - the electron level of the molecule <imath>M_{2}</imath>

                    * ``kappa::models_prob_vv model`` - the model to be used to calculate the probability (default value is ``kappa::models_prob_vv::model_prob_vv_fho``), possible values:

                        * ``kappa::models_prob_vv::model_prob_vv_fho`` - the FHO model
            </doc>
            <doc lang="rus">
                Рассчитывает вероятность VV перехода <insert>i1</insert>

                **Параметры**:
                
                    * ``double rel_vel`` - относительная скорость частиц

                    * ``const kappa::Molecule &molecule1`` - молекула <imath>M_{1}</imath>

                    * ``const kappa::Molecule &molecule2`` - молекула <imath>M_{2}</imath>

                    * ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся частиц

                    * ``int i`` - колебательный уровень молекулы <imath>M_{1}</imath>

                    * ``int k`` - колебательный уровень молекулы <imath>M_{2}</imath>

                    * ``int delta_i`` - изменение колебательного уровня молекулы <imath>M_{1}</imath> (величина <imath>\Delta i</imath>)

                    * ``int e1`` - электронный уровень молекулы <imath>M_{1}</imath>

                    * ``int e2`` - электронный уровень молекулы <imath>M_{2}</imath>

                    * ``kappa::models_prob_vv model`` - модель для расчета вероятности (значение по умолчанию ``kappa::models_prob_vv::model_prob_vv_fho``), возможные значения:

                        * ``kappa::models_prob_vv::model_prob_vv_fho`` - модель FHO
            </doc>
        </fdoc>
        */

        double probability_VT(double rel_vel, const kappa::Molecule &molecule, const kappa::Interaction &interaction, int i,  
                              int delta_i, int e, kappa::models_prob_vt model=kappa::models_prob_vt::model_prob_vt_fho);
        /*
        <fdoc>
            <i1><math>M(e, i) + P \to M(e, i+\Delta i) +P</math></i1>
            <doc lang="eng">
                Calculate the probability of a VT transition <insert>i1</insert>

                **Parameters**:
                
                    * ``double rel_vel`` - the relative velocity of the particles

                    * ``const kappa::Molecule &molecule`` - the molecule <imath>M</imath>

                    * ``const kappa::Interaction &interaction`` - the interaction data for the colliding particles

                    * ``int i`` - the vibrational level of the molecule <imath>M</imath>

                    * ``int delta_i`` - the change in the vibrational level of the molecule <imath>M</imath> (the quantity <imath>\Delta i</imath>)

                    * ``int e`` - the electron level of the molecule <imath>M</imath>

                    * ``kappa::models_prob_vt model`` - the model to be used to calculate the probability (default value is ``kappa::models_prob_vt::model_prob_vt_fho``), possible values:

                        * ``kappa::models_prob_vt::model_prob_vt_fho`` - the FHO model
            </doc>
            <doc lang="rus">
                Рассчитывает вероятность VT перехода <insert>i1</insert>

                **Параметры**:
                
                    * ``double rel_vel`` - относительная скорость частиц

                    * ``const kappa::Molecule &molecule`` - молекула <imath>M</imath>

                    * ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся частиц

                    * ``int i`` - колебательный уровень молекулы <imath>M</imath>

                    * ``int delta_i`` - изменение колебательного уровня молекулы <imath>M</imath> (величина <imath>\Delta i</imath>)

                    * ``int e`` - электронный уровень молекулы <imath>M</imath>

                    * ``kappa::models_prob_vt model`` - модель для расчета вероятности (значение по умолчанию ``kappa::models_prob_vt::model_prob_vt_fho``), возможные значения:

                        * ``kappa::models_prob_vt::model_prob_vt_fho`` - модель FHO
            </doc>
        </fdoc>
        */

        double probability_diss(double rel_vel, const kappa::Molecule &molecule, const kappa::Interaction &interaction, int i, int e,
                                kappa::models_prob_diss model=kappa::models_prob_diss::model_prob_diss_thresh_cmass_vibr);
        /*
        <fdoc>
            <i1><math>AB(e, i) + P \to A + B +P</math></i1>
            <doc lang="eng">
                Calculate the probability of a dissociation reaction <insert>i1</insert>

                **Parameters**:
                
                    * ``double rel_vel`` - the relative velocity of the particles

                    * ``const kappa::Molecule &molecule`` - the molecule <imath>AB</imath>

                    * ``const kappa::Interaction &interaction`` - the interaction data for the colliding particles

                    * ``int i`` - the vibrational level of the molecule <imath>AB</imath>

                    * ``int e`` - the electron level of the molecule <imath>AB</imath>

                    * ``kappa::models_prob_diss model`` - the model to be used to calculate the probability (default value is ``kappa::models_prob_diss::model_prob_diss_thresh_cmass_vibr``), possible values:

                        * ``kappa::models_prob_diss::model_prob_diss_thresh_cmass_vibr`` - threshold model accounting for the vibrational energy and the translational energy along the center-of-mass line

                        * ``kappa::models_prob_diss::model_prob_diss_thresh_vibr`` - threshold model accounting for the vibrational energy and the full translational energy

                        * ``kappa::models_prob_diss::model_prob_diss_thresh_cmass`` - threshold model accounting for the translational energy along the center-of-mass line

                        * ``kappa::models_prob_diss::model_prob_diss_thresh`` - threshold model accounting for the full translational energy
            </doc>
            <doc lang="rus">
                Рассчитывает вероятность реакции диссоциации <insert>i1</insert>

                **Параметры**:
                
                    * ``double rel_vel`` - относительная скорость частиц

                    * ``const kappa::Molecule &molecule`` - молекула <imath>AB</imath>

                    * ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся частиц

                    * ``int i`` - колебательный уровень молекулы <imath>AB</imath>

                    * ``int e`` - электронный уровень молекулы <imath>AB</imath>

                    * ``kappa::models_prob_diss model`` - модель для расчета вероятности (значение по умолчанию ``kappa::models_prob_diss::model_prob_diss_thresh_cmass_vibr``), возможные значения:

                        * ``kappa::models_prob_diss::model_prob_diss_thresh_cmass_vibr`` - пороговая модель, учитывающая колебательную энергию и поступательную энергию вдоль линии центра масс

                        * ``kappa::models_prob_diss::model_prob_diss_thresh_vibr`` - пороговая модель, учитывающая колебательную энергию и полную поступательную энергию

                        * ``kappa::models_prob_diss::model_prob_diss_thresh_cmass`` - пороговая модель, учитывающая поступательную энергию вдоль линии центра масс

                        * ``kappa::models_prob_diss::model_prob_diss_thresh`` - пороговая модель, учитывающая полную поступательную энергию
            </doc>
        </fdoc>
        */

        double probability_diss(double rel_vel, kappa::Molecule const &molecule, kappa::Interaction const &interaction, int i,
                                kappa::models_prob_diss model=kappa::models_prob_diss::model_prob_diss_thresh_cmass_vibr);
        /*
        <fdoc>
            <i1><math>AB(i) + P \to A + B +P</math></i1>

            <doc lang="rus">
                Рассчитывает вероятность реакции диссоциации (предполагается, что молекула <imath>AB</imath> находится в основном электронном состоянии) (см. документацию предыдущей функции, параметр ``e=0``) <insert>i1</insert> 
            </doc>
        </fdoc>
        */

        double crosssection_elastic(double rel_vel, kappa::Interaction const &interaction, kappa::models_cs_elastic model=kappa::models_cs_elastic::model_cs_el_vss);
        /*
        <fdoc>
            <doc lang="eng">
                Calculate the elastic collision cross section

                **Parameters**:
                
                    * ``double rel_vel`` - the relative velocity of the particles

                    * ``const kappa::Interaction &interaction`` - the interaction data for the colliding particles

                    * ``kappa::models_cs_elastic model`` - the model to be used to calculate the cross section (default value is ``model_cs_el_vss``), possible values:

                        * ``model_cs_el_rs`` - the rigid sphere (RS) model

                        * ``model_cs_el_vss`` - the variable soft sphere (VSS) model
            </doc>
            <doc lang="rus">
                Рассчитывает упругое сечение столкновения

                **Параметры**:
                
                    * ``double rel_vel`` - относительная скорость частиц

                    * ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся частиц

                    * ``kappa::models_prob_diss model`` - модель для расчета сечения (значение по умолчанию ``kappa::models_cs_elastic::model_cs_el_vss``), возможные значения:

                        * ``kappa::models_cs_elastic::model_cs_el_rs`` - модель твердых сфер (RS)

                        * ``kappa::models_cs_elastic::model_cs_el_vss`` - модель переменных мягких сфер (VSS)
            </doc>
        </fdoc>
        */

        double crosssection_VT(double rel_vel, kappa::Molecule const &molecule, kappa::Interaction const &interaction,
                               int i, int delta_i, int e, kappa::models_cs_vt model=kappa::models_cs_vt::model_cs_vt_vss_fho);
		/*
		<fdoc>
            <i1><math>M(e, i) + P \to M(e, i+\Delta i) +P</math></i1>
			<doc lang="rus">	
				Расчет сечения VT перехода <insert>i1</insert>
				
				**Параметры**:
					
					* ``double rel_vel`` - относительная скорость частиц
					
					* ``const kappa::Molecule &molecule`` - молекула для которой производится расчет
					
					* ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся частиц
					
					* ``int i`` - колебательный уровень молекулы 
					
					* ``int delta_i - изменение колебательного уровня молекулы 
					
					* ``int e`` - электронный уровень молекулы
					
					* ``kappa::models_cs_vt model`` - модель для расчета сечения VT-перехода (значение по умолчанию ``kappa::models_cs_vt::model_cs_vt_vss_fho``), возможные значения:
						
						* ``kappa::models_cs_vt::model_cs_vt_rs_fho`` - модель твердых сфер (RS) в сочетании с моделью FHO 
						
						* ``kappa::models_cs_vt::model_cs_vt_vss_fho`` - модель сфер переменного диаметра (VSS) в сочетании с моделью FHO 
			</doc>
		</fdoc>
		*/
        

        double crosssection_VT(double rel_vel, kappa::Molecule const &molecule, kappa::Interaction const &interaction,
                                int i, int delta_i, kappa::models_cs_vt model=kappa::models_cs_vt::model_cs_vt_vss_fho);
		/*
		<fdoc>
            <i1><math>M(i) + P \to M(i+\Delta i) +P</math></i1> 
			<doc lang="rus">
				Расчет сечения VT перехода (предполагается, что молекула <imath>AB</imath> находится в основном электронном состоянии) (см. документацию предыдущей функции, параметр ``e=0``) <insert>i1</insert> 
			</doc>
		</fdoc>
		*/

        double crosssection_diss(double rel_vel, kappa::Molecule const &molecule, kappa::Interaction const &interaction, int i, int e, kappa::models_cs_diss model=kappa::models_cs_diss::model_cs_diss_vss_thresh_cmass_vibr);
		/*
        <fdoc>
            <i1><math>AB(e, i) + P \to A + B +P</math></i1>
            <doc lang="eng">
                Calculate the cross section of a dissociation reaction <insert>i1</insert>

                **Parameters**:
                
                    * ``double rel_vel`` - the relative velocity of the particles

                    * ``const kappa::Molecule &molecule`` - the molecule <imath>AB</imath>

                    * ``const kappa::Interaction &interaction`` - the interaction data for the colliding particles

                    * ``int i`` - the vibrational level of the molecule <imath>AB</imath>

                    * ``int e`` - the electron level of the molecule <imath>AB</imath>

                    * ``kappa::models_cs_diss model`` - the model to be used to calculate the cross section (default value is ``kappa::models_cs_diss::model_cs_diss_vss_thresh_cmass_vibr``), possible values:

                        * ``kappa::models_cs_diss::model_cs_diss_rs_thresh_cmass_vibr`` - rigid sphere (RS) model coupled with a threshold model accounting for the vibrational energy
                          and the translational energy along the center-of-mass line

                        * ``kappa::models_cs_diss::model_cs_diss_rs_thresh_vibr`` - rigid sphere (RS) model coupled with a threshold model accounting for the vibrational energy and the full translational energy

                        * ``kappa::models_cs_diss::model_cs_diss_rs_thresh_cmass`` - rigid sphere (RS) model coupled with a threshold model accounting for the translational energy along the center-of-mass line

                        * ``kappa::models_cs_diss::model_cs_diss_rs_thresh`` - rigid sphere (RS) model coupled with a threshold model accounting for the full translational energy

                        * ``kappa::models_cs_diss::model_cs_diss_vss_thresh_cmass_vibr`` - variable soft sphere (VSS) model coupled with a threshold model accounting for the vibrational energy
                          and the translational energy along the center-of-mass line

                        * ``kappa::models_cs_diss::model_cs_diss_vss_thresh_vibr`` - variable soft sphere (VSS) model coupled with a threshold model accounting for the vibrational energy and the full translational energy

                        * ``kappa::models_cs_diss::model_cs_diss_vss_thresh_cmass`` - variable soft sphere (VSS) model coupled with a threshold model accounting for the translational energy along the center-of-mass line

                        * ``kappa::models_cs_diss::model_cs_diss_vss_thresh`` - variable soft sphere (VSS) model coupled with a threshold model accounting for the full translational energy

                        * ``kappa::models_cs_diss::model_cs_diss_ilt`` - a model based on the inverse Laplace transform of QCT calculations in Phys4Entry database

            </doc>
            <doc lang="rus">
                Рассчитывает сечение реакции диссоциации <insert>i1</insert>

                **Параметры**:
                
                    * ``double rel_vel`` - относительная скорость частиц

                    * ``const kappa::Molecule &molecule`` - молекула <imath>AB</imath>

                    * ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся частиц

                    * ``int i`` - колебательный уровень молекулы <imath>AB</imath>

                    * ``int e`` - электронный уровень молекулы <imath>AB</imath>

                    * ``kappa::models_cs_diss model`` - модель для расчета сечения (значение по умолчанию ``model_cs_diss_vss_thresh_cmass_vibr``), возможные значения:

                        * ``kappa::models_cs_diss::model_cs_diss_rs_thresh_cmass_vibr`` - модель твердых сфер (RS) в сочетании с пороговой моделью, учитывающей колебательную энергию и поступательную энергию вдоль линии центра масс

                        * ``kappa::models_cs_diss::model_cs_diss_rs_thresh_vibr`` - модель твердых сфер (RS) в сочетании с пороговой модель, учитывающей колебательную энергию и полную поступательную энергию

                        * ``kappa::models_cs_diss::model_cs_diss_rs_thresh_cmass`` - модель твердых сфер (RS) в сочетании с пороговой моделью, учитывающей поступательную энергию вдоль линии центра масс

                        * ``kappa::models_cs_diss::model_cs_diss_rs_thresh`` - модель твердых сфер (RS) в сочетании с пороговой моделью, учитывающей полную поступательную энергию

                        * ``kappa::models_cs_diss::model_cs_diss_vss_thresh_cmass_vibr`` - модель переменных мягких сфер (VSS) в сочетании с пороговой моделью,
                          учитывающей колебательную энергию и поступательную энергию вдоль линии центра масс

                        * ``kappa::models_cs_diss::model_cs_diss_vss_thresh_vibr`` - модель переменных мягких сфер (VSS) в сочетании с пороговой модель, учитывающей колебательную энергию и полную поступательную энергию

                        * ``kappa::models_cs_diss::model_cs_diss_vss_thresh_cmass`` - модель переменных мягких сфер (VSS) в сочетании с пороговой моделью, учитывающей поступательную энергию вдоль линии центра масс

                        * ``kappa::models_cs_diss::model_cs_diss_vss_thresh`` - модель переменных мягких сфер (VSS) в сочетании с пороговой моделью, учитывающей полную поступательную энергию

                        * ``kappa::models_cs_diss::model_cs_diss_ilt`` - модель, основанная на преобразовании Лапласа от резултьтатов квази-классических траекторных расчетов в базе данных Phys4Entry 

            </doc>
        </fdoc>
        */

        double crosssection_diss(double rel_vel, kappa::Molecule const &molecule, kappa::Interaction const &interaction, int i, kappa::models_cs_diss model=kappa::models_cs_diss::model_cs_diss_vss_thresh_cmass_vibr);
        /*
        /*
        <fdoc>
            <i1><math>AB(i) + P \to A + B + P</math></i1>

            <doc lang="rus">
                Рассчитывает сечения реакции диссоциации (предполагается, что молекула <imath>AB</imath> находится в основном электронном состоянии) (см. документацию предыдущей функции, параметр ``e=0``) <insert>i1</insert> 
            </doc>
        </fdoc>
        */

        double Z_coll(double T, double n, kappa::Interaction const &interaction);
		/*
		<fdoc>
            <doc lang="rus">
                Число столкновений в единицу времени
                
                **Параметры**:
                    
                    * ``double T`` - температура

                    * ``double n`` - числовая плотность

                    * ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся частиц
            </doc>
		</fdoc>
		*/

        double rot_collision_number_parker(double T, kappa::Molecule const &molecule, kappa::Interaction const &interaction);
       /*
	   <fdoc>
			<doc lang="rus">
				Число столкновений, необходимое для установления равновесия по вращательным степеням свободы, вычисленное по формуле Паркера
				
				**Параметры**:
					
					* ``double T`` - температура

					* ``const kappa::Molecule &molecule`` - молекула, для которой производится расчет

					* ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся частиц
			</doc>
		</fdoc>
		*/

		double rot_relaxation_time_parker(double T, double n, kappa::Molecule const &molecule, kappa::Interaction const &interaction, kappa::models_omega model=kappa::models_omega::model_omega_esa);
        /*
        <fdoc>
            <doc lang="rus">
                Время вращательной релаксации, вычисленное по формуле Паркера
                
                **Параметры**:
                    
                    * ``double T`` - температура

                    * ``double n`` - числовая плотность

                    * ``const kappa::Molecule &molecule`` - молекула, для которой производится расчет

                    * ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся частиц

                    * ``kappa::models_omega model`` - модель для вычисления Омега-интеграла (значение модели по умолчанию ``kappa::models_omega::model_omega_esa``):
                        
                        * ``kappa::models_omega::model_omega_rs`` - модель твердых сфер (RS)
                        
                        * ``kappa::models_omega::model_omega_vss`` - модель сфер переменного диаметра (VSS)
                        
                        * ``kappa::models_omega::model_omega_bornmayer`` - модель Борна–Майера
                        
                        * ``kappa::models_omega::model_omega_lennardjones`` - модель Леннарда–Джонса
                        
                        * ``kappa::models_omega::model_omega_esa`` - модель ESA (феноменологическая модель)
            </doc>
        </fdoc>
        */

		double vibr_relaxation_time_MW(double T, double n, kappa::Molecule const &molecule, kappa::Interaction const &interaction);
		/*
		<fdoc>
			<doc lang="rus">
				Время колебательной релаксации, вычисленное по формуле Милликена-Уайта
				
				**Параметры**:
					
					* ``double T`` - температура

					* ``double n`` - числовая плотность

					* ``const kappa::Molecule &molecule`` - молекула, для которой производится расчет

					* ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся частиц
			</doc>
		</fdoc>
		*/

		double vibr_relaxation_time_Park_corr(double T, double n, kappa::Interaction const &interaction, double crosssection=1e-10);
		/*
		<fdoc>
			<doc lang="rus">
                Поправка Парка к формуле Милликена-Уайта для вычисления времени колебательной релаксации 
                
                **Параметры**:
                    
                    * ``double T`` - температура

                    * ``double n`` - числовая плотность

                    * ``const kappa::Molecule &molecule`` - молекула, для которой производится расчет

                    * ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся частиц

                    * ``double crosssection`` - сечение для расчета поправки Парка (значение по умолчанию ``1e-10``)
			</doc>
		</fdoc>
		*/

        double omega_integral(double T, kappa::Interaction const &interaction, int l, int r, kappa::models_omega model=kappa::models_omega::model_omega_esa, bool dimensional = true);
		/*
		<fdoc>
			<doc lang="rus">
				Вычисление Омега-интегралов для взаимодействий нейтральных частиц, взаимодействий ионов и нейтральных частиц, для взаимодействия электронов и нейтральных частиц
				
				**Параметры**:
					
					* ``double T`` - температура

					* ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся частиц

					* ``int l`` - степень интеграла

					* ``int r`` - степень интеграла

					* ``bool dimensional`` - тип возвращаемого интеграла (размерный/безразмерный) (значение по умолчанию ``true``)

					* ``kappa::models_omega model`` - модель для вычисления Омега-интеграла (значение модели по умолчанию ``kappa::models_omega::model_omega_esa``), возможные значения:
						
						* ``kappa::models_omega::model_omega_rs`` - модель твердых сфер (RS)
						
						* ``kappa::models_omega::model_omega_vss`` - модель сфер переменного диаметра (VSS)
						
						* ``kappa::models_omega::model_omega_bornmayer`` - модель Борна–Майера
						
						* ``kappa::models_omega::model_omega_lennardjones`` - модель Леннарда–Джонса
						
						* ``kappa::models_omega::model_omega_esa`` - модель ESA (феноменологическая модель)
			</doc>
		</fdoc>
		*/

        double omega_integral(double T, kappa::Interaction const &interaction, int l, int r, double debye_length, kappa::models_omega model=kappa::models_omega::model_omega_esa, bool dimensional = true);
        /*
        <fdoc>
            <doc lang="rus">
                Вычисление Омега-интегралов для взаимодействий произвольных частиц
                
                **Параметры**:
                    
                    * ``double T`` - температура

                    * ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся частиц

                    * ``int l`` - степень интеграла

                    * ``int r`` - степень интеграла

                    * ``double debye_length` - Дебаевская длина (используется при расчете интегралов для взаимодействий заряженных частиц)

                    * ``bool dimensional`` - тип возвращаемого интеграла (размерный/безразмерный) (значение по умолчанию ``true``)

                    * ``kappa::models_omega model`` - модель для вычисления Омега-интеграла (значение модели по умолчанию ``kappa::models_omega::model_omega_esa``),
                      возможные значения (для столкновений заряженных частиц данный параметр ни на что не влияет):
                        
                        * ``kappa::models_omega::model_omega_rs`` - модель твердых сфер (RS)
                        
                        * ``kappa::models_omega::model_omega_vss`` - модель сфер переменного диаметра (VSS)
                        
                        * ``kappa::models_omega::model_omega_bornmayer`` - модель Борна–Майера
                        
                        * ``kappa::models_omega::model_omega_lennardjones`` - модель Леннарда–Джонса
                        
                        * ``kappa::models_omega::model_omega_esa`` - модель ESA (феноменологическая модель)
            </doc>
        </fdoc>
        */

        double integral_VT(double T, int degree, kappa::Molecule const &molecule, kappa::Interaction const &interaction, int i, int delta_i, int e, kappa::models_cs_vt model = kappa::models_cs_vt::model_cs_vt_vss_fho);
        /*
		<fdoc>
            <i1><math>M(e, i) + P \to M(e, i+\Delta i) + P</math></i1>
			<doc lang="rus">
				Вычисление интегралов по сечению VT-перехода <insert>i1</insert>

				**Параметры**:
					
					* ``double T`` - температура
					
					* ``int degree`` -  степень интеграла
					
					* ``const kappa::Molecule &molecule`` - молекула, для которой производится расчет
					
					* ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся частиц
					
					* ``int i`` - колебательный уровень молекулы
					
					* ``int delta_i`` - изменение колебательного уровня молекулы
					
					* ``int e`` -электронный уровень молекулы
					
					* ``kappa::models_cs_vt model`` - модель для вычисления интегралов по сечению VT-перехода (значение модели по умолчанию ``model_cs_vt_vss_fho``), возможные значения:
						
						* ``kappa::models_cs_vt::model_cs_vt_rs_fho`` - модель твердых сфер (RS) в сочетании с моделью FHO 
						
						* ``kappa::models_cs_vt::model_cs_vt_vss_fho`` - модель сфер переменного диаметра (VSS) в сочетании с моделью FHO 
			</doc>
		</fdoc>
		*/

		double integral_VT(double T, int degree, kappa::Molecule const &molecule, kappa::Interaction const &interaction, int i, int delta_i, kappa::models_cs_vt model = kappa::models_cs_vt::model_cs_vt_vss_fho);
		/*
		<fdoc>
            <i1><math>M(i) + P \to M(i+\Delta i) + P</math></i1>
			<doc lang="rus">
				Вычисление интегралов по сечению VT-перехода (предполагается, что молекула <imath>M</imath> находится в основном электронном состоянии) (см. документацию предыдущей функции, параметр ``e=0``) <insert>i1</insert>
			</doc>
		</fdoc>
		*/

		double k_VT(double T, kappa::Molecule const &molecule, kappa::Interaction const &interaction, int i, int delta_i, int e, kappa::models_k_vt model = kappa::models_k_vt::model_k_vt_vss_fho);
        /*
		<fdoc>
            <i1><math>M(e, i) + P \to M(e, i+\Delta i) + P</math></i1>
			<doc lang="rus">
				Вычисление коэффициента скорости VT-перехода <insert>i1</insert>
				
				**Параметры**:
					
					* ``double T`` - температура
					
					* ``const kappa::Molecule &molecule`` - молекула, для которой производится расчет
					
					* ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся частиц
					
					* ``int i`` - колебательный уровень молекулы
					
					* ``int delta_i`` - изменения колебательного уровня молекулы
					
					* ``int e`` -  электронный уровень молекулы
					
					* ``kappa::models_k_vt model`` - модель для расчета коэффициента скорости (значение модели по умолчанию ``model_k_vt_vss_fho``), возможные значения:
					
						* ``kappa::models_k_vt::model_k_vt_rs_fho`` - модель твердых сфер (RS) в сочетании с моделью FHO 
						
						* ``kappa::models_k_vt::model_k_vt_vss_fho`` - модель переменных мягких сфер (VSS) в сочетании с моделью FHO 
						
						* ``kappa::models_k_vt::model_k_vt_ssh`` - модель SSH 
						
						* ``kappa::models_k_vt::model_k_vt_phys4entry`` -  аппроксимация из базы данных Phys4Entry
						
						* ``kappa::models_k_vt::model_k_vt_billing`` - аппроксимация траекторных расчетов Биллинга 
			</doc>
		</fdoc>
		*/

		double k_VT(double T, kappa::Molecule const &molecule, kappa::Interaction const &interaction, int i, int delta_i, kappa::models_k_vt model = kappa::models_k_vt::model_k_vt_vss_fho);
		/*
		<fdoc>
            <i1><math>M(i) + P \to M(i+\Delta i) + P</math></i1>
			<doc lang="rus">
				Вычисление коэффициента скорости VT-перехода (предполагается, что молекула <imath>M</imath> находится в основном электронном состоянии) (см. документацию предыдущей функции, параметр ``e=0``) <insert>i1</insert>
			</doc>
		</fdoc>
		*/

        double k_exch(double T, kappa::Molecule const &molecule, kappa::Atom const &atom, kappa::Interaction const &interaction,
                      int i, int e, int num_electron_levels=-1, kappa::models_k_exch model=kappa::models_k_exch::model_k_exch_warnatz);
		/*
        <fdoc>
            <doc lang="rus">
                Вычисление коэффициентов скорости обменных реакций (для молекул с электронным возбуждением) без учета колебательного и электронного возбуждения молекулы-продукта реакции
                
                **Параметры**:
                    
                    * ``double T`` - температура
                    
                    * ``const kappa::Molecule &molecule`` - молекула, для которой производится расчет
                    
                    * ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся частиц
                    
                    * ``kappa::Atom const &atom`` - рассматриваемый атом
                    
                    * ``int i`` - колебательный уровень исходной молекулы
                    
                    * ``int k`` - колебательный уровень образующейся молекулы 
                    
                    * ``int e`` - электронный уровень молекулы
                    
                    * ``int num_electron_level`` - число учитываемых электронных уровней молекулы-реагента (значение по умолчанию =-1 - учитываются все имеющиеся уровни)
                    
                    * ``kappa::models_k_exch model`` - модель для расчета коэффициента скорости обменной реакции (значение модели по умолчанию ``model_k_exch_warnatz``), возможные значения:
                        
                        * ``kappa::models_k_exch::model_k_exch_arrh_scanlon - модель Аррениуса в сочетании с данными из Scanlon et al.
                        
                        * ``kappa::models_k_exch::model_k_exch_arrh_park`` - модель Аррениуса в сочетании с данными из Park et al.
                        
                        * ``kappa::models_k_exch::model_k_exch_warnatz`` - модель Варнатца
                        
                        * ``kappa::models_k_exch::model_k_exch_rf`` - модель Русанова-Фридмана
                        
                        * ``kappa::models_k_exch::model_k_exch_polak`` - модель Полака
            </doc>
		</fdoc>
		*/

        double k_exch(double T, kappa::Molecule const &molecule, kappa::Atom const &atom, kappa::Interaction const &interaction,
                      int i, kappa::models_k_exch model=kappa::models_k_exch::model_k_exch_warnatz);
		/*
		<fdoc>
			<doc lang="rus">
                Вычисление коэффициентов скорости обменных реакций (для молекул без электронного возбуждения) без учета колебательного и электронного возбуждения молекулы-продукта реакции
                (см. документацию предыдущей функции, параметр ``e=0``)
			</doc>
		</fdoc>
		*/

        double k_exch(double T, kappa::Molecule const &molecule, kappa::Atom const &atom, kappa::Molecule const &molecule_prod, kappa::Interaction const &interaction,
                      int i, int k, int e, int num_electron_levels=-1, kappa::models_k_exch model=kappa::models_k_exch::model_k_exch_maliat_infty_arrh_scanlon);
        /*
        <fdoc>
            <doc lang="rus">
                Вычисление коэффициентов скорости обменных реакций с учетом колебательного возбуждения молекулы-продукта реакции (для молекул с электронным возбуждением). В случае, если используется не модифицированная модель Алиата,
                колебательное возбуждение продукта реакции учитываться не будет.
                
                **Параметры**:
                
                    * ``double T`` - температура
                
                    * ``const kappa::Molecule &molecule`` - молекула, для которой производится расчет
                
                    * ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся частиц
                
                    * ``kappa::Atom const &atom`` - рассматриваемый атом
                    
                    * ``int i`` - колебательный уровень исходной молекулы
                    
                    * ``int k`` - колебательный уровень образующейся молекулы 
                
                    * ``int e`` - электронный уровень молекулы
                
                    * ``int num_electron_levels`` - число учитываемых электронных уровней молекулы-реагента (значение по умолчанию =-1 - учитываются все имеющиеся уровни)
                
                    * ``kappa::models_k_exch model`` - модель для расчета коэффициента скорости обменной реакции (значение модели по умолчанию ``model_k_exch_maliat_infty_arrh_scanlon``), возможные значения:
                
                        * ``kappa::models_k_exch::model_k_exch_arrh_scanlon - модель Аррениуса с данными из Scanlon et al. (не учитывает электронное возбуждение молекулы-реагента и колебательное возбуждение молекулы-продукта)
                        
                        * ``kappa::models_k_exch::model_k_exch_arrh_park`` - модель Аррениуса с данными из Park et al. (не учитывает электронное возбуждение молекулы-реагента и колебательное возбуждение молекулы-продукта)
                        
                        * ``kappa::models_k_exch::model_k_exch_warnatz`` - модель Варнатца (не учитывает электронное возбуждение молекулы-реагента и колебательное возбуждение молекулы-продукта)
                        
                        * ``kappa::models_k_exch::model_k_exch_rf`` - модель Русанова-Фридмана (не учитывает электронное возбуждение молекулы-реагента и колебательное возбуждение молекулы-продукта)
                        
                        * ``kappa::models_k_exch::model_k_exch_polak`` - модель Полака (не учитывает электронное возбуждение молекулы-реагента и колебательное возбуждение молекулы-продукта)
                
                        * ``model_k_exch_maliat_D6k_arrh_scanlon`` - модифицированная модель Алиэта с данными из Scanlon et al., <imath>U=D/6k</imath>
                
                        * ``model_k_exch_maliat_3T_arrh_scanlon`` - модифицированная модель Алиэта с данными из Scanlon et al., <imath>U=3Tk</imath>
                
                        * ``model_k_exch_maliat_infty_arrh_scanlon`` - модифицированная модель Алиэта с данными из Scanlon et al., <imath>U=\infty</imath>
                
                        * ``model_k_exch_maliat_D6k_arrh_park`` - модифицированная модель Алиэта с данными из Park et al., <imath>U=D/6k</imath>
                
                        * ``model_k_exch_maliat_3T_arrh_park`` - модифицированная модель Алиэта с данными из Park et al., <imath>U=3Tk</imath>
                
                        * ``model_k_exch_maliat_infty_arrh_park`` - модифицированная модель Алиэта с данными из Park et al., <imath>U=\infty</imath>
            </doc>
        </fdoc>
        */

        double k_exch(double T, kappa::Molecule const &molecule, kappa::Atom const &atom, kappa::Molecule const &molecule_prod, kappa::Interaction const &interaction,
                      int i, int k, kappa::models_k_exch model=kappa::models_k_exch::model_k_exch_maliat_infty_arrh_scanlon);
        /*
        <fdoc>
            <doc lang="rus">
                Вычисление коэффициентов скорости обменных реакций с учетом колебательного возбуждения молекулы-продукта реакции (для молекул без электронного возбуждения) (см. документацию предыдущей функции,
                параметры ``e=0``, ``num_electron_levels=1``)
            </doc>
        </fdoc>
        */

        double integral_diss(double T, int degree, kappa::Molecule const &molecule, kappa::Interaction const &interaction, int i, int e, kappa::models_cs_diss model=kappa::models_cs_diss::model_cs_diss_vss_thresh_cmass_vibr);
		/*
		<fdoc>
            <i1><math>AB(e, i) + P \to A + B + P</math></i1>
			<doc lang="rus">
				Расчета интеграла по сечению реакции диссоциации <insert>i1</insert>
                
				**Параметры**:
					
					* ``double T`` - температура
					
					* ``int degree`` - степень интеграла
					
					* ``const kappa::Molecule &molecule`` - молекула, для которой производится расчет
					
					* ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся молекул
					
					* ``int i`` - колебательный уровень молекулы
					
					* ``int e`` - электронный уровень молекулы
					
					* ``kappa::models_cs_diss model`` - модель сечения (значение по умолчанию = ``kappa::models_cs_diss::model_cs_diss_vss_thresh_cmass_vibr``), возможные значение:
						
						* ``kappa::models_cs_diss::model_cs_diss_rs_thresh_cmass_vibr``- модель твердых сфер (RS) в сочетании с пороговой моделью, учитывающей колебательную энергию и поступательную энергию вдоль линии центра масс
						
						* ``kappa::models_cs_diss::model_cs_diss_rs_thresh_vibr`` - модель твердых сфер (RS) в сочетании с пороговой модель, учитывающей колебательную энергию и полную поступательную энергию
						
						* ``kappa::models_cs_diss::model_cs_diss_rs_thresh_cmass`` - модель твердых сфер (RS) в сочетании с пороговой моделью, учитывающей поступательную энергию вдоль линии центра масс
						
						* ``kappa::models_cs_diss::model_cs_diss_rs_thresh`` -  модель твердых сфер (RS) в сочетании с пороговой моделью, учитывающей полную поступательную энергию
					    
						* ``kappa::models_cs_diss::model_cs_diss_vss_thresh_cmass_vibr`` -  модель переменных мягких сфер (VSS) в сочетании с пороговой моделью, учитывающей колебательную энергию
                          и поступательную энергию вдоль линии центра масс
						
						* ``kappa::models_cs_diss::model_cs_diss_vss_thresh_vibr`` -  модель переменных мягких сфер (VSS) в сочетании с пороговой модель, учитывающей колебательную энергию и полную поступательную энергию
						
						* ``kappa::models_cs_diss::model_cs_diss_vss_thresh_cmass`` -  модель переменных мягких сфер (VSS) в сочетании с пороговой моделью, учитывающей поступательную энергию вдоль линии центра масс
						
						* ``kappa::models_cs_diss::model_cs_diss_vss_thresh`` - модель переменных мягких сфер (VSS) в сочетании с пороговой моделью, учитывающей полную поступательную энергию
						
						* ``kappa::models_cs_diss::model_cs_diss_ilt`` -  модель, основанная на обратном преобразовании Лапласа от результатов квази-классических траекторных расчетов в базе данных Phys4Entry
			</doc>
		</fdoc>
		*/

        double integral_diss(double T, int degree, kappa::Molecule const &molecule, kappa::Interaction const &interaction, int i, kappa::models_cs_diss model=kappa::models_cs_diss::model_cs_diss_vss_thresh_cmass_vibr);
		/*
		<fdoc>
            <i1><math>AB(i) + P \to A + B + P</math></i1>
			<doc lang="rus">
				Расчет интеграла по сечению реакции диссоциации (предполагается, что молекула <imath>AB</imath> находится в основном электронном состоянии) (см. документацию предыдущей функции, параметр ``e=0``) <insert>i1</insert>
			</doc>
		</fdoc>
		*/

		double k_diss(double T, kappa::Molecule const &molecule, kappa::Interaction const &interaction, int i, int e, int num_electron_levels=-1, kappa::models_k_diss model=kappa::models_k_diss::model_k_diss_tm_3T_arrh_scanlon);
		/*
		<fdoc>
            <i1><math>AB(e, i) + P \to A + B + P</math></i1>
			<doc lang="rus">
				Вычисление коэффициента скорости диссоциации <insert>i1</insert>
				
				**Параметры**:
					
					* ``double T`` - температура
					
					* ``const kappa::Molecule &molecule`` - молекула, для которой производится расчет
					
					* ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся молекул
					
					* ``int i`` - колебательный уровень молекулы
					
					* ``int e`` - электронный уровень молекулы
					
					* ``int num_electron_levels`` - число учитываемых электронных уровней (значение по умолчанию =-1 - учитываются все имеющиеся уровни)
					
					* ``kappa::models_k_diss model`` - модель для вычисление коэффициента скорости диссоциации (значение по умолчанию = ``model_k_diss_tm_3T_arrh_scanlon``), возможные значения:
						
						* ``kappa::models_k_diss::model_k_diss_rs_thresh_cmass_vibr``- модель твердых сфер (RS) в сочетании с пороговой моделью, учитывающей колебательную энергию и поступательную энергию вдоль линии центра масс
						
						* ``kappa::models_k_diss::model_k_diss_rs_thresh_vibr`` - модель твердых сфер (RS) в сочетании с пороговой модель, учитывающей колебательную энергию и полную поступательную энергию
						
						* ``kappa::models_k_diss::model_k_diss_rs_thresh_cmass`` - модель твердых сфер (RS) в сочетании с пороговой моделью, учитывающей поступательную энергию вдоль линии центра масс
						
						* ``kappa::models_k_diss::model_k_diss_rs_thresh`` - модель твердых сфер (RS) в сочетании с пороговой моделью, учитывающей полную поступательную энергию
					    
						* ``kappa::models_k_diss::model_k_diss_vss_thresh_cmass_vibr`` - модель переменных мягких сфер (VSS) в сочетании с пороговой моделью, учитывающей колебательную энергию и поступательную энергию
                          вдоль линии центра масс
						
						* ``kappa::models_k_diss::model_k_diss_vss_thresh_vibr`` - модель переменных мягких сфер (VSS) в сочетании с пороговой модель, учитывающей колебательную энергию и полную поступательную энергию
						
						* ``kappa::models_k_diss::model_k_diss_vss_thresh_cmass`` - модель переменных мягких сфер (VSS) в сочетании с пороговой моделью, учитывающей поступательную энергию вдоль линии центра масс
						
						* ``kappa::models_k_diss::model_k_diss_vss_thresh`` - модель переменных мягких сфер (VSS) в сочетании с пороговой моделью, учитывающей полную поступательную энергию
					    
						* ``kappa::models_k_diss::model_k_diss_arrh_scanlon`` - модель Аррениуса с данными из Scanlon et al.
						
						* ``kappa::models_k_diss::model_k_diss_arrh_park`` - модель Аррениуса с данными из Park et al.
						
						* ``kappa::models_k_diss::model_k_diss_tm_D6k_arrh_scanlon`` - модель Тринора-Маррона с данными из Scanlon et al., <imath>U=D/6k</imath>
						
						* ``kappa::models_k_diss::model_k_diss_tm_3T_arrh_scanlon`` - модель Тринора-Маррона с данными из Scanlon et al., <imath>U=3T</imath>
						
						* ``kappa::models_k_diss::model_k_diss_tm_infty_arrh_scanlon`` - модель Тринора-Маррона с данными из Scanlon et al., <imath>U=\infty</imath>
						
						* ``kappa::models_k_diss::model_k_diss_tm_D6k_arrh_park`` - модель Тринора-Маррона с данными из Park et al., <imath>U=D/6k</imath>
						
						* ``kappa::models_k_diss::model_k_diss_tm_3T_arrh_park`` - модель Тринора-Маррона с данными из Park et al., <imath>U=3T</imath>
						
						* ``kappa::models_k_diss::model_k_diss_tm_infty_arrh_park`` - модель Тринора-Маррона с данными из Park et al., <imath>U=\infty</imath>
					    
						* ``kappa::models_k_diss::model_k_diss_phys4entry`` - аппроксимация из базы данных Phys4Entry
						
						* ``kappa::models_k_diss::model_k_diss_ilt``- модель, основанная на применении обратного преобразовании Лапласа к результатам квази-классических траекторных расчетов в базе данных Phys4Entry
			</doc>
		</fdoc>
		*/

        double k_diss(double T, kappa::Molecule const &molecule, kappa::Interaction const &interaction, int i, kappa::models_k_diss model=kappa::models_k_diss::model_k_diss_tm_3T_arrh_scanlon);
		/*
		<fdoc>
            <i1><math>AB(i) + P \to A + B + P</math></i1>
            <doc lang="rus">
                Вычисление коэффициента скорости диссоциации (предполагается, что молекула <imath>AB</imath> находится в основном электронном состоянии)
                (см. документацию предыдущей функции, параметры ``e=0``, ``num_electron_levels=1``) <insert>i1</insert>
			</doc>
		</fdoc>
		*/


        arma::vec Boltzmann_distribution(double T, double n, kappa::Molecule const &molecule, int e=0);
		/*
		<fdoc>
			<doc lang="rus">
				Возвращает вектор заселенности колебательных уровней в соответствии с больцмановским распределением
                
				**Параметры**:
                
    				* ``double T`` - температура
                
    				* ``double n`` - числовая плотность
                
    				* ``const kappa::Molecule &molecule`` - молекула, для которой производится расчёт
                
    				* ``int e`` - электронный уровень молекулы (значение по умолчанию = 0)
			</doc>
		</fdoc>
		*/

        double k_bf_VT(double T, kappa::Molecule const &molecule, int i, int delta_i, int e=0); // ratio of backward to forward (k(i+delta_i->i) / k(i->i+delta_i))
        /*
        <fdoc>
            <i1><math>k_{i+\Delta i \to i} / k_{i \to i+\Delta i}</math></i1>
            <doc lang="rus">
                Отношение обратного и прямого коэффициентов скорости VT-перехода

                **Параметры**:
                    
                    * ``double T`` - температура

                    * ``kappa::Molecule const &molecule`` - молекула, для которой производится расчёт

                    * ``int i`` - колебательный уровень молекулы 

                    * ``int delta_i`` - изменение колебательного уровня молекулы

                    * ``int e`` - электронный уровень молекулы (значение по умолчанию = 0)
            </doc>
        </fdoc> 
        */

        double k_bf_VV(double T, kappa::Molecule const &molecule1, kappa::Molecule const &molecule2, int i, int k, int delta_i, int e1=0, int e2=0);
        /*
        <fdoc>
            <doc lang="rus">
                Отношение прямого и обратного коэффициентов скорости VV-перехода

                **Параметры**:
                    
                    * ``double T`` - температура                      
                    
                    * ``const kappa::Molecule &molecule1`` - молекула, для которой производится расчёт, в начальном состоянии 

                    * ``const kappa::Molecule &molecule2`` - молекула, для которой производится расчет, в конечном состоянии

                    * ``int i`` - колебательный уровень молекулы в начальном состоянии

                    * ``int k`` - колебательный уровень молекулы в конечном состоянии

                    * ``int delta_i`` - изменение колебательного уровня молекулы 

                    * ``int e1`` - электронный уровень молекулы в начальном состоянии (значение по умолчанию = 0)

                    * ``int e2`` - электронный уровень второй в конечном состоянии (значение по умолчанию = 0)
            </doc>
        </fdoc>
        */

        double k_bf_exch(double T, kappa::Molecule const &molecule_before, kappa::Atom const &atom_before,
                         kappa::Molecule const &molecule_after, kappa::Atom const &atom_after,
                         kappa::Interaction const &interaction, int i_before, int i_after, int e_before=0, int e_after=0);
        /*
        <fdoc>
            <doc lang="rus">
                Отношение прямого и обратного коэффициентов скорости обменной реакции

                **Параметры**:
                    
                    * ``double T`` - температура                      
                    
                    * ``const kappa::Molecule &molecule_before`` - молекула, для которой производится расчёт, в начальном состоянии 

                    * ``const kappa::Molecule &molecule_after`` - молекула, для которой производится расчет, в конечном состоянии

                    * ``kappa::Atom const &atom_before`` - рассматриваемый атом в начальном состоянии

                    * ``kappa::Atom const &atom_after`` - рассматриваемый атом в конечном состоянии

                    * ``const kappa::Interaction &interaction`` - данные о взаимодействии сталкивающихся частиц
                    
                    * ``int i_before`` - колебательный уровень молекулы в начальном состоянии

                    * ``int i_after`` - колебательный уровень молекулы в конечном состоянии

                    * ``int e_before`` - электронный уровень молекулы в начальном состоянии (значение по умолчанию = 0)

                    * ``int e_after`` - электронный уровень молекулы в конечном состоянии (значение по умолчанию = 0)
            </doc>
        </fdoc>
        */

        double k_bf_diss(double T, kappa::Molecule const &molecule, kappa::Atom const &atom1, kappa::Atom const &atom2, int i, int e=0);
        /*
        <fdoc>
            <doc lang="rus">
                Отношение прямого и обратного коэффициентов скорости диссоциации

                **Параметры**:

                    * ``double T`` - температура

                    * ``const kappa::Molecule &molecule`` - молекула, для которой производится расчёт

                    * ``kappa::Atom const &atom1`` - рассматриваемый атом в начальном состоянии

                    * ``kappa::Atom const &atom2`` - рассматриваемый атом в конечном состоянии

                    * ``int i`` - колебательный уровень молекулы

                    * ``int e`` - электронный уровень молекулы в начальном состоянии (значение по умолчанию = 0)
            </doc>
        </fdoc>
        */

        double k_bf_diss(double T, kappa::Molecule const &molecule, kappa::Atom const &atom, int i, int e=0);
        /*
        <fdoc>
            <doc lang="rus">
                Отношение прямого и обратного коэффициентов скорости диссоциации

                **Параметры**:

                    * ``double T`` - температура

                    * ``const kappa::Molecule &molecule`` - молекула, для которой производится расчёт

                    * ``kappa::Atom const &atom`` - рассматриваемый атом

                    * ``int i`` - колебательный уровень молекулы

                    * ``int e`` - электронный уровень молекулы (значение по умолчанию = 0)
            </doc>
        </fdoc>
        */
        // </method_docs>

        // double k_bf_bimolecular(double T, kappa::Molecule const &molecule1_before, kappa::Molecule const &molecule2_before,
        //                         kappa::Molecule const &molecule1_after, kappa::Molecule const &molecule2_after,
        //                         kappa::Interaction const &interaction, int i1_before, int i1_after, int i2_before,
        //                         int i2_after, int e1_before, int e1_after, int e2_before, int e2_after);


    };
}

#endif /* kappa_approximation_hpp */
