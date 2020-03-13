//
//  approximation_onet.hpp
//

#ifndef kappa_approximation_onet_hpp
#define kappa_approximation_onet_hpp

#include "approximation.hpp"

namespace kappa {
    class ApproximationOneT : public Approximation {
        // one-temperature approximation class
    protected:
        double Z_int(double T, const arma::vec &electron_energy, const arma::Col<unsigned long> &statistical_weight, int n_electron_levels);
        // calculate the internal partition function of an atom given a temperature T, values of the energies of the electron levels,
        // the statistical weights of the electron levels, the amount of electron levels to take into account

		double Z_int(double T, const arma::vec &electron_energy, const arma::Col<unsigned long> &statistical_weight, int n_electron_levels,
			         const std::vector<arma::vec> &vibr_energy, const std::vector<int> &num_vibr_levels,
			         const std::vector<std::vector<arma::vec> > &rot_energy, const std::vector<std::vector<int> > &num_rot_levels, int rot_symmetry);
        // calculate the internal partition function given a temperature T, values of the energies of the electron levels,
        // the statistical weights of the electron levels, the amount of electron levels, the vibrational energies, the amount of
        // vibrational levels, the rotational energies, the amount of rotational levels and
        // whether the molecule is symmetric

        double avg_energy(double T, const arma::vec &electron_energy, const arma::Col<unsigned long> &statistical_weight, int n_electron_levels);
        // calculate the averaged internal energy of an atom given a temperature T, values of the energies of the electron levels,
        // the statistical weights of the electron levels, the amount of electron levels

        double avg_energy(double T, const arma::vec &electron_energy, const arma::Col<unsigned long> &statistical_weight, int n_electron_levels,
                          const std::vector<arma::vec> &vibr_energy, const std::vector<int> &num_vibr_levels,
                          const std::vector<std::vector<arma::vec> > &rot_energy, const std::vector<std::vector<int> > &num_rot_levels, int rot_symmetry);
        // calculate the averaged internal energy given a temperature T, values of the energies of the electron levels,
        // the statistical weights of the electron levels, the amount of electron levels, the vibrational energies, the amount of
        // vibrational levels, the rotational energies, the amount of rotational levels and
        // whether the molecule is symmetric

        double avg_energy_sq(double T, const arma::vec &electron_energy, const arma::Col<unsigned long> &statistical_weight, int n_electron_levels);
        // calculate the averaged squared internal energy of an atom given a temperature T, values of the energies of the electron levels,
        // the statistical weights of the electron levels, the amount of electron levels

        double avg_energy_sq(double T, const arma::vec &electron_energy, const arma::Col<unsigned long> &statistical_weight, int n_electron_levels,
                             const std::vector<arma::vec> &vibr_energy, const std::vector<int> &num_vibr_levels,
                             const std::vector<std::vector<arma::vec> > &rot_energy, const std::vector<std::vector<int> > &num_rot_levels, int rot_symmetry);
        // calculate the averaged squared internal energy given a temperature T, values of the energies of the electron levels,
        // the statistical weights of the electron levels, the amount of electron levels, the vibrational energies, the amount of
        // vibrational levels, the rotational energies, the amount of rotational levels and
        // whether the molecule is symmetric

        double c_int(double T, double mass, const arma::vec &electron_energy, const arma::Col<unsigned long> &statistical_weight, int n_electron_levels);
        // calculate the specific heat capacity of internal degrees of freedom of an atom given a temperature T, the mass of the molecule, values of the energies of the electron levels,
        // the statistical weights of the electron levels, the amount of electron levels

        double c_int(double T, double mass, const arma::vec &electron_energy, const arma::Col<unsigned long> &statistical_weight, int n_electron_levels,
                     const std::vector<arma::vec> &vibr_energy, const std::vector<int> &num_vibr_levels,
                     const std::vector<std::vector<arma::vec> > &rot_energy, const std::vector<std::vector<int> > &num_rot_levels, int rot_symmetry);
        // calculate the specific heat capacity of internal degrees of freedom given a temperature T, the mass of the molecule, values of the energies of the electron levels,
        // the statistical weights of the electron levels, the amount of electron levels, the vibrational energies, the amount of
        // vibrational levels, the rotational energies, the amount of rotational levels and
        // whether the molecule is symmetric

    public:
        ApproximationOneT();

		double Z_int(double T, const kappa::Atom &atom, double Delta_E=-1);
		// calculate the internal partition function for an atom for a temperature T

        double Z_int(double T, const kappa::Molecule &molecule, int num_electron_levels=1);
        // calculate the internal partition function for a molecule for a temperature T

		double avg_energy(double T, const kappa::Atom &atom, double Delta_E=-1);
		// calculate the averaged internal energy for an atom for a temperature T

        double avg_energy(double T, const kappa::Molecule &molecule, int num_electron_levels=1);
        // calculate the averaged internal energy for a molecule for a temperature T

		double avg_energy_sq(double T, const kappa::Atom &atom, double Delta_E=-1);
		// calculate the averaged squared internal energy for an atom for a temperature T

        double avg_energy_sq(double T, const kappa::Molecule &molecule, int num_electron_levels=1);
        // calculate the averaged squared internal energy for a molecule for a temperature T

		double c_int(double T, const kappa::Atom &atom, double Delta_E=-1);
		// calculate the specific heat capacity of internal degrees of freedom for an atom for a temperature T

        double c_int(double T, const kappa::Molecule &molecule, int num_electron_levels=1);
        // calculate the specific heat capacity of internal degrees of freedom for a molecule for a temperature T
    };
}

#endif /* kappa_approximation_sts_hpp */
