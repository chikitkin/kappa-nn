// cvibr_multit.cpp
#include <cstdlib>
#include <iostream>
#include "../../src/kappa.hpp"

using namespace std;
using namespace kappa;
using namespace arma; 

int main(){
	double T, T1;
	Molecule Molecule1_H("N2", false, true );
	Molecule Molecule1_AnH("N2", true, true);
	Molecule Molecule2_H("O2", false, true);
	Molecule Molecule2_AnH("O2", true, true);
	ApproximationMultiT approx{};
	cout << Molecule1_H.name << endl;
	
	std::cout << endl;
	std::cout << endl;
	

	T1 = 2000.0; T = 100.0;
	cout << "Calculation of the vibrational heat capacity at T1 ="<<T1<<" and different values of T [500, 3000]" << endl;
	std::cout << "-c_vibr_T(T, T1, AnHarmonic)" << endl;
	while ( T <= 3000.0) {
		std::cout   << -approx.c_vibr_T(T, T1, Molecule1_AnH) << endl;
		T += 100;
	}
	std::cout << endl;
	std::cout << endl;
	T1 = 2000.0; T = 100.0;
	std::cout << "c_vibr_T1(T, T1, AnHarmonic)" << endl;
	while (T <= 3000.0) {
		std::cout  << approx.c_vibr_T1(T, T1, Molecule1_AnH) << endl;
		T += 100;
	}
	std::cout << endl;
	std::cout << "***************************************************************" << endl;
	T1 = 2000.0; T = 100.0;
	std::cout << "c_vibr_T1(T, T1, Harmonic)" << endl;
	while (T <= 3000.0) {
		std::cout  << approx.c_vibr_T1(T, T1, Molecule1_H) << endl;
		T += 100;
	}

	std::cout << endl;
	std::cout << endl;
	

	std::cout << "********************************************************************************************************" << endl;
	std::cout << "********************************************************************************************************" << endl;
	std::cout << "My chart" << endl;
	std::cout << "Molecule name = " << Molecule1_AnH.name << endl;
	T = 1000; T1 = 1000;
	cout << "Calculation of the vibrational heat capacity at T1 =" << T1 << " and different values of T [1000, 20000], anharmonic spectrum." << endl;
	cout << "- c_vibr_T(Anharmonic, N2)  " << endl;
	while (T <= 20000) {
		std::cout  << - approx.c_vibr_T(T, T1, Molecule1_AnH) << endl;
		T += 500;
	}

	std::cout << endl;
	T = 1000; T1 = 1000;
	cout << "Calculation of the vibrational heat capacity at T1 =" << T1 << " and different values of T [1000, 20000], harmonic spectrum." << endl;
	cout << " c_vibr_T(Harmonic, N2)  " << endl;
	while (T <= 20000) {
		std::cout  << approx.c_vibr_T(T, T1, Molecule1_H) << endl;
		T += 500;
	}

	std::cout << endl;
	std::cout << "***************************************************************" << endl;
	T = 1000; T1 = 1000;
	cout << "Calculation of the vibrational heat capacity at T1 =" << T1 << " and different values of T [1000, 20000], anharmonic spectrum." << endl;
	cout << " c_vibr_T1(Anharmonic, N2)  " << endl;
	while (T <= 20000) {
		std::cout  << approx.c_vibr_T1(T, T1, Molecule1_AnH) << endl;
		T += 500;
	}

	std::cout << endl;
	T = 1000; T1 = 1000;
	cout << "Calculation of the vibrational heat capacity at T1 =" << T1 << " and different values of T [1000, 20000], harmonic spectrum." << endl;
	cout << " c_vibr_T1(Harmonic, N2)  " << endl;
	while (T <= 20000) {
		std::cout << approx.c_vibr_T1(T, T1, Molecule1_H) << endl;
		T += 500;
	}
	std::cout << endl;
	std::cout << endl;
	std::cout << "***************************************************************" << endl;
	std::cout << "***************************************************************" << endl;



	T = 5000; T1 = 1000;
	cout << "Calculation of the vibrational heat capacity at T =" << T << " and different values of T1 [1000, 20000], anharmonic spectrum." << endl;
	cout << "- c_vibr_T(Anharmonic, N2)  " << endl;
	while (T1 <= 20000) {
		std::cout << -approx.c_vibr_T(T, T1, Molecule1_AnH) << endl;
		T1 += 500;
	}
	std::cout << endl;
	T = 5000; T1 = 1000;
	cout << "Calculation of the vibrational heat capacity at T =" << T << " and different values of T1 [1000, 20000], harmonic spectrum." << endl;
	cout << " c_vibr_T(Harmonic, N2)  " << endl;
	while (T1 <= 20000) {
		std::cout  << approx.c_vibr_T(T, T1, Molecule1_H) << endl;
		T1 += 500;
	}
	std::cout << endl;
	std::cout << "***************************************************************" << endl;
	T = 5000; T1 = 1000;
	cout << "Calculation of the vibrational heat capacity at T =" << T << " and different values of T1 [1000, 20000], anharmonic spectrum." << endl;
	cout << " c_vibr_T1(Anharmonic, N2)  " << endl;
	while (T1 <= 20000) {
		std::cout  << approx.c_vibr_T1(T, T1, Molecule1_AnH) << endl;
		T1 += 500;
	}
	std::cout << endl;
	T = 5000; T1 = 1000;
	cout << "Calculation of the vibrational heat capacity at T =" << T << " and different values of T1 [1000, 20000], harmonic spectrum." << endl;
	cout << " c_vibr_T1(Harmonic, N2)  " << endl;
	while (T1 <= 20000) {
		std::cout  << approx.c_vibr_T1(T, T1, Molecule1_H) << endl;
		T1 += 500;
	}
	std::cout << endl;
	std::cout << endl;
	std::cout << "***************************************************************" << endl;
	std::cout << "***************************************************************" << endl;
	std::cout << endl;
	std::cout << endl;
	std::cout << "Molecule name = " << Molecule2_AnH.name << endl;
	T = 1000; T1 = 1000;
	cout << "Calculation of the vibrational heat capacity at T1 =" << T1 << " and different values of T [1000, 20000], anharmonic spectrum." << endl;
	cout << "- c_vibr_T(Anharmonic, O2)  " << endl;
	while (T <= 20000) {
		std::cout << -approx.c_vibr_T(T, T1, Molecule2_AnH) << endl;
		T += 500;
	}
	std::cout << endl;
	T = 1000; T1 = 1000;
	cout << "Calculation of the vibrational heat capacity at T1 =" << T1 << " and different values of T [1000, 20000], harmonic spectrum." << endl;
	cout << " c_vibr_T(Harmonic, O2)  " << endl;
	while (T <= 20000) {
		std::cout << approx.c_vibr_T(T, T1, Molecule2_H) << endl;
		T += 500;
	}

	std::cout << endl;
	std::cout << "***************************************************************" << endl;
	T = 1000; T1 = 1000;
	cout << "Calculation of the vibrational heat capacity at T1 =" << T1 << " and different values of T [1000, 20000], anharmonic spectrum." << endl;
	cout << " c_vibr_T1(Anharmonic, O2)  " << endl;
	while (T <= 20000) {
		std::cout << approx.c_vibr_T1(T, T1, Molecule2_AnH) << endl;
		T += 500;
	}
	std::cout << endl;
	T = 1000; T1 = 1000;
	cout << "Calculation of the vibrational heat capacity at T1 =" << T1 << " and different values of T [1000, 20000], harmonic spectrum." << endl;
	cout << " c_vibr_T1(Harmonic, O2)  " << endl;
	while (T <= 20000) {
		std::cout  << approx.c_vibr_T1(T, T1, Molecule2_H) << endl;
		T += 500;
	}
	std::cout << endl;
	std::cout << endl;
	std::cout << "***************************************************************" << endl;
	std::cout << "***************************************************************" << endl;
	T = 5000; T1 = 1000;
	cout << "Calculation of the vibrational heat capacity at T =" << T << " and different values of T1 [1000, 20000], anharmonic spectrum." << endl;
	cout << "- c_vibr_T(Anharmonic, O2)  " << endl;
	while (T1 <= 20000) {
		std::cout  <<- approx.c_vibr_T(T, T1, Molecule2_AnH) << endl;
		T1 += 500;
	}
	std::cout << endl;
	T = 5000; T1 = 1000;
	cout << "Calculation of the vibrational heat capacity at T =" << T << " and different values of T1 [1000, 20000], harmonic spectrum." << endl;
	cout << " c_vibr_T(Harmonic, O2)  " << endl;
	while (T1 <= 20000) {
		std::cout << "T1=" << T1 << " c_vibr_T=" << approx.c_vibr_T(T, T1, Molecule2_H) << endl;
		T1 += 500;
	}
	std::cout << endl;
	std::cout << "***************************************************************" << endl;
	T = 5000; T1 = 1000;
	cout << "Calculation of the vibrational heat capacity at T =" << T << " and different values of T1 [0, 20000], anharmonic spectrum." << endl;
	cout << " c_vibr_T1(Anharmonic, O2)  " << endl;
	while (T1 <= 20000) {
		std::cout << approx.c_vibr_T1(T, T1, Molecule2_AnH) << endl;
		T1 += 500;
	}
	std::cout << endl;

	T = 5000; T1 = 1000;
	cout << "Calculation of the vibrational heat capacity at T =" << T << " and different values of T1 [0, 20000], harmonic spectrum." << endl;
	cout << " c_vibr_T1(Harmonic, O2)  " << endl;
	while (T1 <= 20000) {
		std::cout  << approx.c_vibr_T1(T, T1, Molecule2_H) << endl;
		T1 += 500;
	}

	system("pause");
}