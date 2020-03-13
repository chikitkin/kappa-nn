/* 
 * File:   omegaTest.cpp
 */

#include <iostream>
#include "kappa.hpp"
#include <fstream>

using namespace std;
using namespace kappa;


int main(int argc, char** argv) {
    std::cout << "Start Test for Omega integrals, loading particle data" << endl;
	Approximation ApproximationTest{};

	Molecule N2("N2");
	Molecule O2("O2");
	Molecule NO("NO");

	Atom N("N");
	Atom O("O");

    std::vector<Molecule> molecules = {N2, O2, NO};
    std::vector<Atom> atoms = {N, O};

	std::vector<models_omega> omega_integral_models = {models_omega::model_omega_vss, models_omega::model_omega_lennardjones, models_omega::model_omega_esa };
    
	std::map<models_omega, std::string> omega_model_names;

	omega_model_names.insert(std::pair<models_omega, std::string>(models_omega::model_omega_vss, "VSS"));
	omega_model_names.insert(std::pair<models_omega, std::string>(models_omega::model_omega_bornmayer, "Born-Mayer"));
	omega_model_names.insert(std::pair<models_omega, std::string>(models_omega::model_omega_lennardjones, "LJ"));
	omega_model_names.insert(std::pair<models_omega, std::string>(models_omega::model_omega_esa, "ESA"));

	int l = 1;
    int r = 1;
    std::vector<double> T_vals = {200., 300., 400., 500., 600., 700., 800., 900., 1000., 1200.};
    // std::vector<double> T_vals = {500., 1000., 2000., 5000., 10000., 15000., 20000., 25000., 30000., 40000.};
    std::vector<double> res_vals;
    std::string curr_model;
	ofstream out_m1_m2;

    for (int i=0; i<molecules.size(); i++) {
        for (int k=i; k<molecules.size(); k++) {
			std::cout << molecules[i].name + ", " + molecules[k].name << endl;
			Interaction int_m1_m2(molecules[i], molecules[k]);
            out_m1_m2.open(molecules[i].name + "_" + molecules[k].name + "_" + to_string(l) + "_" + to_string(r) + ".txt");
			out_m1_m2 << "T;collision diameter;";
			for (auto model : omega_integral_models) {
				curr_model = omega_model_names[model];
				out_m1_m2 << curr_model << "(" << to_string(l) + "," + to_string(r) + ")" << ";";
			}
			out_m1_m2 << endl;
			for (auto T : T_vals) {
				out_m1_m2 << T << ";" << int_m1_m2["collision diameter"] << "; ";
				for (auto model : omega_integral_models) {
					try {
						out_m1_m2 << ApproximationTest.omega_integral(T, int_m1_m2, l, r, model, true) << ";";
					}
					catch (const ModelParameterException &e) {
						std::cout << e.what() << endl;
						out_m1_m2 << ";";
					}
				}
				out_m1_m2 << endl;
			}
			out_m1_m2.close();
		}
    }

	for (int i = 0; i<molecules.size(); i++) {
		for (int k = 0; k<atoms.size(); k++) {
			Interaction int_m1_a1(molecules[i], atoms[k]);
			std::cout << molecules[i].name + ", " + atoms[k].name << endl;
			out_m1_m2.open(molecules[i].name + "_" + atoms[k].name + "_" + to_string(l) + "_" + to_string(r) + ".txt");
			out_m1_m2 << "T;collision diameter;";
			for (auto model : omega_integral_models) {
				curr_model = omega_model_names[model];
				out_m1_m2 << curr_model << "(" << to_string(l) + "," + to_string(r) + ")" << ";";
			}
			out_m1_m2 << endl;
			for (auto T : T_vals) {
				out_m1_m2 << T << ";" << int_m1_a1["collision diameter"] << ";";
				for (auto model : omega_integral_models) {
					try {
						out_m1_m2 << ApproximationTest.omega_integral(T, int_m1_a1, l, r, model, true) << ";";
					}
					catch (const ModelParameterException &e) {
						std::cout << e.what() << endl;
						out_m1_m2 << ";";
					}
				}
				out_m1_m2 << endl;
			}
			out_m1_m2.close();
		}
	}

	for (int i = 0; i<atoms.size(); i++) {
		for (int k = i; k<atoms.size(); k++) {
			Interaction int_a1_a2(atoms[i], atoms[k]);
			std::cout << atoms[i].name + ", " + atoms[k].name << endl;
			out_m1_m2.open(atoms[i].name + "_" + atoms[k].name + "_" + to_string(l) + "_" + to_string(r) + ".txt");
			out_m1_m2 << "T;collision diameter;";
			for (auto model : omega_integral_models) {
				curr_model = omega_model_names[model];
				out_m1_m2 << curr_model << "(" << to_string(l) + "," + to_string(r) + ")" << ";";
			}
			out_m1_m2 << endl;
			for (auto T : T_vals) {
				out_m1_m2 << T << ";" << int_a1_a2["collision diameter"] << ";";
				for (auto model : omega_integral_models) {
					try {
						out_m1_m2 << ApproximationTest.omega_integral(T, int_a1_a2, l, r, model, true) << ";";
					}
					catch (const ModelParameterException &e) {
						std::cout << e.what() << endl;
						out_m1_m2 << ";";
					}
				}
				out_m1_m2 << endl;
			}
			out_m1_m2.close();
		}
	}
    
    string a;
    cout << "Enter anything to quit: ";
    cin >> a;
    return 0;
}

