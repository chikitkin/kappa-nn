/* 
 * File:   omegaTest.cpp
 */

#include <iostream>
#include "kappa.hpp"
#include <fstream>

using namespace std;
using namespace kappa;

int main(int argc, char** argv) {
    std::cout << "Start Test for k_exch coefficients, loading particle data" << endl;

    Molecule Molecule1("N2");
    Atom Atom1("O");

    std::cout << "Loading interaction parameters" << endl;

    Interaction inter(Molecule1, Atom1);
    Approximation ApproximationTest{};

    ofstream outf;
	
	std::vector<double> T_vals = { 500., 1000., 2000., 5000., 10000., 15000., 20000., 25000., 30000., 40000. };
	std::vector<int> i_vals = { 0, 10, 15, 20, 25 };
    std::vector<std::string> model_names = {"Warnatz", "R-F", "Polak"};
    std::vector<kappa::models_k_exch> model_vals = {models_k_exch::model_k_exch_warnatz, models_k_exch::model_k_exch_rf, models_k_exch::model_k_exch_polak};
	double res;
    int j =0;
    for (auto exch_model: model_vals) {
        outf.open(Molecule1.name + "_" + Atom1.name + "_" + model_names[j] + ".txt");
		outf << "T;";
		for (auto i : i_vals) {
			outf << "i=" << i << ",energy=" << Molecule1.vibr_energy[0][i] / K_CONST_EV << ";";
		}
		outf << endl;
        for (auto T : T_vals) {
			outf << T << ";";
            for (auto i : i_vals) {
                res = ApproximationTest.k_exch(T, Molecule1, Atom1, inter, i, exch_model);
				outf << res << ";";
            }
            outf << endl;
        }
        outf.close();
        j++;
    }
    
    string a;
    cout << "Enter anything to quit: ";
    cin >> a;
    return 0;
}

