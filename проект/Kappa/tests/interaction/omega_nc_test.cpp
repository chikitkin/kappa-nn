#include <iostream>
#include "kappa.hpp"
#include <fstream>

using namespace std;
using namespace kappa;


void calculate(Interaction interaction, int r, int l) {
    Approximation ApproximationTest{};

    ofstream out_m1_m2;


    out_m1_m2.open(interaction.particle1_name + "_" 
        + interaction.particle2_name + "_" + to_string(l) + "_" + to_string(r) + ".txt");
    out_m1_m2 << "T;omega" << endl;

    vector<double> T_vals = {1000., 10000., 19000., 28000.};

    for (auto T : T_vals) {
        out_m1_m2 << T << "; ";
        out_m1_m2 << ApproximationTest.omega_integral(T, interaction, l, r, kappa::models_omega::model_omega_esa, false);
        out_m1_m2 << endl;
    }

    out_m1_m2.close();
}


int main(int argc, char** argv) {
    cout << "Start Test for Omega integrals, loading particle data..." << endl;
    
    Molecule mol_n2("N2");

    Atom at_n_ion("N+");
    Atom at_n("N");

    Atom at_o_ion("O+");
    Atom at_o("O");

    
    vector<double> res_vals;
    string curr_model;

    Interaction int_o_o_ion(at_o_ion, at_o);
    Interaction int_o_o(at_o, at_o);
    Interaction int_n_n_ion(at_n_ion, at_n);
    Interaction int_n_n(at_n, at_n);

    int l = 1;
    int r = 1;

    calculate(int_o_o_ion, r, l);
    calculate(int_o_o, r, l);
    calculate(int_n_n_ion, r, l);
    calculate(int_n_n, r, l);

    l = 2;
    r = 2;

    calculate(int_o_o_ion, r, l);
    calculate(int_o_o, r, l);
    calculate(int_n_n_ion, r, l);
    calculate(int_n_n, r, l);

    string a;
    cout << "Press enter to quit...";
    cin.ignore();

    return 0;
}