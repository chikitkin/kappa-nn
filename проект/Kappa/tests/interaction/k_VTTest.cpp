/* 
 * File:   omegaTest.cpp
 */

#include <iostream>
#include "kappa.hpp"
#include <fstream>

using namespace std;
using namespace kappa;

void calculate(Molecule m, Particle a, Approximation appr) {
    std::cout << "Loading interaction parameters for " << m.name << "+" << a.name << endl;

    Interaction inter(m, a, "c:/Users/st024385/Documents/DB/Interaction/interaction.yaml");

    ofstream outf;
    
    std::cout << "Calculating k_VT" << endl;

    std::vector<double> T_vals = { 500., 1000., 2000., 5000., 10000., 15000., 20000., 25000., 30000., 40000. };
    std::vector<int> i_vals;
    i_vals = { 1, 5, 10, 20, 30 };

    std::vector<std::string> model_names = {"Billing", "FHO+RS", "FHO+VSS"};
	std::vector<std::string> cs_model_names = {"FHO+RS"};
    std::vector<kappa::models_k_vt> VT_models = {models_k_vt::model_k_vt_billing, models_k_vt::model_k_vt_rs_fho, models_k_vt::model_k_vt_vss_fho};
    std::vector<kappa::models_cs_vt> VT_cs_models = {models_cs_vt::model_cs_vt_rs_fho};
    // std::vector<std::string> TM_model_names = {"Treanor-Marrone, U=inf, Scanlon", "Treanor-Marrone, U=D/6K, Scanlon", "Treanor-Marrone, U=3T, Scanlon"};
    // std::vector<kappa::models_k_diss> TM_model_vals = {model_k_diss_tm_infty_arrh_scanlon, model_k_diss_tm_D6k_arrh_scanlon,  model_k_diss_tm_3T_arrh_scanlon};

    double res, tmp;
    outf.open("c:/Users/st024385/Documents/Code/kappa-tests/VT-FHO/" + m.name + "_" + a.name + ".txt");

    outf << "T;";

    int counter;

    for (auto i : i_vals) {
        counter = 0;
        for (auto VT_model: VT_models) {
            outf << model_names[counter] << ",i=" << i << "(" << m.vibr_energy[0][i] / K_CONST_EV << " eV);"; 
            counter++;
        }
        outf << ";";
    }
    outf << endl;

    for (auto T : T_vals) {
        outf << T << ";";

        for (auto i : i_vals) {
            for (auto VT_model: VT_models) {
                outf << appr.k_VT(T, m, inter, i, -1, VT_model) << ";"; 
                counter++;
            }
            outf << ";";
        }
        outf << endl;
    }

    outf.close();
    double g;
    outf.open("c:/Users/st024385/Documents/Code/kappa-tests/VT-FHO/" + m.name + "_" + a.name + "_CS.txt");
    outf << "t;";
    for (auto i : i_vals) {
        counter = 0;
        for (auto VT_cs_model: VT_cs_models) {
            outf << cs_model_names[counter] << ",i=" << i << "(" << m.vibr_energy[0][i] / K_CONST_EV << " eV);";
            counter++;
        }
    }
    outf << endl;

    for (int k=0; k<40; k++) {
        g = k * 5000; // energy of relative motion in Kelvins
        outf << g << ";";
        g = sqrt(2 * g * K_CONST_K / inter.collision_mass);
        for (auto i : i_vals) {
            for (auto VT_cs_model: VT_cs_models) {
                outf << appr.crosssection_VT(g, m, inter, i, -1, VT_cs_model) << ";"; 
                counter++;
            }
        }
        outf << endl;
    }
    outf.close();
}

int main(int argc, char** argv) {
    std::cout << "Start Test for k_VT coefficients, loading particle data" << endl;

    Molecule N2("N2", false, true, "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
    // Molecule O2("O2", false);
    Atom N("N", "c:/Users/st024385/Documents/DB/Particles/particles.yaml");
    // Atom O("O");
    Approximation ApproximationTest{};

    calculate(N2, N, ApproximationTest);
    calculate(N2, N2, ApproximationTest);
    // calculate(O2, O, ApproximationTest);

    string a;
    cout << "Enter anything to quit: ";
    cin >> a;
    return 0;
}

