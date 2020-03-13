/* 
 * File:   omegaTest.cpp
 */

#include <iostream>
#include "kappa.hpp"
#include <fstream>

using namespace std;
using namespace kappa;

void calculate(Molecule m, Atom a, Approximation appr) {
    std::cout << "Loading interaction parameters for " << m.name << "+" << a.name << endl;

    Interaction inter(m, a);

    ofstream outf;
    
    std::cout << "Calculating k_diss" << endl;

    std::vector<double> T_vals = { 500., 1000., 2000., 5000., 10000., 15000., 20000., 25000., 30000., 40000. };
    std::vector<int> i_vals;

    if (m.name == "O2") {
        i_vals = { 0, 5, 10, 20 };
    }
    else {
        i_vals = { 0, 5, 10, 20, 30 };
    }
    std::vector<std::string> non_TM_model_names = {"RS, center-of-mass line energy + vibr energy", "RS, full energy + vibr energy",
                                                   "VSS, center-of-mass line energy + vibr energy", "VSS, full energy + vibr energy",
                                                   "Phys4Entry approximation", "ILT"};
    std::vector<kappa::models_k_diss> non_TM_model_vals = {models_k_diss::model_k_diss_rs_thresh_cmass_vibr, models_k_diss::model_k_diss_rs_thresh_vibr,
                                                           models_k_diss::model_k_diss_vss_thresh_cmass_vibr, models_k_diss::model_k_diss_vss_thresh_vibr,
                                                           models_k_diss::model_k_diss_phys4entry, models_k_diss::model_k_diss_ilt};
    // std::vector<std::string> TM_model_names = {"Treanor-Marrone, U=inf, Scanlon", "Treanor-Marrone, U=D/6K, Scanlon", "Treanor-Marrone, U=3T, Scanlon"};
    // std::vector<kappa::models_k_diss> TM_model_vals = {model_k_diss_tm_infty_arrh_scanlon, model_k_diss_tm_D6k_arrh_scanlon,  model_k_diss_tm_3T_arrh_scanlon};
    std::vector<std::string> TM_model_names = {"Treanor-Marrone, U=inf, Park", "Treanor-Marrone, U=D/6K, Park", "Treanor-Marrone, U=3T, Park"};
    std::vector<kappa::models_k_diss> TM_model_vals = {models_k_diss::model_k_diss_tm_infty_arrh_park, models_k_diss::model_k_diss_tm_D6k_arrh_park, models_k_diss::model_k_diss_tm_3T_arrh_park};

    double res, tmp;
    outf.open(m.name + "_" + a.name + ".txt");

    outf << "T;Z(U=inf);Z(U=D/6K);Z(U=3T);;";

    int counter;

    for (auto i : i_vals) {
        counter = 0;
        for (auto diss_model: TM_model_vals) {
            outf << TM_model_names[counter] << ",i=" << i << "(" << m.vibr_energy[0][i] / K_CONST_EV << " eV);"; 
            counter++;
        }
        outf << ";";
    }

    for (auto i : i_vals) {
        counter = 0;
        for (auto diss_model: non_TM_model_vals) {
            outf << non_TM_model_names[counter] << ",i=" << i << "(" << m.vibr_energy[0][i] / K_CONST_EV << " eV);"; 
            counter++;
        }
    }

    outf << endl;

    for (auto T : T_vals) {
        outf << T << ";";

        outf << m.num_vibr_levels[0] << ";" << appr.Z_vibr(-m.diss_energy[0] / (6 * K_CONST_K), m) << ";" << appr.Z_vibr(-3 * T, m) << ";;";

        for (auto i : i_vals) {
            for (auto diss_model: TM_model_vals) {
                outf << appr.k_diss(T, m, inter, i, diss_model) << ";"; 
                counter++;
            }
            outf << ";";
        }
        for (auto i : i_vals) {
            for (auto diss_model: non_TM_model_vals) {
                outf << appr.k_diss(T, m, inter, i, diss_model) << ";"; 
                counter++;
            }
        }
        outf << endl;
    }
}

int main(int argc, char** argv) {
    std::cout << "Start Test for k_diss coefficients, loading particle data" << endl;

    Molecule N2("N2", false);
    Molecule O2("O2", false);
    Atom N("N");
    Atom O("O");
    Approximation ApproximationTest{};

    calculate(N2, N, ApproximationTest);
    calculate(O2, O, ApproximationTest);

    string a;
    cout << "Enter anything to quit: ";
    cin >> a;
    return 0;
}

