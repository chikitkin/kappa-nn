#include <iostream>
#include "kappa.hpp"
#include <ctime>
#include <chrono>

using namespace std;
using namespace kappa;


int main(int argc, char** argv) {
    std::cout << "Start Test for Omega integrals, loading particle data" << endl;
    Approximation ApproximationTest{};

    double n = 1.25 / (2 * 23.2651287e-27); //particle number density of N2

    Molecule N2("N2");
    Interaction int_m1_m2(N2, N2);
    double res1=0., res2=0;
    // std::vector<double> T_vals = {200., 300., 400., 500., 600., 700., 800., 900., 1000.};
	std::string curr_model;
    ofstream out_m1_m2;
    
    clock_t t;
    // t = clock();
    auto begin = std::chrono::high_resolution_clock::now();
    for (int i=0; i<20000; i++) {
    	res1 += ApproximationTest.Z_coll(2000. + i, 1e23, int_m1_m2);
    }
    // t = clock()-t;
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
    cout << duration << "ns total, average : " << duration / 20000 << "ns." << std::endl;
    // cout << "Res: " << res1 << "clicks: " << t << ", time: " << float(t) / CLOCKS_PER_SEC << " seconds" << std::endl;

    // t = clock();
    begin = std::chrono::high_resolution_clock::now();
    for (int i=0; i<20000; i++) {
        res2 += ApproximationTest.Z_coll(2000. + i, 1e23, int_m1_m2);
    }
	end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
    cout << duration << "ns total, average : " << duration / 20000 << "ns." << std::endl;
    // t = clock()-t;
    
    //cout << "New optimised Res: " << res2 << "clicks: " << t << ", time: " << float(t) / CLOCKS_PER_SEC << " seconds" << std::endl;
    cout << "Z=" <<  ApproximationTest.Z_coll(2500, n, int_m1_m2) << ", Z_opt=" <<  ApproximationTest.Z_coll2(2500, n, int_m1_m2) << std::endl;

    string a;
    cout << "Enter anything to quit: ";
    cin >> a;
    return 0;
}

