/* 
 * File:   interactionTest.cpp
 * Author: aspera
 *
 * Created on 30 марта 2016 г., 14:57
 */

#include <cstdlib>
#include <iostream>
#include "kappa.hpp"

using namespace std;
using namespace kappa;

int main(int argc, char** argv) {
    cout << "Start Test for Interaction classes" << endl;


    Molecule N2("N2");
    Molecule CO("CO");
	Atom N("N");


    Interaction int_N2_N2(N2, N2);
    Interaction int_N2_CO(N2, CO);
    Interaction int_CO_N2(CO, N2);


    cout << "N2+N2 phi_zero Born-Mayer = " << int_N2_N2["phi_zero Born-Mayer"] << endl;
    cout << "N2+N2 Bzowski, sigma = " << int_N2_N2["Bzowski, sigma"] << endl;
    
    //
    cout << "N2+CO Bzowski, sigma = " << int_N2_CO["Bzowski, sigma"] << endl;
    cout << "N2+CO Bzowski, V_0* = " << int_N2_CO["Bzowski, V_0*"] << endl;

    cout << "CO+N2 Bzowski, sigma = " << int_N2_CO["Bzowski, sigma"] << endl;
    cout << "CO+N2 Bzowski, V_0* = " << int_N2_CO["Bzowski, V_0*"] << endl;

	
	cout << int_N2_N2["diss," + N2.name + ",Arrh_A,Scanlon"] << endl;

    try {
        double value = int_N2_CO["C (VSS (1 1))"];
    }
    catch (const DataNotFoundException &e) {
        std::cout << e.what() << endl;
    }

	try {
		Interaction int_wrong(CO, N); // currently no data available
	}
	catch (const DataNotFoundException &e) {
		std::cout << e.what() << endl;
	}

    try {
        Interaction int_N2_N2_badfilename(N2, N2, "badfilename.yaml");  // test that trying to read a non-existent database file throws the correct error
    }
    catch (const UnopenedFileException &e) {
        std::cout << e.what() << endl;
    }
    
    string a;
    cout << "Enter anything to quit: ";
    cin >> a;
    return 0;
}

