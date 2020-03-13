//
//  numeric.hpp
//  Various numeric methods and helper functions go here
//

#ifndef kappa_numeric_hpp
#define kappa_numeric_hpp

#include <functional>
#include <cmath>
#include <armadillo>

namespace kappa {
    double integrate_semi_inf(std::function<double(double) > f, double a = 0, int subdivisions = 5, double* error_estimate = nullptr);
    // integrate a function f on the semi-infinite interval [a,+infinity)
    // it performs a change of variables and then integrates over the interval [0,1]
    // the interval is split into equal-sized subintervals and the integral is computed over each
    // subinterval, the amount of subintervals is determined by the subdivisions parameter
    // if error_estimate is not a NULL pointer, write the error_estimate to the passed variable
    //
    // default values: int subdivisions=5, double a = 0.0, double* error_estimate = nullptr
    //
    // usage: to integrate a function f(double x1, double x2, int x3)
    // over x1 from 0 to infinity (x2, x3 are fixed values defined somewhere in the code):
    //
    // auto integrand = [x2, x3](double x1) {return f(x1, x2, x3)};
    // double result = kappa::integrate_semi_inf(integrand);

    double integrate_interval(std::function<double (double) > f, double a, double b, double* error_estimate = nullptr);
    // integrate a function f on the finite interval [a, b], where b > a
    // if error_estimate is not a NULL pointer, write the error_estimate to the passed variable
    //
    // default values: double* error_estimate = nullptr
    //
    // usage: to integrate a function f(double x1, double x2, int x3)
    // over x1 from -1 to 1 (x2, x3 are fixed values defined somewhere in the code):
    //
    // auto integrand = [x2, x3](double x1) {return f(x1, x2, x3)};
    // double result = kappa::integrate_semi_interval(integrand, -1, 1);

	int find_max_value(std::function<double(int) > f, double max_value, int start=0);
	// find the value of a parameter i of type int such that f(i) < max_value, f(i+1) > max_value
	// the search starts from i=start and then i is increased by 1 on each iteration (up to INT_MAX-2)
	// returns -1 if no value is found up to INT_MAX-2
	//
	// default values: int start=0

    double factorial(int n);
    // Calculate n!

    double fact_div_fact(int start, int end);
    // Calculate end! / start!, where start and end are integers

    double convert_cm_to_Joule(double x);
    // Convert a value given in cm^-1 to Joules

    const extern arma::vec::fixed<70> factorial_table; // precomputed values of factorials from 0! to 69!, factorial_table[i] = i!

    arma::vec::fixed<70> compute_factorial_table(void);

    constexpr double Born_Mayer_coeff_array[6][4] = {-267.0, 201.570, 174.672, 54.305,
                          26700, -19226.5, -27693.8, -10860.9,
                          -8.9e5, 6.3201e5, 1.0227e6, 5.4304e5,
                          -33.0838, 20.0862, 72.1059, 68.5001,
                          101.571, -56.4472, -286.393, -315.4531,
                          -87.7036, 46.3130, 277.146, 363.1807 }; // coefficients for Born-Mayer potential, Kustova-Nagnibeda (5.100)
}

#endif /* kappa_numeric_hpp */
