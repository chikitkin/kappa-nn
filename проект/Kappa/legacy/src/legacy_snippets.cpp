

double kappa::ApproximationSTS::c_rot(double mass) {
	// simplified formula for specific heat capacity of rotational degrees of freedom for main electron state
	return K_CONST_K / mass;
}

double kappa::ApproximationSTS::Z_rot(double T, double rot_inertia, int rot_symmetry) {
	// simplified formula for rotational partition function for main electron state
    return 8 * K_CONST_PI * K_CONST_PI * rot_inertia * K_CONST_K * T / (rot_symmetry * K_CONST_H * K_CONST_H);
}