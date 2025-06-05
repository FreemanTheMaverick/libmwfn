class MwfnShell{ public:
	int Type = -114; // Angular momentum. Positive of cartesian and negative for read spherical harmonics.
	std::vector<double> Exponents = {}; // Exponents of the primitive gaussians
	std::vector<double> Coefficients = {}; // Coefficients of the primitive gaussians. The unnormalized values given in the basis set file.
	std::vector<double> NormalizedCoefficients = {}; // For Chinium.
	int getSize(); // Number of components, 1 for s, 3 for p, 5 for spherical d and 6 for cartesian d.
	int getNumPrims(); // Number of primitive gaussians.
	void Print();
};
