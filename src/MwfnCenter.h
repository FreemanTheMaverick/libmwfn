class MwfnCenter{ public:
	int Index = -114;
	double Nuclear_charge = -114;
	std::vector<double> Coordinates = {114, 514, 1919810}; // a.u.
	std::vector<MwfnShell> Shells = {};
	MwfnCenter() = default;
	int getNumShells(); // Number of shells originating from the atom.
	int getNumBasis(); // Number of basis functions originating from the atom.
	std::string getSymbol(); // The atomic symbol.
	void Print();
};
