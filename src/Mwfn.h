class Mwfn{ public:
	// Field 1
	int Wfntype = -114;
	double E_tot = -114;
	double VT_ratio = -114;

	// Fields 2 & 3
	std::vector<MwfnCenter> Centers = {};

	// Field 4
	std::vector<MwfnOrbital> Orbitals = {};

	// Field 5
	EigenMatrix Overlap;

	double getCharge();
	double getNumElec(int spin = 0);

	int getNumCenters();

	int getNumBasis();
	int getNumIndBasis();
	int getNumPrims();
	int getNumShells();
	int getNumPrimShells();

	EigenMatrix getCoefficientMatrix(int spin = 0);
	void setCoefficientMatrix(EigenMatrix matrix, int spin = 0);

	EigenVector getEnergy(int spin = 0);
	void setEnergy(EigenVector energies, int spin = 0);
	EigenVector getOccupation(int spin = 0);
	void setOccupation(EigenVector occupancies, int spin = 0);
	EigenMatrix getFock(int spin = 0);
	EigenMatrix getDensity(int spin = 0);
	EigenMatrix getEnergyDensity(int spin = 0);

	std::vector<int> Shell2Atom();
	std::vector<int> Atom2Shell();
	std::vector<std::vector<int>> Atom2ShellList();
	std::vector<int> Basis2Atom();
	std::vector<int> Atom2Basis();
	std::vector<std::vector<int>> Atom2BasisList();
	std::vector<int> Basis2Shell();
	std::vector<int> Shell2Basis();
	std::vector<std::vector<int>> Shell2BasisList();
	
	MwfnShell& getShell(int ishell);
	std::vector<int> getSpins();
	void Orthogonalize(std::string scheme);

	std::unique_ptr<Mwfn> Clone();
	Mwfn() = default;
	Mwfn(std::string mwfn_filename);
	void Export(std::string mwfn_filename);
	void PrintCenters();
	void PrintOrbitals();
	void setBasis(std::string basis_filename);
	void setCenters(std::vector<std::vector<double>> atoms);
	std::tuple<double, EigenMatrix, EigenMatrix> NuclearRepulsion();
};
