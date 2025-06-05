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
	EigenMatrix Overlap; // Overlap matrix among (normalized) basis functions.

	double getCharge();
	double getNumElec(int spin = 0);

	int getNumCenters(); // Number of atoms.

	int getNumBasis(); // Number of basis functions.
	int getNumIndBasis(); // Number of indenpendent basis functions, equivalent to number of molecular orbitals.
	int getNumPrims();
	int getNumShells(); // Number of shells.
	int getNumPrimShells();

	EigenMatrix getCoefficientMatrix(int spin = 0);
	void setCoefficientMatrix(EigenMatrix matrix, int spin = 0);

	EigenVector getEnergy(int spin = 0);
	void setEnergy(EigenVector energies, int spin = 0);
	EigenVector getOccupation(int spin = 0);
	void setOccupation(EigenVector occupancies, int spin = 0);
	EigenMatrix getFock(int spin = 0); // Fock matrix.
	EigenMatrix getDensity(int spin = 0); // Density matrix.
	EigenMatrix getEnergyDensity(int spin = 0); // Energy-weighted density matrix.

	std::vector<int> Shell2Atom(); // The i^th element is the index of the atom which the i^th shell originates from.
	std::vector<int> Atom2Shell(); // The i^th element is the index of the first shell that orignates from the i^th atom.
	std::vector<std::vector<int>> Atom2ShellList(); // The i^th list consists of the indeces of the shells that originate from the i^th atom.
	std::vector<int> Basis2Atom(); // The i^th element is the index of the atom which the i^th basis function originates from.
	std::vector<int> Atom2Basis(); // The i^th element is the index of the first basis function that orignates from the i^th atom.
	std::vector<std::vector<int>> Atom2BasisList(); // The i^th list consists of the indeces of the basis functions that originate from the i^th atom.
	std::vector<int> Basis2Shell(); // The i^th element is the index of the shell which the i^th basis function belongs to.
	std::vector<int> Shell2Basis(); // The i^th element is the index of the first basis function that belongs to the i^th shell.
	std::vector<std::vector<int>> Shell2BasisList(); // The i^th list consists of the indeces of the basis functions that belong to the i^th shell.

	MwfnShell& getShell(int ishell); // The reference to the i^th shell.
	std::vector<int> getSpins(); // A list of spin indeces of the current wavefunction type. [0] for Wfntype = 0; [1, 2] for Wfntype = 1.
	void Orthogonalize(std::string scheme); // Orthogonalize the orbitals in the scheme of "GramSchmidt" or "Lowdin".

	std::unique_ptr<Mwfn> Clone();
	Mwfn() = default;
	Mwfn(std::string mwfn_filename);
	void Export(std::string mwfn_filename);
	void PrintCenters();
	void PrintOrbitals();
	void setBasis(std::string basis_filename); // Advanced.
	void setCenters(std::vector<std::vector<double>> atoms); // Advanced.
	std::tuple<double, EigenMatrix, EigenMatrix> NuclearRepulsion(); // Nuclear repulsion energy and its first two nuclear derivatives.
};
