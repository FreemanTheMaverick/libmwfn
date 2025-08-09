#include <Eigen/Dense>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <map>
#include <regex>
#include <numeric>

#include "NecessaryHeaders.h"
#include "Macro.h"
#include "MwfnShell.h"
#include "MwfnCenter.h"
#include "MwfnOrbital.h"
#include "Mwfn.h"

// By default spin = 0
#define __Check_Spin_Type_Shift__\
	if ( spin != 0 && spin != 1 && spin != 2 ) throw std::runtime_error("Invalid spin type!");\
	const int shift = ( this->Wfntype == 1 && spin == 2 ) ? this->getNumIndBasis() : 0;

double Mwfn::getNumElec(int spin){
	double nelec = 0;
	for ( int i = 0; i < this->getNumBasis(); i++ ){
		if ( spin == 0 ) nelec += this->Orbitals[i].Occ; // Total number of electrons if spin == 0.
		else{
			if ( this->Wfntype == 0 ) nelec += this->Orbitals[i].Occ / 2;
			else if ( this->Wfntype == 1 ){
				if ( spin == 1 ) nelec += this->Orbitals[i].Occ;
				else if ( spin == 2 ) nelec += this->Orbitals[i + this->getNumIndBasis()].Occ;
			}else if ( this->Wfntype == 2 ){
				if ( spin == 1 ) nelec += this->Orbitals[i].Occ > 0 ? 1 : 0;
				else if ( spin == 2) nelec += this->Orbitals[i].Occ > 1 ? 1 : 0;
			}
		}
	}
	return nelec;
}

double Mwfn::getCharge(){
	double nuclear_charge = 0;
	for ( MwfnCenter& center : this->Centers )
		nuclear_charge += center.Nuclear_charge;
	double nelec = this->getNumElec();
	return nuclear_charge - nelec;
}

int Mwfn::getNumCenters(){
	return this->Centers.size();
}

int Mwfn::getNumBasis(){
	int nbasis = 0;
	for ( MwfnCenter& center : this->Centers ) for ( MwfnShell& shell : center.Shells )
		nbasis += shell.getSize();
	return nbasis;
}

int Mwfn::getNumIndBasis(){
	return this->Orbitals.size() / ( this->Wfntype == 1 ? 2 : 1);
}

int Mwfn::getNumPrims(){
	int nprims = 0;
	for ( MwfnCenter& center : this->Centers ) for ( MwfnShell& shell : center.Shells ){
		int l = std::abs(shell.Type);
		nprims += ( l + 1 ) * ( l + 2 ) / 2 * shell.getNumPrims();
	}
	return nprims;
}

int Mwfn::getNumShells(){
	int nshells = 0;
	for ( MwfnCenter& center : this->Centers )
		nshells += center.Shells.size();
	return nshells;
}

int Mwfn::getNumPrimShells(){
	int nprimshells = 0;
	for ( MwfnCenter& center : this->Centers ) for ( MwfnShell& shell : center.Shells )
		nprimshells += shell.getNumPrims();
	return nprimshells;
}

Eigen::MatrixXd Mwfn::getCoefficientMatrix(int spin){
	Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(this->getNumBasis(), this->getNumIndBasis());
	for ( int jcol = 0; jcol < this->getNumIndBasis(); jcol++ ){
		if ( spin == 0 ) throw std::runtime_error("You must specify a spin index for \"Mwfn::getCoefficientMatrix(int spin)\".");
		else if ( spin == 1 || this->Wfntype == 0 || this->Wfntype == 2 ) matrix.col(jcol) = this->Orbitals[jcol].Coeff;
		else if ( spin == 2 ) matrix.col(jcol) = this->Orbitals[jcol + this->getNumIndBasis()].Coeff;
	}
	return matrix;
}

void Mwfn::setCoefficientMatrix(Eigen::MatrixXd matrix, int spin){
	for ( int jcol = 0; jcol < this->getNumIndBasis(); jcol++ ){
		if ( spin == 0 ) throw std::runtime_error("You must specify a spin index for \"Mwfn::setCoefficientMatrix(Eigen::MatrixXd matrix, int spin)\".");
		else if ( spin == 1 || this->Wfntype == 0 || this->Wfntype == 2 ) this->Orbitals[jcol].Coeff = matrix.col(jcol);
		else if ( spin == 2 ) this->Orbitals[jcol + this->getNumIndBasis()].Coeff = matrix.col(jcol);
	}
}

Eigen::VectorXd Mwfn::getEnergy(int spin){
	Eigen::VectorXd energies(this->getNumIndBasis());
	for ( int iorbital = 0; iorbital < this->getNumIndBasis(); iorbital++ ){
		if ( spin == 0 ) throw std::runtime_error("You must specify a spin index for \"Mwfn::getEnergy(int spin)\".");
		if ( spin == 1 || this->Wfntype == 0 || this->Wfntype == 2 ) energies(iorbital) = this->Orbitals[iorbital].Energy;
		if ( spin == 2 ) energies(iorbital) = this->Orbitals[iorbital + this->getNumIndBasis()].Energy;
	}
	return energies;
}

void Mwfn::setEnergy(Eigen::VectorXd energies, int spin){
	for ( int iorbital = 0; iorbital < this->getNumIndBasis(); iorbital++ ){
		if ( spin == 0 ) throw std::runtime_error("You must specify a spin index for \"Mwfn::setEnergy(Eigen::VectorXd energies, int spin)\".");
		if ( spin == 1 || this->Wfntype == 0 || this->Wfntype == 2 ) this->Orbitals[iorbital].Energy = energies(iorbital);
		if ( spin == 2 ) this->Orbitals[iorbital + this->getNumIndBasis()].Energy = energies(iorbital);
	}
}

Eigen::VectorXd Mwfn::getOccupation(int spin){
	Eigen::VectorXd occupancies(this->getNumIndBasis());
	for ( int iorbital = 0; iorbital < this->getNumIndBasis(); iorbital++ ){
		if ( spin == 0 ) occupancies(iorbital) = this->Orbitals[iorbital].Occ;
		else{
			if ( this->Wfntype == 0 ) occupancies(iorbital) = this->Orbitals[iorbital].Occ / 2;
			else if ( this->Wfntype == 1 ){
				if ( spin == 1 ) occupancies(iorbital) = this->Orbitals[iorbital].Occ;
				else if ( spin == 2 ) occupancies(iorbital) = this->Orbitals[iorbital + this->getNumIndBasis()].Occ;
			}else if ( this->Wfntype == 2 ){
				if ( spin == 1 ) occupancies(iorbital) = this->Orbitals[iorbital].Occ > 0 ? 1 : 0;
				else if ( spin == 2 ) occupancies(iorbital) = this->Orbitals[iorbital].Occ > 1 ? 1 : 0;
			}
		}
	}
	return occupancies;
}

void Mwfn::setOccupation(Eigen::VectorXd occupancies, int spin){
	for ( int iorbital = 0; iorbital < this->getNumIndBasis(); iorbital++ ){
		if ( spin == 0 ) this->Orbitals[iorbital].Occ = occupancies(iorbital);
		else{
			if ( this->Wfntype == 0 ) this->Orbitals[iorbital].Occ = occupancies(iorbital) * 2;
			else if ( this->Wfntype == 1 ){
				if ( spin == 1 ) this->Orbitals[iorbital].Occ = occupancies(iorbital);
				else if ( spin == 2 ) this->Orbitals[iorbital + this->getNumIndBasis()].Occ = occupancies(iorbital);
			}else if ( this->Wfntype == 2 ) throw std::runtime_error("You cannot set occupation numbers with spin != 0!");
		}
	}
}

Eigen::MatrixXd Mwfn::getFock(int spin){
	const Eigen::MatrixXd S = this->Overlap;
	if ( S.size() == 0 ) throw std::runtime_error("Overlap matrix is missing!");
	const Eigen::DiagonalMatrix<double,-1,-1> E = this->getEnergy(spin).asDiagonal();
	const Eigen::MatrixXd C = this->getCoefficientMatrix(spin);
	return S * C * E * C.transpose() * S;
}

Eigen::MatrixXd Mwfn::getDensity(int spin){
	const Eigen::DiagonalMatrix<double,-1,-1> N = this->getOccupation(spin).asDiagonal();
	const Eigen::MatrixXd C = this->getCoefficientMatrix(spin);
	return C * N * C.transpose();
}

Eigen::MatrixXd Mwfn::getEnergyDensity(int spin){
	const Eigen::DiagonalMatrix<double,-1,-1> N = this->getOccupation(spin).asDiagonal();
	const Eigen::DiagonalMatrix<double,-1,-1> E = this->getEnergy(spin).asDiagonal();
	const Eigen::MatrixXd C = this->getCoefficientMatrix(spin);
	return C * N * E * C.transpose();
}

std::vector<int> Mwfn::Shell2Atom(){
	std::vector<int> shell2atom = {};
	shell2atom.reserve(this->getNumShells());
	for ( int icenter = 0; icenter < this->getNumCenters(); icenter++ ){
		for ( int jshell = 0; jshell < this->Centers[icenter].getNumShells(); jshell++ ){
			shell2atom.push_back(icenter);
		}
	}
	return shell2atom;
}

std::vector<int> Mwfn::Atom2Shell(){
	std::vector<int> atom2shell = {};
	atom2shell.reserve(this->getNumCenters());
	int ishell = 0;
	for ( MwfnCenter& center : this->Centers ){
		atom2shell.push_back(ishell);
		ishell += center.getNumShells();
	}
	return atom2shell;
}

std::vector<std::vector<int>> Mwfn::Atom2ShellList(){
	std::vector<std::vector<int>> shell_indices_by_center = {};
	std::vector<int> atom2shell = this->Atom2Shell();
	for ( int iatom = 0; iatom < this->getNumCenters(); iatom++ ){
		int shell_head = atom2shell[iatom];
		int shell_length = this->Centers[iatom].getNumShells();
		std::vector<int> this_center(shell_length);
		std::iota(this_center.begin(), this_center.end(), shell_head);
		shell_indices_by_center.push_back(this_center);
	}
	return shell_indices_by_center;
}

std::vector<int> Mwfn::Basis2Atom(){
	std::vector<int> bf2atom = {}; bf2atom.reserve(this->getNumBasis());
	for ( int icenter = 0; icenter < this->getNumCenters(); icenter++ ){
		for ( int jbasis = 0; jbasis < this->Centers[icenter].getNumBasis(); jbasis++ ){
			bf2atom.push_back(icenter);
		}
	}
	return bf2atom;
}

std::vector<int> Mwfn::Atom2Basis(){
	std::vector<int> atom2bf = {}; atom2bf.reserve(this->getNumCenters());
	int ibasis = 0;
	for ( MwfnCenter& center : this->Centers ){
		atom2bf.push_back(ibasis);
		ibasis += center.getNumBasis();
	}
	return atom2bf;
}

std::vector<std::vector<int>> Mwfn::Atom2BasisList(){
	std::vector<std::vector<int>> basis_indices_by_center = {};
	std::vector<int> atom2bf = this->Atom2Basis();
	for ( int iatom = 0; iatom < this->getNumCenters(); iatom++ ){
		int basis_head = atom2bf[iatom];
		int basis_length = this->Centers[iatom].getNumBasis();
		std::vector<int> this_center(basis_length);
		std::iota(this_center.begin(), this_center.end(), basis_head);
		basis_indices_by_center.push_back(this_center);
	}
	return basis_indices_by_center;
}

std::vector<int> Mwfn::Basis2Shell(){
	std::vector<int> bf2shell = {}; bf2shell.reserve(this->getNumBasis());
	for ( int icenter = 0, jshell = 0; icenter < this->getNumCenters(); icenter++ ){
		for ( int _ = 0; _ < this->Centers[icenter].getNumShells(); _++, jshell++ ){
			for ( int kbasis = 0; kbasis < this->Centers[icenter].Shells[_].getSize(); kbasis++ ){
				bf2shell.push_back(jshell);
			}
		}
	}
	return bf2shell;
}

std::vector<int> Mwfn::Shell2Basis(){
	std::vector<int> shell2bf = {}; shell2bf.reserve(this->getNumShells());
	int ibasis = 0;
	for ( MwfnCenter& center : this->Centers ) for ( MwfnShell& shell : center.Shells ){
		shell2bf.push_back(ibasis);
		ibasis += shell.getSize();
	}
	return shell2bf;
}

MwfnShell& Mwfn::getShell(int ishell){
	if ( ishell < 0 ) throw std::runtime_error("The shell index must be >= 0!");
	if ( ishell >= this->getNumShells() ) throw std::runtime_error("The shell index exceeds the total number!");
	for ( MwfnCenter& center : this->Centers ) for ( MwfnShell& shell : center.Shells ){
		if ( ishell == 0 ) return shell;
		else ishell--;
	}
	throw std::runtime_error("You shouldn't be here!");
}

std::vector<std::vector<int>> Mwfn::Shell2BasisList(){
	std::vector<std::vector<int>> basis_indices_by_shell = {};
	std::vector<int> shell2bf = this->Shell2Basis();
	for ( int ishell = 0; ishell < this->getNumShells(); ishell++ ){
		int basis_head = shell2bf[ishell];
		int basis_length = this->getShell(ishell).getSize();
		std::vector<int> this_shell(basis_length);
		std::iota(this_shell.begin(), this_shell.end(), basis_head);
		basis_indices_by_shell.push_back(this_shell);
	}
	return basis_indices_by_shell;
}

std::vector<int> Mwfn::getSpins(){
	switch ( this->Wfntype ){
		case 0: return std::vector<int>{1};
		case 1:
		case 2: return std::vector<int>{1, 2};
		default: throw std::runtime_error("Invalid wavefunction type!");
	}
}

static Eigen::MatrixXd GramSchmidt(Eigen::MatrixXd C, Eigen::MatrixXd S){
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
	const Eigen::MatrixXd Ssqrt = es.operatorSqrt();
	const Eigen::MatrixXd Sinvsqrt = es.operatorInverseSqrt();
	const Eigen::MatrixXd Cprime = Ssqrt * C;
	Eigen::HouseholderQR<Eigen::MatrixXd> qr(Cprime);
	const Eigen::MatrixXd Cprime_new = qr.householderQ();
	const Eigen::MatrixXd C_new = Sinvsqrt * Cprime_new;
	return C_new;
}

static Eigen::MatrixXd Lowdin(Eigen::MatrixXd C, Eigen::MatrixXd S){
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C.transpose() * S * C);
	const Eigen::MatrixXd X = es.operatorInverseSqrt();
	const Eigen::MatrixXd C_new = C * X;
	return C_new;
}

void Mwfn::Orthogonalize(std::string scheme){
	const Eigen::MatrixXd S = this->Overlap;
	if ( S.size() == 0 ) throw std::runtime_error("Overlap matrix is missing!");
	std::transform(scheme.begin(), scheme.end(), scheme.begin(), ::toupper);
	for  ( int spin : this->getSpins() ){
		Eigen::MatrixXd C = this->getCoefficientMatrix(spin);
		if ( scheme == "GRAMSCHMIDT" || scheme == "GS" || scheme == "GRAM" ) C = GramSchmidt(C, S);
		else if ( scheme == "LOWDIN" || scheme == "SYMMETRIC" || scheme == "SYM" ) C = Lowdin(C, S);
		this->setCoefficientMatrix(C, spin);
	}
}

std::unique_ptr<Mwfn> Mwfn::Clone(){
	return std::make_unique<Mwfn>(*this);
}
