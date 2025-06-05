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
	if ( spin == 0 ){ // Total number of electrons if spin == 0.
		for ( int i = 0; i < (int)this->Orbitals.size(); i++ )
			nelec += this->Orbitals[i].Occ;
	}else{
		__Check_Spin_Type_Shift__
		for ( int i = 0; i < this->getNumIndBasis(); i++ )
			nelec += this->Orbitals[i + shift].Occ;
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
	return this->Orbitals.size() / ( this->Wfntype == 0 ? 1 : 2);
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

EigenMatrix Mwfn::getCoefficientMatrix(int spin){
	__Check_Spin_Type_Shift__
	EigenMatrix matrix = EigenZero(this->getNumBasis(), this->getNumIndBasis());
	for ( int jcol = 0; jcol < this->getNumIndBasis(); jcol++ )
		matrix.col(jcol) = this->Orbitals[jcol + shift].Coeff;
	return matrix;
}

void Mwfn::setCoefficientMatrix(EigenMatrix matrix, int spin){
	__Check_Spin_Type_Shift__
	for ( int jcol = 0; jcol < this->getNumIndBasis(); jcol++ )
		this->Orbitals[jcol + shift].Coeff = matrix.col(jcol);
}

EigenVector Mwfn::getEnergy(int spin){
	__Check_Spin_Type_Shift__
	EigenVector energies(this->getNumIndBasis());
	for ( int iorbital = 0; iorbital < this->getNumIndBasis(); iorbital++ )
		energies(iorbital) = this->Orbitals[iorbital + shift].Energy;
	return energies;
}

void Mwfn::setEnergy(EigenVector energies, int spin){
	__Check_Spin_Type_Shift__
	for ( int iorbital = 0; iorbital < this->getNumIndBasis(); iorbital++ )
		this->Orbitals[iorbital + shift].Energy = energies(iorbital);
}

EigenVector Mwfn::getOccupation(int spin){
	__Check_Spin_Type_Shift__
	EigenVector occupancies(this->getNumIndBasis());
	for ( int iorbital = 0; iorbital < this->getNumIndBasis(); iorbital++ )
		occupancies(iorbital) = this->Orbitals[iorbital + shift].Occ;
	return occupancies;
}

void Mwfn::setOccupation(EigenVector occupancies, int spin){
	__Check_Spin_Type_Shift__
	for ( int iorbital = 0; iorbital < std::min((int)occupancies.size(), this->getNumIndBasis()); iorbital++ )
		this->Orbitals[iorbital + shift].Occ = occupancies(iorbital);
}

EigenMatrix Mwfn::getFock(int spin){
	const EigenMatrix S = this->Overlap;
	if ( S.size() == 0 ) throw std::runtime_error("Overlap matrix is missing!");
	const EigenDiagonal E = this->getEnergy(spin).asDiagonal();
	const EigenMatrix C = this->getCoefficientMatrix(spin);
	return S * C * E * C.transpose() * S;
}

EigenMatrix Mwfn::getDensity(int spin){
	const EigenDiagonal N = this->getOccupation(spin).asDiagonal();
	const EigenMatrix C = this->getCoefficientMatrix(spin);
	return C * N * C.transpose();
}

EigenMatrix Mwfn::getEnergyDensity(int spin){
	const EigenDiagonal N = this->getOccupation(spin).asDiagonal();
	const EigenDiagonal E = this->getEnergy(spin).asDiagonal();
	const EigenMatrix C = this->getCoefficientMatrix(spin);
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
		case 0: return std::vector<int>{0};
		case 1: return std::vector<int>{1, 2};
		default: throw std::runtime_error("Invalid wavefunction type!");
	}
}

static EigenMatrix GramSchmidt(EigenMatrix C, EigenMatrix S){
	Eigen::SelfAdjointEigenSolver<EigenMatrix> es(S);
	const EigenMatrix Ssqrt = es.operatorSqrt();
	const EigenMatrix Sinvsqrt = es.operatorInverseSqrt();
	const EigenMatrix Cprime = Ssqrt * C;
	Eigen::HouseholderQR<EigenMatrix> qr(Cprime);
	const EigenMatrix Cprime_new = qr.householderQ();
	const EigenMatrix C_new = Sinvsqrt * Cprime_new;
	return C_new;
}

static EigenMatrix Lowdin(EigenMatrix C, EigenMatrix S){
	Eigen::SelfAdjointEigenSolver<EigenMatrix> es(C.transpose() * S * C);
	const EigenMatrix X = es.operatorInverseSqrt();
	const EigenMatrix C_new = C * X;
	return C_new;
}

void Mwfn::Orthogonalize(std::string scheme){
	const EigenMatrix S = this->Overlap;
	if ( S.size() == 0 ) throw std::runtime_error("Overlap matrix is missing!");
	std::transform(scheme.begin(), scheme.end(), scheme.begin(), ::toupper);
	for  ( int spin : this->getSpins() ){
		EigenMatrix C = this->getCoefficientMatrix(spin);
		if ( scheme == "GRAMSCHMIDT" || scheme == "GS" || scheme == "GRAM" ) C = GramSchmidt(C, S);
		else if ( scheme == "LOWDIN" || scheme == "SYMMETRIC" || scheme == "SYM" ) C = Lowdin(C, S);
		this->setCoefficientMatrix(C, spin);
	}
}

std::unique_ptr<Mwfn> Mwfn::Clone(){
	return std::make_unique<Mwfn>(*this);
}
