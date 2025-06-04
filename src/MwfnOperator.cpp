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
#include <memory>

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
