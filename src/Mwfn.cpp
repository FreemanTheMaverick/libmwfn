#ifdef PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <cstdio>

#include "Macro.h"
#include "MwfnShell.h"
#include "MwfnCenter.h"
#include "MwfnOrbital.h"
#include "Mwfn.h"

void Init_Mwfn(pybind11::module_& m){
	pybind11::class_<Mwfn>(m ,"Mwfn")
		.def_readwrite("Wfntype", &Mwfn::Wfntype)
		.def_readwrite("E_tot", &Mwfn::E_tot)
		.def_readwrite("VT_ratio", &Mwfn::VT_ratio)
		.def_readwrite("Centers", &Mwfn::Centers)
		.def_readwrite("Orbitals", &Mwfn::Orbitals)
		.def_readwrite("Overlap", &Mwfn::Overlap)
		.def("getCharge", &Mwfn::getCharge)
		.def("getNumElec", &Mwfn::getNumElec)
		.def("getNumCenters", &Mwfn::getNumCenters)
		.def("getNumBasis", &Mwfn::getNumBasis)
		.def("getNumIndBasis", &Mwfn::getNumIndBasis)
		.def("getNumPrims", &Mwfn::getNumPrims)
		.def("getNumPrimShells", &Mwfn::getNumPrimShells)
		.def("getCoefficientMatrix", &Mwfn::getCoefficientMatrix)
		.def("setCoefficientMatrix", &Mwfn::setCoefficientMatrix)
		.def("getEnergy", &Mwfn::getEnergy)
		.def("setEnergy", &Mwfn::setEnergy)
		.def("getOccupation", &Mwfn::getOccupation)
		.def("setOccupation", &Mwfn::setOccupation)
		.def("getFock", &Mwfn::getFock)
		.def("getDensity", &Mwfn::getDensity)
		.def("getEnergyDensity", &Mwfn::getEnergyDensity)
		.def("Basis2Atom", &Mwfn::Basis2Atom)
		.def("Atom2Basis", &Mwfn::Atom2Basis)
		.def(pybind11::init<>())
		.def(pybind11::init<std::string>())
		.def("Export", &Mwfn::Export)
		.def("PrintCenters", &Mwfn::PrintCenters)
		.def("PrintOrbitals", &Mwfn::PrintOrbitals)
		.def("setBasis", &Mwfn::setBasis)
		.def("setCenters", &Mwfn::setCenters)
		.def("NuclearRepulsion", &Mwfn::NuclearRepulsion);
}
#endif
