#ifdef PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <cstdio>
#include <memory>

#include "NecessaryHeaders.h"
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
		.def("getNumElec", &Mwfn::getNumElec, pybind11::arg("spin") = 0)
		.def("getNumCenters", &Mwfn::getNumCenters)
		.def("getNumBasis", &Mwfn::getNumBasis)
		.def("getNumIndBasis", &Mwfn::getNumIndBasis)
		.def("getNumPrims", &Mwfn::getNumPrims)
		.def("getNumShells", &Mwfn::getNumShells)
		.def("getNumPrimShells", &Mwfn::getNumPrimShells)
		.def("getCoefficientMatrix", &Mwfn::getCoefficientMatrix, pybind11::arg("spin") = 0)
		.def("setCoefficientMatrix", &Mwfn::setCoefficientMatrix, pybind11::arg("matrix"), pybind11::arg("spin") = 0)
		.def("getEnergy", &Mwfn::getEnergy, pybind11::arg("spin") = 0)
		.def("setEnergy", &Mwfn::setEnergy, pybind11::arg("energies"), pybind11::arg("spin") = 0)
		.def("getOccupation", &Mwfn::getOccupation, pybind11::arg("spin") = 0)
		.def("setOccupation", &Mwfn::setOccupation, pybind11::arg("occupancies"), pybind11::arg("spin") = 0)
		.def("getFock", &Mwfn::getFock, pybind11::arg("spin") = 0)
		.def("getDensity", &Mwfn::getDensity, pybind11::arg("spin") = 0)
		.def("getEnergyDensity", &Mwfn::getEnergyDensity, pybind11::arg("spin") = 0)
		.def("Shell2Atom", &Mwfn::Shell2Atom)
		.def("Atom2Shell", &Mwfn::Atom2Shell)
		.def("Atom2ShellList", &Mwfn::Atom2ShellList)
		.def("Basis2Atom", &Mwfn::Basis2Atom)
		.def("Atom2Basis", &Mwfn::Atom2Basis)
		.def("Atom2BasisList", &Mwfn::Atom2BasisList)
		.def("Basis2Shell", &Mwfn::Basis2Shell)
		.def("Shell2Basis", &Mwfn::Shell2Basis)
		.def("Shell2BasisList", &Mwfn::Shell2BasisList)
		.def("getShell", &Mwfn::getShell)
		.def("getSpins", &Mwfn::getSpins)
		.def("Orthogonalize", &Mwfn::Orthogonalize)
		.def("Clone", &Mwfn::Clone)
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
