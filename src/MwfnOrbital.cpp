#ifdef PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <string>

#include "Macro.h"
#include "MwfnOrbital.h"

void Init_MwfnOrbital(pybind11::module_& m){
	pybind11::class_<MwfnOrbital>(m, "MwfnOrbital")
		.def_readwrite("Type", &MwfnOrbital::Type)
		.def_readwrite("Energy", &MwfnOrbital::Energy)
		.def_readwrite("Occ", &MwfnOrbital::Occ)
		.def_readwrite("Sym", &MwfnOrbital::Sym)
		.def_readwrite("Coeff", &MwfnOrbital::Coeff);
}
#endif
