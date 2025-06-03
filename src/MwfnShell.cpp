#ifdef PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#endif
#include <vector>
#include <string>
#include <cstdio>

#include "Macro.h"
#include "MwfnShell.h"

int MwfnShell::getSize(){
	if ( this->Type >= 0 )
		return ( this->Type + 1 ) * ( this->Type + 2 ) / 2;
	else
		return -2 * this->Type + 1;
}

int MwfnShell::getNumPrims(){
	return this->Exponents.size();
}

void MwfnShell::Print(){
	std::printf("Type: %d\n", this->Type);
	std::printf("Exponents and Coefficients:\n");
	for ( int i = 0; i < this->getNumPrims(); i++ )
		std::printf("%f %f\n", this->Exponents[i], this->Coefficients[i]);
}

#ifdef PYTHON
void Init_MwfnShell(pybind11::module_& m){
	pybind11::class_<MwfnShell>(m, "MwfnShell")
		.def_readwrite("Type", &MwfnShell::Type)
		.def_readwrite("Exponents", &MwfnShell::Exponents)
		.def_readwrite("Coefficients", &MwfnShell::Coefficients)
		.def_readwrite("NormalizedCoefficients", &MwfnShell::NormalizedCoefficients)
		.def("getSize", &MwfnShell::getSize)
		.def("getNumPrims", &MwfnShell::getNumPrims)
		.def("Print", &MwfnShell::Print);
}
#endif
