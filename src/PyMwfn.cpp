#ifdef PYTHON
#include <pybind11/pybind11.h>

void Init_MwfnShell(pybind11::module_& m);
void Init_MwfnCenter(pybind11::module_& m);
void Init_MwfnOrbital(pybind11::module_& m);
void Init_Mwfn(pybind11::module_& m);
PYBIND11_MODULE(libmwfn, m){
	Init_MwfnShell(m);
	Init_MwfnCenter(m);
	Init_MwfnOrbital(m);
	Init_Mwfn(m);
}
#endif
