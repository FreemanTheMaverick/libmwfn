# libmwfn
> A flexible .mwfn format handler

## Overview
`mwfn` is a text file format for gaussian-based wavefunctions of molecules and solids described [here](https://doi.org/10.26434/chemrxiv-2021-lt04f-v6).
This library, *libmwfn*, is a flexible handler of the `mwfn` format.
With *libmwfn*, the users are able to manipulate gaussian-based wavefunctions on their demand.
*libmwfn* is written in C++ and provides Python interface by `pybind11`.

The motive of this library is my need for a tool to deal with wavefunction file in my other two projects, Chinium and Orbaplaw, written in C++ and Python, respectively.

## Installation
### Prerequisites
* A C++ compiler that supports C++20 standard
* GNU make for C++
* `setuptools` and `pybind11` for Python
### C++
* Obtain the source code via `git clone`ing the repository or `wget`ing the newest release.
* Specify the C++ compiler and the `Eigen3` directory in `makefile`.
* Build the library with `make -j[N]` and you will find the newly created directories `include/` and `lib/`.
* Tryout. (Reading in an existent file `test.mwfn` and exporting its information to another `test_new.mwfn`)
```
// test.cpp
#include <libmwfn.h>
int main(){
	Mwfn mwfn("test.mwfn");
	mwfn.Export("test_new.mwfn");
	return 0;
}
// end of test.cpp
$ gcc test.cpp -I$(LIBMWFN)/include/ -L$(LIBMWFN)/lib/ -lmwfn # or -l:libmwfn.a
```
### Python
* Install via pip.
```
pip install libmwfn
```
* Tryout. (Reading in an existent file `test.mwfn` and exporting its information to another `test_new.mwfn`)
```
$ python
>>> import libmwfn as lm
>>> mwfn = lm.Mwfn("test.mwfn")
>>> mwfn.Export("test_new.mwfn")
```

