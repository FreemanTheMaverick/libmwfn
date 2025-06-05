# libmwfn
> A flexible .mwfn format handler

## Overview
`mwfn` is a text file format for gaussian-based wavefunctions of molecules and solids (described [here](https://doi.org/10.26434/chemrxiv-2021-lt04f-v6)).
This library, **libmwfn**, is a flexible handler of the `mwfn` format.
With **libmwfn**, the users are able to manipulate gaussian-based wavefunctions on their demand.
**libmwfn** is written in C++ and provides Python interface by `pybind11`.

The motive of this library is my need for a tool to deal with wavefunction files in my other two projects, [Chinium](https://github.com/FreemanTheMaverick/Chinium.git) and [Orbaplaw](https://github.com/FreemanTheMaverick/Orbaplaw.git), written in C++ and Python, respectively.

## Installation
### Prerequisites
* A C++ compiler that supports C++20 standard
* GNU `make` and `Eigen3` for C++
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
$ pip install libmwfn
```
* Tryout. (Reading in an existent file `test.mwfn` and exporting its information to another `test_new.mwfn`)
```
$ python
>>> import libmwfn as lm
>>> mwfn = lm.Mwfn("test.mwfn")
>>> mwfn.Export("test_new.mwfn")
```

## Usage
Here is a typical procedure to deal with wavefunctions in `Orbaplaw`.

### *Ab initio* quantum chemistry computation
A wavefunction file given by an *ab initio* calculation is needed.
This can be done by a variety of computational chemistry packages.
Here we use `Gaussian 16` as an example.
```
! Gaussian input file
%nprocshared=40
%mem=60GB
%chk=job.chk
# b3lyp 6-31g(d) 5d 7f

Title Card Required

0 1
[Geometry]
```

### Converting wavefunction file
The resultant wavefunction is stored in `job.chk`, an `chk` format file which is not supported by **libmwfn**.
We need to transform `job.chk` to `fchk` with `formchk` and then to `mwfn` format with [Multiwfn](http://sobereva.com/multiwfn).

```
$ # Shell
$ formchk job.chk # Now we have job.fchk
$ Multiwfn job.fchk
100 # Other functions (Part 1)
2 # Export various files (mwfn/pdb/xyz/wfn/wfx/molden/fch/47/mkl...) or generate input file of quantum chemistry programs
32 # Output current wavefunction as .mwfn file
# Default name: job.mwfn
2 # Export wavefunction, density matrix and overlap matrix
0 # Return
q # Exit Multiwfn.
```
Now we have `job.mwfn`.

### Wavefunction processing
+ Loading the `libmwfn` module.
```c++
// C++
#include <libmwfn.h>
```
```python
# Python
import libmwfn as lm
```

+ Reading and exporting wavefunction information from `job.mwfn` to `job_new.mwfn`.
```c++
// C++
Mwfn job_mwfn("job.mwfn");
job_mwfn.Export("job_new.mwfn");
```
```python
# Python
job_mwfn = lm.Mwfn("job.mwfn")
job_mwfn.Export("job_new.mwfn")
```

+ Getting the numbers of electrons.
```c++
// C++
int spin = 0; // 0 - default spin (depending on wavefunction types and functions); 1 - alpha; 2 - beta.
double n = job_mwfn.getNumElec(spin); // spin is optional. Default is 0, the total number of electrons of two spin types.
double n_alpha = mwfn.getNumElec(1); // The number of alpha electrons.
double n_beta = mwfn.getNumElec(2); // The number of beta electrons.
double charge = mwfn.getCharge(); // Total nuclear charges minus the number of electrons.
```
```python
# Python
n = job_mwfn.getNumElec() # or mwfn.getNumElec(0)
n_alpha = job_mwfn.getNumElec(1)
n_beta = job_mwfn.getNumElec(2)
charge = job_mwfn.getCharge()
```
Note that the returned values are all floating-point numbers for the general case of fractional occupation.

+ Getting and setting the occupation numbers.
```c++
// C++
int homo = std::round(job_mwfn.getNumElec(0) / 2) - 1; // The indices of the HOMO and LUMO.
int lumo = std::round(job_mwfn.getNumElec(0) / 2);
Eigen::VectorXd N = job_mwfn.getOccupation(0); // 0 for spin-restricted, 1 and 2 for alpha and beta in spin-unrestricted.
N(homo) = 0;
N(lumo) = 2; // Swapping the occupation numbers of the HOMO and LUMO.
mwfn.setOccupation(N, 0);
```

+ Getting and setting the orbital energies.
```c++
// C++
Eigen::VectorXd E = job_mwfn.getEnergy(0);
E.setZero();
job_mwfn.setEnergy(E, 0); // Setting the orbital energies to zeros.
```

+ Getting and setting the coefficient matrix.
```c++
// C++
Eigen::MatrixXd C = job_mwfn.getCoefficientMatrix(0);
Eigen::VectorXd C_0 = C.col(0); // The coefficient vector of the first orbital.
Eigen::VectorXd C_1 = C.col(1); // The coefficient vector of the second orbital.
C.col(0) = (C_0 + C_1) / 1.414213562; // Mixing the first two orbitals.
C.col(1) = (C_0 - C_1) / 1.414213562;
job_mwfn.setCoefficientMatrix(C, 0);
```

+ Make a (deep) copy of a `Mwfn` instance.
```c++
// C++
std::unique_ptr<Mwfn> another_mwfn = job_mwfn.Clone();
std::cout << another_mwfn->getNumElec() << std::endl; // The copy is a pointer.
```
```python
# Python
another_mwfn = job_mwfn.Clone()
print(another_mwfn.getNumElec()) # No need to worry about the pointer thing.
```

+ Check the header `libmwfn.h` for more functions.

## Caution
+ Ordering of basis functions.
The ordering of basis functions of $l \ge 2$ in **libmwfn** is [d-2, d-1, d0, d+1, d+2] and [f-3, f-2, f-1, f0, f+1, f+2, f+3], etc., different from [d0, d+1, d-1, d+2, d-2] and [f0, f+1, f-1, f+2, f-2, f+3, f-3, f+4, f-4], etc., in the `mwfn` file format.
Upon I/O of `mwfn` file, a matrix transformation is applied to accommodate this difference.
