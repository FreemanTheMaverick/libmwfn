# libmwfn
> A flexible .mwfn format handler

## Overview

## Installation
### Prerequisites
* A C++ compiler that supports C++20 standard
* GNU make for C++
* `setuptools` and `pybind11` for Python
### C++
* Obtain the source code via `git clone`ing the repository or `wget`ing the newest release.
* Specify the C++ compiler and the `Eigen3` directory in `makefile`.
* Build the library with `make -j[N]` and you will find the newly created directories `include/` and `lib/`.
* Tryout.
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
* Tryout.
```
$ python
>>> import libmwfn as lm
>>> mwfn = lm.Mwfn("test.mwfn")
>>> mwfn.Export("test_new.mwfn")
```

