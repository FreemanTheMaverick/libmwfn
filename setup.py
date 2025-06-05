import os
import urllib.request
import tarfile
from glob import glob
from setuptools import setup
from setuptools.command.build import build
import pybind11
from pybind11.setup_helpers import Pybind11Extension, ParallelCompile, naive_recompile

__version__ = "0.0.2"

# Downloading Eigen3
pwd = os.path.dirname(__file__)
EIGEN3 = pwd + "/eigen-3.4.0/"
class CustomBuild(build):
	def run(self):
		url = "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz"
		dest = pwd + "/eigen-3.4.0.tar.gz"
		print("Downloading Eigen3 from %s to %s ..." % (url, dest))
		urllib.request.urlretrieve(url, dest)
		print("Extracting %s to %s ..." % (dest, EIGEN3))
		with tarfile.open(dest) as tar:
			tar.extractall(path = pwd) # Directory: eigen-3.4.0
		super().run()

ParallelCompile(
	"NPY_NUM_BUILD_JOBS",
	needs_recompile = naive_recompile
).install()

CPP = sorted(glob("src/*.cpp"))
HEADER = sorted(glob("src/*.h"))
ext_module = Pybind11Extension(
	"libmwfn",
	CPP,
	include_dirs = [EIGEN3],
	depends = HEADER,
	extra_compile_args = ["-O3", "-DPYTHON", "-DEIGEN_INITIALIZE_MATRICES_BY_ZERO"],
	language = "c++",
	cxx_std = 20
)

setup(
		name = "libmwfn",
		version = __version__,
		author = "FreemanTheMaverick",
		description = "A flexible .mwfn format handler",
		long_description = open("README.md").read(),
		long_description_content_type = "text/markdown",
		url = "https://github.com/FreemanTheMaverick/libmwfn.git",
		cmdclass = {"build": CustomBuild},
		ext_modules = [ext_module],
		classifiers = ["Programming Language :: Python :: 3"]
)
