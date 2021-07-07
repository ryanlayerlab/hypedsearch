from setuptools import Extension, setup
from Cython.Build import cythonize
import os 

os.environ["CC"] = "g++"

cpp_source_dir = './cppModules'

extensions = [
    Extension(
            name="cppModules.gen_spectra", 
            sources=[f"{cpp_source_dir}/spectraGeneration/gen_spectra.pyx", f"{cpp_source_dir}/spectraGeneration/genSpectra.cpp"],
            language="c++", 
            extra_compile_args=["-std=c++11"]
        )
]

setup(
    name="cppModules",
    ext_modules=cythonize(
        extensions
    )
)
