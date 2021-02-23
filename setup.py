from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(
    name='cyadd',
    ext_modules=cythonize("cyadd.pyx"),
    include_dirs=[numpy.get_include()],
    zip_safe=False,
)
