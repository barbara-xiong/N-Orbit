from distutils.core import setup
from Cython.Build import cythonize
import numpy
import pandas as pd
from collections import Counter

setup(
    ext_modules=cythonize("permutation_test.pyx", annotate=True, compiler_directives={'language_level': "3"}),
    include_dirs=[numpy.get_include()],
    script_args=["build_ext", "--inplace", "--verbose"]
)
