from setuptools import setup
from Cython.Build import cythonize
import numpy as np
import os

print(np.get_include())

os.environ["C_INCLUDE_PATH"] = np.get_include()
setup(
    ext_modules = cythonize(curr_dir +'functions/fastFunctions.pyx')
)

