import sys

import numpy as np
from setuptools import setup, Extension

import Cython
from Cython.Distutils import build_ext
from Cython.Build import cythonize

try:
    numpy_include = np.get_include()
except AttributeError:
    numpy_include = np.get_numpy_include()

# check Python version
if not (sys.version_info[0] == 2 and sys.version_info[1] >= 6):
    sys.exit("You need Python 2.6.x or Python 2.7.x to install the DockBox package!")

ext_modules = [Extension(
    name='dockbox.pyqcprot',
    sources=["dockbox/pyqcprot.pyx"],
    include_dirs=[numpy_include])]

setup(name='dockbox',
    version='1.4',
    packages=['dockbox'],
    scripts=['bin/rundbx', 'bin/extract_top_poses'],
    install_requires=['mdkit', 'pandas<=0.24.2', 'nwalign', 'oldnumeric'],
    ext_modules = cythonize(ext_modules),
    license='LICENSE.txt',
    description='Platform package to simplify the use of docking programs and consensus methods',
    long_description=open('README.rst').read(),
)
