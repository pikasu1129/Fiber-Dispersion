# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 19:57:42 2018

@author: arne
"""

#from distutils.core import setup
#from distutils.extension import Extension
#from Cython.Distutils import build_ext
import numpy

#ext_modules=[ Extension("speedfiber",
              #["speedfiber.pyx"],
              #include_dirs=[numpy.get_include()],
              #libraries=["m","fftw3"],
              #extra_compile_args = ["-ffast-math", "-fopenmp", "-lfftw3 -lm"],
              #extra_link_args=['-fopenmp'])]

#setup(
  #name = "speedfiber",
  #cmdclass = {"build_ext": build_ext},
#  ext_modules = ext_modules)

#from distutils.core import setup
from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize

#from distutils.extension import Extension
#from Cython.Distutils import build_ext
#ext_modules

extensions = [Extension("speedfiber",
                         ["speedfiber.pyx"],
                           include_dirs=[numpy.get_include(), '/usr/include/', '/usr/local/include/include'], # not needed for fftw unless it is installed in an unusual place
                            libraries=['m', 'fftw3','fftw3_omp'],
                        library_dirs=['/usr/lib/'], # not needed for fftw unless it is installed in an unusual place 
                         extra_compile_args=['-ffast-math', '-lfftw3', '-lm', '-fopenmp'], # '-fopenmp', lpthread
                        extra_link_args=['-fopenmp'])]


setup(
  name = 'Speedfiber',
  #cmdclass = {'build_ext': build_ext},
  packages = find_packages(),
  #ext_modules = ext_modules
  ext_modules = cythonize(extensions)
)
