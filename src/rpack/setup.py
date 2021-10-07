
"""
setup.py file
"""

from setuptools import setup, Extension
from distutils import sysconfig
import os

packages = ['numpy']
for package in packages:
    try:
        __import__(package)
    except ImportError:
        print("Please install", package)
        exit()

import numpy

sources_spglib = ['arithmetic.c',
                  'cell.c',
                  'delaunay.c',
                  'determination.c',
                  'hall_symbol.c',
                  'kgrid.c',
                  'kpoint.c',
                  'mathfunc.c',
                  'niggli.c',
                  'overlap.c',
                  'pointgroup.c',
                  'primitive.c',
                  'refinement.c',
                  'sitesym_database.c',
                  'site_symmetry.c',
                  'spacegroup.c',
                  'spin.c',
                  'spg_database.c',
                  'spglib.c',
                  'symmetry.c']

source_dir = "../spglib_src"
include_dirs = [source_dir, ]
for i, s in enumerate(sources_spglib):
    sources_spglib[i] = "%s/%s" % (source_dir, s)

rigid_press = Extension('_rigid_press',
                  include_dirs= ['./', numpy.get_include() ],
                  sources=["rigid_press.i", "rigid_press.c", "symmetrization.c", "d_algebra.c", "../randomgen.c"]+sources_spglib,
                  extra_compile_args=["-std=gnu99", "-fPIC", "-O3", "-DROPT_DEBUG"], extra_link_args=["-llapack", "-lblas"])

setup (name = 'rigid_press',
       version = '0.1',
       author      = "Rithwik Tom",
       description = """email:rtom@andrew.cmu.edu""",
       ext_modules = [rigid_press],
       py_modules = ["rigid_press"],
       )
