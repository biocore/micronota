#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
import os
import re
import ast
from setuptools import find_packages, setup
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext as _build_ext

# Bootstrap setup.py with numpy
# Huge thanks to coldfix's solution
# http://stackoverflow.com/a/21621689/579416
class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())

# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('micronota/__init__.py', 'rb') as f:
    hit = _version_re.search(f.read().decode('utf-8')).group(1)
    version = str(ast.literal_eval(hit))

classifiers = [
    'Development Status :: 2 - Pre-Alpha',
    'License :: OSI Approved :: BSD License',
    'Environment :: Console',
    'Topic :: Software Development :: Libraries :: Application Frameworks',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.3',
    'Programming Language :: Python :: 3.4',
    'Operating System :: Unix',
    'Operating System :: POSIX',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: Microsoft :: Windows']


description = 'MICRONOTA: annotation pipeline for microbial (meta)genomes'
with open('README.rst') as f:
    long_description = f.read()

keywords = 'genome metagenome gene annotation RNA',

# Dealing with Cython
USE_CYTHON = os.environ.get('USE_CYTHON', False)
ext = '.pyx' if USE_CYTHON else '.c'

extensions = [
    Extension("micronota.lib.intersection",
              ["micronota/lib/intersection" + ext])
]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)


setup(name='micronota',
      version=version,
      license='BSD',
      description=description,
      long_description=long_description,
      keywords=keywords,
      classifiers=classifiers,
      author="micronota development team",
      author_email="zhenjiang.xu@gmail.com",
      maintainer="micronota development team",
      maintainer_email="zhenjiang.xu@gmail.com",
      url='http://microbio.me/micronota',
      test_suite='nose.collector',
      packages=find_packages(),
      ext_modules=extensions,
      package_data={},
      install_requires=[
          'click >= 6',
          'scikit-bio >= 0.4.0',
          'burrito >= 0.9.1'
      ],
      extras_require={'test': ["nose", "pep8", "flake8"],
                      'coverage': ["coverage"],
                      'doc': ["Sphinx == 1.3.3"]},
      entry_points={
          'console_scripts': [
              'micronota=micronota.cli:cmd',
          ]})
