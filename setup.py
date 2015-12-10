#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re
import ast
from setuptools import find_packages, setup


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
      package_data={},
      install_requires=[
          'click >= 5.1',
          'scikit-bio >= 0.4.0',
          'burrito >= 0.9.1'
      ],
      extras_require={'test': ["nose", "pep8", "flake8"],
                      'doc': ["Sphinx == 1.2.2", "sphinx-bootstrap-theme"]},
      entry_points={
          'console_scripts': [
              'micronota=micronota.cli:cmd',
          ]})
