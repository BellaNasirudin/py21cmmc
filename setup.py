#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function

import io
import os
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import relpath
from os.path import splitext
from os.path import expanduser

from setuptools import Extension
from setuptools import find_packages
from setuptools import setup

from shutil import copyfile, move

def read(*names, **kwargs):
    return io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ).read()

# ======================================================================================================================
# Create a user-level config directory for py21cmmc, just so we can have reasonable default storage areas.
try:
    pkgdir = expanduser(join("~", ".py21cmmc"))
    os.mkdir(pkgdir)
    os.mkdir(join(pkgdir, "Boxes"))
except:
    pass

try:
    move(join(pkgdir, "example_config.yml"), join(pkgdir, "example_config.yml.bk"))
except:
    pass

copyfile("example_config.yml", join(pkgdir, "example_config.yml"))
# ======================================================================================================================


# Enable code coverage for C code: we can't use CFLAGS=-coverage in tox.ini, since that may mess with compiling
# dependencies (e.g. numpy). Therefore we set SETUPPY_CFLAGS=-coverage in tox.ini and copy it to CFLAGS here (after
# deps have been safely installed).
if 'TOXENV' in os.environ and 'SETUPPY_CFLAGS' in os.environ:
    os.environ['CFLAGS'] = os.environ['SETUPPY_CFLAGS']

setup(
    name='py21cmmc',
    version='0.1.0',
    license='MIT license',
    description='A set of python bindings for 21cmFAST allowing native Python plugins for 21cmMC',
    long_description='%s\n%s' % (
        re.compile('^.. start-badges.*^.. end-badges', re.M | re.S).sub('', read('README.rst')),
        re.sub(':[a-z]+:`~?(.*?)`', r'``\1``', read('CHANGELOG.rst'))
    ),
    author='Steven Murray',
    author_email='steven.murray@curtin.edu.au',
    url='https://github.com/steven-murray/py21cmmc',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
        # uncomment if you test on these interpreters:
        # 'Programming Language :: Python :: Implementation :: IronPython',
        # 'Programming Language :: Python :: Implementation :: Jython',
        # 'Programming Language :: Python :: Implementation :: Stackless',
        'Topic :: Utilities',
    ],
    keywords=[
        # eg: 'keyword1', 'keyword2', 'keyword3',
    ],
    install_requires=[
        'click',
        'tqdm',
        'numpy',
        'pyyaml'
        # eg: 'aspectlib==1.1.1', 'six>=1.7',
    ],
    extras_require={
        # eg:
        #   'rst': ['docutils>=0.11'],
        #   ':python_version=="2.6"': ['argparse'],
    },
    entry_points={
        'console_scripts': [
            'py21cmmc = py21cmmc.cli:main',
        ]
    },
    ext_modules=[
        Extension(
            'py21cmmc._21cmfast.drive_21cmMC_streamlined',
            sources=['src/py21cmmc/_21cmfast/drive_21cmMC_streamlined.c'],
            libraries=['m', 'gsl', 'gslcblas', 'fftw3f_omp', 'fftw3f'],
            include_dirs=['/usr/local/include', 'src/py21cmmc/_21cmfast'],
            extra_compile_args = ['-fopenmp', '-Ofast']
        ),
        # Extension(
        #     'py21cmmc._21cmfast.init',
        #     sources=['src/py21cmmc/_21cmfast/init.c'],
        #     libraries=['m', 'gsl', 'gslcblas', 'fftw3f_omp', 'fftw3f'],
        #     include_dirs=['/usr/local/include', 'src/py21cmmc/_21cmfast'],
        #     extra_compile_args = ['-fopenmp', '-Ofast']
        # ),
        # Extension(
        #     'py21cmmc._21cmfast.perturb_field',
        #     sources=['src/py21cmmc/_21cmfast/perturb_field.c'],
        #     libraries=['m', 'gsl', 'gslcblas', 'fftw3f_omp', 'fftw3f'],
        #     include_dirs=['/usr/local/include', 'src/py21cmmc/_21cmfast'],
        #     extra_compile_args = ['-fopenmp', '-Ofast']
        # )

    ],
)
