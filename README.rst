========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |travis|
        | |coveralls|
        | |codacy|
    * - package
      - | |version| |wheel| |supported-versions| |supported-implementations|
        | |commits-since|

.. |docs| image:: https://readthedocs.org/projects/py21cmmc/badge/?style=flat
    :target: https://readthedocs.org/projects/py21cmmc
    :alt: Documentation Status

.. |travis| image:: https://travis-ci.org/steven-murray/py21cmmc.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/steven-murray/py21cmmc

.. |coveralls| image:: https://coveralls.io/repos/steven-murray/py21cmmc/badge.svg?branch=master&service=github
    :alt: Coverage Status
    :target: https://coveralls.io/r/steven-murray/py21cmmc

.. |codacy| image:: https://img.shields.io/codacy/REPLACE_WITH_PROJECT_ID.svg
    :target: https://www.codacy.com/app/steven-murray/py21cmmc
    :alt: Codacy Code Quality Status

.. |version| image:: https://img.shields.io/pypi/v/py21cmmc.svg
    :alt: PyPI Package latest release
    :target: https://pypi.python.org/pypi/py21cmmc

.. |commits-since| image:: https://img.shields.io/github/commits-since/steven-murray/py21cmmc/v0.1.0.svg
    :alt: Commits since latest release
    :target: https://github.com/steven-murray/py21cmmc/compare/v0.1.0...master

.. |wheel| image:: https://img.shields.io/pypi/wheel/py21cmmc.svg
    :alt: PyPI Wheel
    :target: https://pypi.python.org/pypi/py21cmmc

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/py21cmmc.svg
    :alt: Supported versions
    :target: https://pypi.python.org/pypi/py21cmmc

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/py21cmmc.svg
    :alt: Supported implementations
    :target: https://pypi.python.org/pypi/py21cmmc


.. end-badges

A set of python bindings for 21cmFAST allowing native Python plugins for 21cmMC.

Note: this package has just started. It doesn't do any MCMC at the moment.

* Free software: MIT license

Installation
============

First, you'll need to have the required C libraries: ``gsl``, ``fftw`` (make sure you install the floating-point version!)
``openmp`` and ``gslcblas``.

Then just do (not actually online yet!)::

    pip install py21cmmc

Or to get the bleeding edge::

    pip install git+git://github.com/steven-murray/py21cmmc.git

For development, it is easiest to do (from top-level directory of this package)::

    pip install -e .

Quick Usage
===========

We try to support two methods of using py21cmmc:

CLI
~~~
The CLI interface always starts with the command ``py21cmmc``, and has a number of subcommands. To list the available
subcommands, use::

    $ py21cmmc --help

To get help on any subcommand, simply use::

    $ py21cmmc <subcommand> --help

The primary subcommand that actually runs things at this point is ``single``, which runs a single instance of 21cmFAST.
There are also (for simplicity) subcommands ``init`` and ``perturb_field``, which operate almost exactly as if those
programs from 21cmFAST had been run. Note that they are redundant however, as ``single`` will run these automatically.

The data files (for init and perturb_field) are by default saved in a ``Boxes/`` directory within the current directory.
This can be changed with the ``outdir`` argument to the CLI script (more generally, it can be set as the ``DIREC`` value
in the ``BoxDimStruct`` to any function which computes 21cmFAST boxes).

There is also a subcommand ``defaults`` which will print out all the default values of all parameters used within a
run of ``single``. For more information on these parameters, look at the classes ``AstroParamStruct``,
``CosmoParamStruct``, ``FlagOptionStruct`` and ``BoxDimStruct`` from the ``wrapper`` module.

By default, the ``single`` command runs almost exactly the same script as in Brad's March21st version 21CMMC.
It optionally outputs some plots of the obtained fields/lightcones.

There is also an ``mcmc`` subcommand, which runs MCMC, using an input config file (by default, it is the
``example_config_mcmc.yml`` file provided in the top-level directory of this repository). This allows a range of options
to be set, and a range of likelihoods to be used, each of which can access the core outputs of 21cmFAST to build the
likelihood. This enables one to insert code which can perform arbitrary analyses on the lightcone, for example, and
use that for the likelihood.
At this point, upon installation, the code creates a ``.py21cmmc`` directory in the user's home. In there is the default
config file (you can copy it anywhere and change it if you wish), and the default defaults for all C options and parameters.


Library
~~~~~~~
Typically the user will want to use ``py21cmmc`` as a library -- calling underlying C routines, and obtaining nicely
wrapped results that are ready for further analysis/plotting. The main wrappers around the underlying C are contained
in the ``wrapper`` module, which are imported into the main namespace::

    >>> from py21cmmc import run_21cmfast, AstroParamStruct,...

This is where the wrappers of the various structures that are passed in/out of the C live. Typically, one should only
really need to use the ``run_21cmfast`` function, unless working in development. The output of the ``run_21cmfast``
function is either a ``LightCone`` object  or a ``CoEval`` object. The contents of this object should be reasonably
self-describing.

The ``likelihood`` module contains various Likelihood classes, written to be useable within the CosmoHammer framework.
One can easily build their own CosmoHammer MCMC framework using these classes, or use the provided function in the ``mcmc``
module to perform MCMC on a select subset of the likelihoods. Refer to the documentation of the ``likelihood`` module
for more details.

TODO (Before Release/Merging with master)
=========================================
- change all input parameters for driver to actually just access the structs.
- functions in drive_21cmmc_streamlined need to use less globals, so they can be called individually.
- remove GeneratePS from main drive_21cmMC function (make it callable by Python).
- Documentation
- Tests
- All the other little bits that polish the package off.

Documentation
=============

To view the docs, install the ``requirements_dev.txt`` packages, go to the docs/ folder, and type "make html", then
open the ``index.html`` file in the ``_build/html`` directory.

Development
===========

The layout of the package
~~~~~~~~~~~~~~~~~~~~~~~~~
I should try to explain how I've gone about modifying this to its current state, from Brad's original version.
This version is installable, and the Makefile stuff happens in ``setup.py``. The CLI commands live in cli.py.
Funnily enough, plotting functions live in ``plotting.py``,
and at the moment, each of them takes as its main argument a ``Lightcone`` or ``CoEval`` object, which is basically what is
passed back from the C driver. The ``_utils.py`` module contains a couple of functions for writing out the parameter
files which can be read in by the C driver. I have not made use of these in the rest of the code, however.

The actual Python wrappers of the C, at its basic level, are found in ``wrapper.py``. All the C code lives in the ``_21cmfast``
folder and is compiled by ``setup.py`` from here (this required changing some of the includes in the C files).

The wrapping is done with CFFI, rather than the native ctypes. This allows for less redundant specification of types
etc. The things to watch out for, when using CFFI, is the memory management. If an array is created in Python, and a
pointer to it is set to a C variable, then that Python variable has to stick around otherwise the memory is effectively
free'd, and weird stuff happens. This is usually obvious, but is sometimes obscured when setting a C variable to the
result of a function call, for which no Python variable has ever been specified (and so it quickly gets garbage collected).

The building of the C code is done in ``build_cffi.py``. At the moment, it's a bit rough, due to the number of global
defines that are used. However, the overall structure is such that ``set_source`` literally just includes the main
source code that needs to be there to run. The ``cdef`` defines the signatures of all global parameters and functions
which ought to be wrapped. This *should* be as easy as including a header file, but #defines only get captured if you
specify them manually as static const, and furthermore, there *is* no header file which contains the main functions we
care about. So they are copied in at this point.

As for input parameters to the functions, I've used a series of Structure classes (I've subclassed each of them to give
defaults for each parameter, so the user doesn't have to worry about most of them). How these work should hopefully be
reasonably clear from the code.


Meta-development stuff
~~~~~~~~~~~~~~~~~~~~~~
I'm using a git-flow git system, where we can create features and fixes etc. If you don't like that, feel free to change
it or discuss it. I think we should use the Github issue system to handle all of our "todo's" and then we can each pick
them off easily, and comment on their viability.

To run the all tests run (no tests as yet...)::

    tox

Note, to combine the coverage data from all the tox environments run:

.. list-table::
    :widths: 10 90
    :stub-columns: 1

    - - Windows
      - ::

            set PYTEST_ADDOPTS=--cov-append
            tox

    - - Other
      - ::

            PYTEST_ADDOPTS=--cov-append tox
