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

Then just do::

    pip install py21cmmc

Or to get the bleeding edge::

    pip install git+git://github.com/steven-murray/py21cmmc.git

Quick Usage
===========

Either use the functions from ``py21cmmc._utils``, or use the provided CLI::

    $ py21cmmc --help

At this point, the basic workflow is to first run ``./init`` and ``./perturb_field`` using a separate working ``21cmFAST``
installation. Do all the redshifts that are required (with ``./perturb_field``). Then you need to create a config file -- an
example is given in this repo as ``example_config.yml``. Make sure all the settings in that file are as you want them
(especially the directories containing the boxes you just made). Then use the CLI as written above, with the config
file as an argument. If the ``output`` option is specified, the full deltaT box will be written to that file (in fact, it
will be several deltaT boxes, one for each redshift).

TODO (Before Release)
=====================
- Work out how to change global variables from Python, if possible (otherwise, its hard to re-compile when its installed).
    - This can be done by changing the ``#define VAR (type) val;`` in the appropriate .H file to be ``type VAR=val;``,
      which means it can be changed globally). I have done this for obvious variables like HII_DIM and BOX_LEN in INIT_PARAMS.H,
      which can now be changed. We can add whatever variables we want to this list. *But* we should remember that variables
      that need to be modified in MCMC will need to be more local than that, because each chain shares the object, so one
      of them updating a global var will screw the rest up.
- Add all the MCMC stuff.
- Change to newer version of 21CMFAST
- Add ./init and ./perturb_field (at least) so the whole process can be done without installing anything else.
- Figure out how to use git submodules to do that^.
- Documentation
- Tests
- All the other little bits that polish the package off.

Documentation
=============

https://py21cmmc.readthedocs.io/

Development
===========

To run the all tests run::

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
