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

The only subcommand that actually runs anything at this point is ``single``, which runs a single instance of 21cmFAST.
There is also a subcommand ``defaults`` which will print out all the default values of all parameters used within a
run of ``single``. For more information on these parameters, look at the classes ``AstroParams``, ``CosmoParams`` and
``FlagOptions``.

By default, the ``single`` command runs almost exactly the same script as in Brad's March21st version 21CMMC.
It optionally outputs some plots of the obtained fields/lightcones.

At this point, upon installation, the code creates a ``.py21cmmc`` directory in the user's home. In there is the default
config file (you can copy it anywhere and change it if you wish), and a ``Boxes/`` directory which will be where the code
looks for the boxes it needs (this is not true at present).

You can alternatively just run (not yet!)::

    $ py21cmmc init
    $ py21cmmc perturb_field

which will run basically the same things as ``./init`` and ``./perturb_field`` respectively, except that the latter
doesn't need the redshift specification, since this is done in the config file.

Library
~~~~~~~
Typically the user will want to use ``py21cmmc`` as a library -- calling underlying C routines, and obtaining nicely
wrapped results that are ready for further analysis/plotting. The main wrappers around the underlying C are contained
in the ``_21cmfast`` package::

    >>> from py21cmmc._21cmfast import drive_21cmMC

This is where the wrappers of the various structures that are passed in/out of the C live. Typically, these will not
need to be accessed directly. The output of the ``drive_21cmMC`` function is a ``Lightcone`` object (though it need not
contain a lightcone, it could contain a series of co-eval boxes). The contents of this object should be reasonably
self-describing.

TODO (Before Release/Merging with master)
=========================================
- change all input parameters for driver to actually just access the structs.
- functions in drive_21cmmc_streamlined need to use less globals, so they can be called individually.
- users should be able to just call init, or perturb_field etc. from command line.
- the places where Boxes are accessed by default in C code needs to be changed.
- really need to make stuff like BOX_LEN variable so the user doesn't have to recompile!
- Add all the MCMC stuff.
- Documentation
- Tests
- All the other little bits that polish the package off.

Documentation
=============

https://py21cmmc.readthedocs.io/

Development
===========

The layout of the package
~~~~~~~~~~~~~~~~~~~~~~~~~
I should try to explain how I've gone about modifying this to its current state, from Brad's original version.
This version is installable, and the Makefile stuff happens in ``setup.py``. The CLI commands live in cli.py.
At the moment, the ``likelihoods`` subpackage is not being used, but eventually it should hold the code skeletons
that compute likelihoods given the lightcone objects. Funnily enough, plotting functions live in ``plotting.py``,
and at the moment, each of them takes as its main argument a ``Lightcone`` object, which is basically what is
passed back from the C driver. The ``_utils.py`` module contains higher-level wrappers around the C code, but actually
at the moment it is not really very useful, as I have kept the wrapping quite tight, and most of the work is done
in the ``_21cmfast`` package.

The actual Python wrappers of the C, at its basic level, are found in ``_21cmfast/__init__.py`` (so as to make them
importable directly by ``from _21cmfast import wrappers`` etc.). We could maybe change this in future and have a
dedicated ``wrappers.py`` module. All the C code lives in this folder and is compiled by ``setup.py`` from here (this
required changing some of the includes in the C files).

As for input parameters to the functions, I've used a series of Structure classes (I've subclassed each of them to give
defaults for each parameter, so the user doesn't have to worry about most of them). How these work should hopefully be
reasonably clear from the code. The output is also a Structure (I think this could be better). My overall goal is to wrap
as small bits of the C code as possible, in a modular way, and in this module, do nothing fancy with them except return
them in a sensible fashion. The higher-level analysis of these objects should be done outside of this sub-package (say
in ``_utils`` or ``plotting``).

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
