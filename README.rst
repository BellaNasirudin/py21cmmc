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

A set of python bindings for 21cmFAST allowing native Python plugins for 21cmMC

* Free software: MIT license

Installation
============

::

    pip install py21cmmc

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
