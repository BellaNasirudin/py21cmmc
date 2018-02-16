"""
The likelihoods package provides a number of pluggable likelihood functions, which take the
output of 21cmFAST and generate a likelihood based on some prescribed data.

Every module defines at least one function, called `get_likelihood`, which must have the following signature::

    get_likelihood(delta_T, box_properties, physical_properties, redshifts, **ll_kwargs)

See the documentation for ``get_likelihood`` in :mod:`~.traditional_c_based` for details on these parameters.
"""
