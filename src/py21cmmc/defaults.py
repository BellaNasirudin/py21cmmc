"""
Provides default values for various parameters of 21CMMC. These can be temporarily reset for a given context by
modifying them in place in this module, after import, eg::

    >>> from py21cmmc import defaults
    >>> defaults.BoxDim['HII_DIM'] = 150
    >>> ...
"""

import yaml
from os import path

with open(path.expanduser(path.join("~",'.py21cmmc','defaults.yml'))) as f:
    df = yaml.load(f)

BoxDim = df['BoxDim']
AstroParams = df['AstroParams']
CosmoParams = df['CosmoParams']
FlagOptions = df['FlagOptions']
