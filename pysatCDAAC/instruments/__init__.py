# -*- coding: utf-8 -*-
"""Collection of instruments for the pysatCDAAC library.

Each instrument is contained within a subpackage of this set.

"""

__all__ = ['cosmic_gps', 'cosmic2_ivm']

for inst in __all__:
    exec("from pysatCDAAC.instruments import {x}".format(x=inst))

# Remove dummy variable
del inst
