# -*- coding: utf-8 -*-
# Full author list can be found in .zenodo.json file
# DOI:10.5281/zenodo.3475493
#
# DISTRIBUTION STATEMENT A: Approved for public release. Distribution is
# unlimited.
# ----------------------------------------------------------------------------
"""Collection of instruments for the pysatCDAAC library.

Each instrument is contained within a subpackage of this set.

"""

__all__ = ['cosmic_gps', 'cosmic2_ivm']

for inst in __all__:
    exec("from pysatCDAAC.instruments import {x}".format(x=inst))

# Remove dummy variable
del inst
