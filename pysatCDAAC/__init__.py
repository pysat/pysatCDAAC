# -*- coding: utf-8 -*-
# Full author list can be found in .zenodo.json file
# DOI:10.5281/zenodo.3475493
#
# Review Status for Classified or Controlled Information by NRL
# -------------------------------------------------------------
# DISTRIBUTION STATEMENT A: Approved for public release. Distribution is
# unlimited.
# ----------------------------------------------------------------------------
"""Core library for pysatCDAAC.

This is a library of `pysat` instrument modules and methods designed to support
COSMIC instruments and missions archived at the CDAAC web portal.

"""

import importlib
from pysatCDAAC import instruments  # noqa F401

# Set version
__version__ = importlib.metadata.version('pysatCDAAC')
