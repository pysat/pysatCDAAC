"""Core library for pysatCDAAC.

This is a library of `pysat` instrument modules and methods designed to support
COSMIC instruments and missions archived at the CDAAC web portal.

"""

import importlib
import importlib_metadata

from pysatCDAAC import instruments  # noqa F401

# set version
try:
    __version__ = importlib.metadata.version('pysatNASA')
except AttributeError:
    # Python 3.6 requires a different version
    __version__ = importlib_metadata.version('pysatNASA')
