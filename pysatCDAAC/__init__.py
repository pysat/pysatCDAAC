"""Core library for pysatCDAAC.

This is a library of `pysat` instrument modules and methods designed to support
COSMIC instruments and missions archived at the CDAAC web portal.

"""

import importlib
from pysatCDAAC import instruments  # noqa F401

# Set version
try:
    __version__ = importlib.metadata.version('pysatCDAAC')
except AttributeError:
    # Python 3.6 requires a different version
    import importlib_metadata
    __version__ = importlib_metadata.version('pysatCDAAC')
