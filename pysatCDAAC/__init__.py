"""Core library for pysatCDAAC.

This is a library of `pysat` instrument modules and methods designed to support
COSMIC instruments and missions archived at the CDAAC web portal.

"""

import pkg_resources

from pysatCDAAC import instruments  # noqa F401

__version__ = pkg_resources.get_distribution('pysatCDAAC').version
