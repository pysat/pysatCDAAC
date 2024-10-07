<div align="center">
        <img height="0" width="0px">
        <img width="20%" src="https://raw.githubusercontent.com/pysat/pysatCDAAC/main/docs/figures/logo.png" alt="pysat" title="pysatCDAAC"</img>
</div>

# pysatCDAAC: pysat support for COSMIC Data Analysis and Archive Center instruments
[![PyPI Package latest release](https://img.shields.io/pypi/v/pysatCDAAC.svg)](https://pypi.python.org/pypi/pysatCDAAC)
[![Build Status](https://github.com/pysat/pysatCDAAC/actions/workflows/main.yml/badge.svg)](https://github.com/pysat/pysatCDAAC/actions/workflows/main.yml/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/pysat/pysatCDAAC/badge.svg?branch=main)](https://coveralls.io/github/pysat/pysatCDAAC?branch=main)
[![DOI](https://zenodo.org/badge/287322558.svg)](https://zenodo.org/badge/latestdoi/287322558)

### Prerequisites

pysatCDAAC uses common Python modules, as well as modules developed by
and for the Space Physics community.  This module officially supports
Python 3.7+.

| Common modules | Community modules   |
| -------------- | ------------------- |
| netCDF4        | pysat>=3.1.0        |
| numpy          |                     |
| pandas         |                     |
| requests       |                     |
| xarray         |                     |

# Installation

Currently, the main way to get pysatCDAAC is through github.

```
git clone https://github.com/pysat/pysatCDAAC.git
```

Change directories into the repository folder and run the setup.py file.  For
a local install use the "--user" flag after "install".

```
cd pysatCDAAC/
pip install .
```

Note: pre-1.0.0 version
-----------------------
pysatCDAAC is currently in an initial development phase and requires pysat 3.0.0+.


# Using with pysat

The instrument modules are portable and designed to be run like any pysat instrument.

```
import pysat
from pysatCDAAC.instruments import cosmic_gps

cosmic = pysat.Instrument(inst_module=cosmic_gps, tag='ionprf')
```
Another way to use the instruments in an external repository is to register the instruments.  This only needs to be done the first time you load an instrument.  Afterward, pysat will identify them using the `platform` and `name` keywords.

```
import pysat

pysat.utils.registry.register(‘pysatCDAAC.instruments.cosmic_gps’)
cosmic = pysat.Instrument('cosmic', 'gps', tag='ionprf')
```
