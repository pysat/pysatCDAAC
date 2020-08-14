# pysatCDAAC: pysat support for COSMIC Data Analysis and Archive Center instruments
[![Build Status](https://travis-ci.org/pysat/pysatCDAAC.svg?branch=main)](https://travis-ci.org/pysat/pysatCDAAC)
[![Coverage Status](https://coveralls.io/repos/github/pysat/pysatCDAAC/badge.svg?branch=main)](https://coveralls.io/github/pysat/pysatCDAAC?branch=main)

# Installation

Currently, the main way to get pysatCDAAC is through github.

```
git clone https://github.com/pysat/pysatCDAAC.git
```

Change directories into the repository folder and run the setup.py file.  For
a local install use the "--user" flag after "install".

```
cd pysatCDAAC/
python setup.py install
```

Note: pre-1.0.0 version
------------------
pysatCDAAC is currently in an initial development phase.  Much of the API is being built off of the upcoming pysat 3.0.0 software in order to streamline the usage and test coverage.  This version of pysat is planned for release later this year.  Currently, you can access the develop version of this through github:
```
git clone https://github.com/pysat/pysat.git
cd pysat
git checkout develop-3
python setup.py install
```
It should be noted that this is a working branch and is subject to change.

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
