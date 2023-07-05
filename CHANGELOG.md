# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.0.X] - 2023-XX-XX
* Bug fixes
  * Update metadata transfer for cosmic gps
* Maintenance
  * Update Github Actions standards
  * Added support for readthedocs
  * Cleaned up contributing guidelines based on latest project standards

## [0.0.3] - 2022-12-12
* Updated `cosmic_gps` to support xarray Datasets
* Added altitude binning profile support for all datasets with altitude
  as a variable.
* Added `download` to general methods
* Update links in documentation
* Added support for `cosmic2_ivm`
* Maintenance
  * Update links in documentation
  * Update testing compliance
  * Add windows testing
  * Add workflow for testing with pysat RC
  * Updated instrument tests to use pysat 3.0.2 syntax
  * Added clean warning routine to `methods.general`
  * Implement pyproject.toml

## [0.0.2] - 2021-07-07
* Update instrument style for pysat 3.0.0
* `cosmic_gps` now supports downloads from new public data location
* Update `cosmic_gps` data directory structure
* Updated parsing of `cosmic_gps` filenames to ensure only the latest
  version on disk is retained.

## [0.0.1] - 2020-08-13
* initial port of existing routines from pysat
