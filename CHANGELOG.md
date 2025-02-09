# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.0.5] - 2025-02-06
* Maintenance
  * Update GitHub Actions standards, including SPEC-0 tests
  * Add compatibility for numpy version>=3.2.0
  * Update usage of 'Dataset.dims' to 'Dataset.sizes'
  * Update compatibility with pysat 3.2.1
  * Set minimum pysat version to 3.2.1
  * Set minimum python version to 3.9
  * Update operational environment
  * Update controlled information statement for accuracy and clarity

## [0.0.4] - 2023-08-11
* Bug fixes
  * Update metadata transfer for COSMIC GPS
* Maintenance
  * Update GitHub Actions standards
  * Added support for readthedocs
  * Cleaned up contributing guidelines based on latest project standards
  * Update links in documentation
  * Implement pyproject.toml

## [0.0.3] - 2022-12-12
* Updated `cosmic_gps` to support xarray Datasets
* Added altitude binning profile support for all datasets with altitude
  as a variable.
* Added `download` to general methods
* Added support for `cosmic2_ivm`
* Maintenance
  * Update links in documentation
  * Update testing compliance
  * Add windows testing
  * Add workflow for testing with pysat RC
  * Updated instrument tests to use pysat 3.0.2 syntax
  * Added clean warning routine to `methods.general`

## [0.0.2] - 2021-07-07
* Update instrument style for pysat 3.0.0
* `cosmic_gps` now supports downloads from new public data location
* Update `cosmic_gps` data directory structure
* Updated parsing of `cosmic_gps` filenames to ensure only the latest
  version on disk is retained.

## [0.0.1] - 2020-08-13
* initial port of existing routines from pysat
