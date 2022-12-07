# -*- coding: utf-8 -*-
"""Supports IVM data from the COSMIC2 satellite.

The Constellation Observing System for Meteorology, Ionosphere, and Climate
(COSMIC) is comprised of six satellites in LEO with GPS receivers. The
occultation of GPS signals by the atmosphere provides a measurement of
atmospheric parameters. Data downloaded from the COSMIC Data Analaysis
and Archival Center.

More info about the dataset can be found at https://doi.org/10.5065/t353-c093

Properties
----------
platform
    'cosmic2'
name
    'ivm'
tag
    None supported
inst_id
    ['e1', 'e2', 'e3', 'e4', 'e5', 'e6']

Warnings
--------
- Routine was not produced by COSMIC team

"""

import datetime as dt
import functools
import os

import pysat
from pysat.instruments.methods import general as mm_gen
from pysatCDAAC.instruments.methods import general as mm_cdaac

# ----------------------------------------------------------------------------
# Instrument attributes

platform = 'cosmic2'
name = 'ivm'

tags = {'': 'Ion Velocity Meter data'}

inst_ids = {inst_id: list(tags.keys()) for inst_id in ['e1', 'e2', 'e3', 'e4',
                                                       'e5', 'e6']}

# Because all data products are stored in one tar file, inst_id not used
directory_format = os.path.join('{platform}', '{name}', '{tag}')

# ----------------------------------------------------------------------------
# Instrument test attributes

_test_dates = {inst_id: {'': dt.datetime(2021, 1, 1)}
               for inst_id in inst_ids.keys()}

# ----------------------------------------------------------------------------
# Instrument methods


def init(self):
    """Initialize the Instrument object with instrument specific values.

    Note
    ----
    Runs once upon instantiation.

    """
    ack = ' '.join(('FORMOSAT-7/COSMIC-2 is a partnership between the National',
                    'Space Organization in Taiwan and NOAA, the U.S. Air',
                    'Force, and the University Corporation for Atmospheric',
                    'Research (UCAR) in the United States.'))
    refs = ' '.join(('UCAR COSMIC Program, 2019: COSMIC-2 Data Products [Data',
                     'set]. UCAR/NCAR - COSMIC, Access date [insert date],'
                     'https://doi.org/10.5065/T353-C093'))
    self.acknowledgements = ack
    self.references = refs
    pysat.logger.info(ack)
    pysat.logger.info(' '.join(('All spacecraft data stored in single tar',
                                'files on the CDAAC website. Pysat will',
                                'download all and update accordingly.')))

    return


# No cleaning, use standard warning function instead
clean = mm_cdaac.clean_warn


def load(fnames, tag='', inst_id=''):
    """Load COSMIC2 IVM data into `pandas.DataFrame` and `pysat.Meta` objects.

    This routine is called as needed by pysat. It is not intended
    for direct user interaction.

    Parameters
    ----------
    fnames : array-like
        iterable of filename strings, full path, to data files to be loaded.
        This input is nominally provided by pysat itself.
    tag : str
        tag name used to identify particular data set to be loaded.
        This input is nominally provided by pysat itself. (default='')
    inst_id : str
        Satellite ID used to identify particular data set to be loaded.
        This input is nominally provided by pysat itself. (default='')

    Returns
    -------
    data : pds.DataFrame
        A pandas DataFrame with data prepared for the pysat.Instrument
    meta : pysat.Meta
        Metadata formatted for a pysat.Instrument object.

    Note
    ----
    Any additional keyword arguments passed to pysat.Instrument
    upon instantiation are passed along to this routine.

    Examples
    --------
    ::

        inst = pysat.Instrument('cosmic2', 'ivm', inst_id='e1')
        inst.load(2020, 1)

    """

    data, meta = pysat.utils.io.load_netcdf(
        fnames, epoch_name='time', epoch_unit='s',
        epoch_origin=dt.datetime(1980, 1, 6),
        meta_translation={})

    return data, meta


# Use general pysat routine for list files.
fname = 'ivmL2m_C2{id}.{{year:04d}}.{{day:03d}}.??_????.{{version:04d}}_nc'
supported_tags = {inst_id: {'': fname.format(id=inst_id.upper())}
                  for inst_id in inst_ids.keys()}
list_files = functools.partial(mm_gen.list_files, supported_tags=supported_tags)

# Use general CDAAC routine for download.
download_tags = {
    inst_id: {'': {
        'remote_dir': 'gnss-ro/cosmic2/postProc/level2/{year:4d}/{day:03d}/',
        'tar_name': 'ivmL2m_postProc_{year:4d}_{day:03d}.tar.gz',
        'fname': supported_tags[inst_id]['']}}
    for inst_id in inst_ids.keys()}
download = functools.partial(mm_cdaac.download, supported_tags=download_tags)
