# -*- coding: utf-8 -*-
"""Supports IVM data from the COSMIC2 satellite.

The Constellation Observing System for Meteorology, Ionosphere, and Climate
(COSMIC) is comprised of six satellites in LEO with GPS receivers. The
occultation of GPS signals by the atmosphere provides a measurement of
atmospheric parameters. Data downloaded from the COSMIC Data Analaysis
and Archival Center.

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
import pandas as pds
import warnings

import pysat
from pysat.instruments.methods import general as mm_gen
from pysatCDAAC.instruments.methods import general as mm_cdaac

# ----------------------------------------------------------------------------
# Instrument attributes

platform = 'cosmic2'
name = 'ivm'

tags = {'': 'Ion Velocity Meter data'}

inst_ids = {'e1': [''], 'e2': [''], 'e3': [''], 'e4': [''], 'e5': [''],
            'e6': ['']}

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
    ack = ' '.join(('Add acknowledgements here')).upper()
    refs = ' '.join(('Add refs here')).upper()
    self.acknowledgements = ack
    self.references = refs
    pysat.logger.info(ack)
    pysat.logger.info(' '.join(('All spacecraft data stored in single tar',
                                'files on the CDAAC website. Pysat will',
                                'download all and update accordingly.')))

    return


def clean(self):
    """Clean COSMIC2 IVM data to the specified level.

    Note
    ----
    'clean' - Not specified
    'dusty' - Not specified
    'dirty' - Not specified
    'none'  No cleaning applied, routine not called in this case.

    """
    warnings.warn('Cleaning not yet implemented')

    return


def load(fnames, tag=None, inst_id=None):
    """Load COSMIC2 IVM data into `pandas.DataFrame` and `pysat.Meta` objects.

    This routine is called as needed by pysat. It is not intended
    for direct user interaction.

    Parameters
    ----------
    fnames : array-like
        iterable of filename strings, full path, to data files to be loaded.
        This input is nominally provided by pysat itself.
    tag : string
        tag name used to identify particular data set to be loaded.
        This input is nominally provided by pysat itself.
    inst_id : string
        Satellite ID used to identify particular data set to be loaded.
        This input is nominally provided by pysat itself.

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

        inst = pysat.Instrument('cosmic2', 'ivm', inst_id='e1', tag='')
        inst.load(2020, 1)

    """

    data, meta = pysat.utils.load_netcdf4(fnames, epoch_name='uts')

    # COSMIC2 uses GPS time, not yet supported in standard pysat load.
    # Update the index seconds since Jan 6, 1980.
    index = (pds.to_timedelta(data['time'].values, unit='s')
             + dt.datetime(1980, 1, 6))
    data.index = index

    return data, meta


fname = 'ivmL2m_C2{id}.{{year:04d}}.{{day:03d}}.??_????.{{version:04d}}_nc'
supported_tags = {inst_id: {'': fname.format(id=inst_id.upper())}
                  for inst_id in inst_ids.keys()}
list_files = functools.partial(mm_gen.list_files, supported_tags=supported_tags)

download_tags = {
    inst_id: {'': {
        'remote_dir': 'gnss-ro/cosmic2/postProc/level2/{year:4d}/{day:03d}/',
        'tar_name': 'ivmL2m_postProc_{year:4d}_{day:03d}.tar.gz',
        'fname': supported_tags[inst_id]['']}}
    for inst_id in inst_ids.keys()}
download = functools.partial(mm_cdaac.download, supported_tags=download_tags)
