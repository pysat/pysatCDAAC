# -*- coding: utf-8 -*-
"""Provides default routines for CDAAC instruments into pysat."""

import os
import requests
import tarfile
import tempfile
import warnings

import pysat


def download(date_array, tag, inst_id, supported_tags=None,
             data_path=None, sub_path=False, user=None, password=None):
    """Download data from CDAAC https server.

    Parameters
    ----------
    date_array : array-like
        list of datetimes to download data for. The sequence of dates need not
        be contiguous.
    tag : str
        Tag identifier used for particular dataset. This input is provided by
        pysat. (default='')
    inst_id : str
        Instrument ID string identifier used for particular dataset. This input
        is provided by pysat. (default='')
    supported_tags : dict
        dict of dicts. Keys are supported tag names for download. Value is
        a dict with 'remote_dir', 'tar_name'. Intended to be pre-set with
        functools.partial then assigned to new instrument code.
        (default=None)
    data_path : str
        Path to directory to download data to. (default='')
    sub_path : bool
        If True, break up data into further subdirectories based on date.
        (default=False)
    user : str or NoneType
        User string input used for download. Provided by user and passed via
        pysat. If an account is required for downloads this routine here must
        error if user not supplied. (default=None)
    password : str or NoneType
        Password for data download. (default=None)

    Note
    ----
    This routine is invoked by pysat and is not intended for direct use by
    the end user.

    """

    # Set up temporary directory for tar files
    temp_dir = tempfile.TemporaryDirectory()

    inst_dict = supported_tags[inst_id][tag]

    for date in date_array:
        pysat.logger.info('Downloading COSMIC data for ' + date.strftime('%D'))
        yr, day = pysat.utils.time.getyrdoy(date)
        yrdoystr = '{year:04d}/{day:03d}'.format(year=yr, day=day)

        # Try re-processed data (preferred). Construct a path string for the
        # online file.
        dwnld = ''.join(("https://data.cosmic.ucar.edu/",
                         inst_dict['remote_dir'],
                         inst_dict['tar_name']))
        dwnld = dwnld.format(year=yr, day=day)
        try:
            # Make online connection.
            with requests.get(dwnld) as req:
                req.raise_for_status()
        except requests.exceptions.HTTPError:
            # If response is negative, try post-processed data. Construct
            # a path string for the online file
            if 'backup' in inst_dict.keys():
                # If a backup exists, try alternate form
                dwnld = dwnld.replace(inst_dict['backup'][0],
                                      inst_dict['backup'][1])
                try:
                    # Make online connection
                    with requests.get(dwnld) as req:
                        req.raise_for_status()
                except requests.exceptions.HTTPError as err:
                    estr = ''.join((str(err), '\n', 'Data not found'))
                    pysat.logger.info(estr)

        # Copy request info to tarball file with generated name in `fname`.
        fname = os.path.join(temp_dir.name, inst_dict['tar_name'])
        fname = fname.format(year=yr, day=day)
        with open(fname, "wb") as local_file:
            local_file.write(req.content)
            local_file.close()
        try:
            # Uncompress files and remove tarball
            tar = tarfile.open(fname)
            if sub_path:
                # Send to subdirectory.
                tar.extractall(path=os.path.join(data_path, yrdoystr))
            else:
                # Send to top level.
                tar.extractall(path=data_path)
            tar.close()

        except tarfile.ReadError:
            # If file cannot be read as a tarfile, then data does not exist.
            # Skip this day since there is nothing left to do.
            pass

    # Remove the temporary directory (even if download fails)
    temp_dir.cleanup()

    return


def clean_warn(self):
    """Warn user that cleaning not yet available for this data set.

    Note
    ----
    'clean' - Not specified
    'dusty' - Not specified
    'dirty' - Not specified
    'none'  No cleaning applied, routine not called in this case.

    """
    warnings.warn(' '.join(('No cleaning routines available for',
                            self.platform, self.name)))

    return
