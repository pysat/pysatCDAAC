# -*- coding: utf-8 -*-
"""Provides default routines for CDAAC instruments into pysat."""

import os
import requests
import tarfile

import pysat


def download(date_array, tag, inst_id, supported_tags=None, data_path=None,
             user=None, password=None):
    """Download COSMIC GPS data.

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
        a dict with 'remote_dir', 'fname'. Inteded to be pre-set with
        functools.partial then assigned to new instrument code.
        (default=None)
    data_path : str
        Path to directory to download data to. (default=None)
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

    level_str = supported_tags[tag]['level']
    sub_str = supported_tags[tag]['substr']

    for date in date_array:
        pysat.logger.info('Downloading COSMIC data for ' + date.strftime('%D'))
        yr, doy = pysat.utils.time.getyrdoy(date)
        yrdoystr = '{year:04d}/{doy:03d}'.format(year=yr, doy=doy)

        # Try re-processed data (preferred).
        # Construct path string for online file.
        dwnld = ''.join(("https://data.cosmic.ucar.edu/gnss-ro/cosmic1",
                         "/repro2013/", level_str, "/", yrdoystr, "/",
                         sub_str, '_repro2013',
                         '_{year:04d}_{doy:03d}.tar.gz'.format(year=yr,
                                                               doy=doy)))
        try:
            # Make online connection.
            with requests.get(dwnld) as req:
                req.raise_for_status()
        except requests.exceptions.HTTPError:
            # If response is negative, try post-processed data
            # Construct path string for online file
            dwnld = ''.join(("https://data.cosmic.ucar.edu/gnss-ro/cosmic1",
                             "/postProc/", level_str, "/", yrdoystr, "/",
                             sub_str, '_postProc',
                             '_{year:04d}_{doy:03d}.tar.gz'))
            dwnld = dwnld.format(year=yr, doy=doy)
            try:
                # Make online connection
                with requests.get(dwnld) as req:
                    req.raise_for_status()
            except requests.exceptions.HTTPError as err:
                estr = ''.join((str(err), '\n', 'Data not found'))
                pysat.logger.info(estr)

        # Copy request info to tarball file with generated name in `fname`.
        fname = os.path.join(data_path,
                             ''.join(('cosmic_', sub_str,
                                      '_{year:04d}.{doy:03d}.tar')))
        fname = fname.format(year=yr, doy=doy)
        with open(fname, "wb") as local_file:
            local_file.write(req.content)
            local_file.close()
        try:
            # Uncompress files and remove tarball
            tar = tarfile.open(fname)
            tar.extractall(path=os.path.join(data_path, yrdoystr))
            tar.close()

        except tarfile.ReadError:
            # If file cannot be read as a tarfile, then data does not exist.
            # Skip this day since there is nothing left to do.
            pass
        # tar file must be removed (even if download fails)
        os.remove(fname)

    return
