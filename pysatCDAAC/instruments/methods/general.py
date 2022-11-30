# -*- coding: utf-8 -*-
"""Provides default routines for CDAAC instruments into pysat."""

import fnmatch
import importlib
import os
import requests
import tarfile
import tempfile

import pysat
from pysat.utils.files import construct_searchstring_from_format


def download(date_array, tag, inst_id, supported_tags=None, data_path=None,
             sub_path=False, sort_files=False, user=None, password=None):
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
        Path to directory to download data to. (default=None)
    sub_path : bool
        If True, break up data into further subdirectories based on date.
        (default=False)
    sort_files : bool
        If True, sorts files by inst_id after all are downloaded.
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

    if tag is None:
        tag = ''
    if inst_id is None:
        inst_id = ''
    try:
        inst_dict = supported_tags[inst_id][tag]
    except KeyError:
        raise ValueError('inst_id / tag combo unknown.')

    for date in date_array:
        pysat.logger.info('Downloading COSMIC data for ' + date.strftime('%D'))
        yr, day = pysat.utils.time.getyrdoy(date)
        yrdoystr = '{year:04d}/{day:03d}'.format(year=yr, day=day)

        # Try re-processed data (preferred).
        # Construct path string for online file.
        dwnld = ''.join(("https://data.cosmic.ucar.edu/",
                         inst_dict['remote_dir'],
                         inst_dict['tar_name']))
        dwnld = dwnld.format(year=yr, day=day)
        try:
            # Make online connection.
            with requests.get(dwnld) as req:
                req.raise_for_status()
        except requests.exceptions.HTTPError:
            # If response is negative, try post-processed data
            # Construct path string for online file
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
            # Note: sort_files only used for cosmic2 ivm
            if sort_files:
                inst_ids = ['e1', 'e2', 'e3', 'e4', 'e5', 'e6']
                fnames = tar.getnames()
                for id in inst_ids:
                    search_str, final_path = get_instrument_data_info(
                        'cosmic2_ivm', tag=tag, inst_id=id)
                    matched_files = fnmatch.filter(fnames, search_str)
                    for fname in matched_files:
                        tar.extract(fname, path=final_path)
            else:
                if sub_path:
                    final_path = os.path.join(data_path, yrdoystr)
                else:
                    final_path = data_path
                tar.extractall(path=final_path)
            tar.close()

        except tarfile.ReadError:
            # If file cannot be read as a tarfile, then data does not exist.
            # Skip this day since there is nothing left to do.
            pass

    # Remove the temporary directory (even if download fails)
    temp_dir.cleanup()

    return


def get_instrument_data_info(inst_mod_name, tag='', inst_id='', **kwargs):
    """Get the search string and `data_path` attribute from an Instrument module.

    Parameters
    ----------
    inst_mod_name : str
        pysatCDAAC Instrument module name
    tag : str
        String specifying the Instrument tag (default='')
    inst_id : str
        String specifying the instrument identification (default='')
    kwargs : dict
        Optional additional kwargs that may be used to initialize an Instrument

    Returns
    -------
    search_str : str
        Pattern for instrument file names
    data_path : str
        Path where the Instrument data is stored

    """

    # Import the desired instrument module by name
    inst_mod = importlib.import_module(".".join(["pysatCDAAC",
                                                 "instruments", inst_mod_name]))

    # Construct a search string to identify data files
    search_dict = construct_searchstring_from_format(
        inst_mod.supported_tags[inst_id][tag])
    search_str = search_dict['search_string']

    # Initialize a temporary instrument to obtain pysat configuration
    temp_inst = pysat.Instrument(inst_module=inst_mod, tag=tag, inst_id=inst_id,
                                 **kwargs)

    # Save the data path for this Instrument down to the inst_id level
    data_path = temp_inst.files.data_path

    # Delete the temporary instrument
    del temp_inst

    return search_str, data_path
