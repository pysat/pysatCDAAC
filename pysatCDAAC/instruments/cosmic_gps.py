# -*- coding: utf-8 -*-
"""
Loads data from the COSMIC satellite.

The Constellation Observing System for Meteorology, Ionosphere, and Climate
(COSMIC) is comprised of six satellites in LEO with GPS receivers. The
occultation of GPS signals by the atmosphere provides a measurement of
atmospheric parameters. Data downloaded from the COSMIC Data Analaysis
and Archival Center.

Default behavior is to search for the 2013 re-processed data first, then the
post-processed data as recommended on
https://cdaac-www.cosmic.ucar.edu/cdaac/products.html

Properties
----------
platform
    'cosmic'
name
    'gps' for Radio Occultation profiles
tag
    Select profile type, or scintillation, one of:
    {'ionprf', 'sonprf', 'wetprf', 'atmprf', 'scnlv1'}
inst_id
    None supported
altitude_bin
    Number of kilometers to bin altitude profiles by when loading.
    Currently only supported for tag='ionprf'.

Note
----
- 'ionprf: 'ionPrf' ionosphere profiles
- 'sonprf': 'sonPrf' files
- 'wetprf': 'wetPrf' files
- 'atmprf': 'atmPrf' files
- 'scnlv1': 'scnLv1' files

Warnings
--------
- Routine was not produced by COSMIC team
- More recent versions of netCDF4 and numpy limit the casting of some variable
  types into others. This issue could prevent data loading for some variables
  such as 'MSL_Altitude' in the 'sonprf' and 'wetprf' files. The default
  UserWarning when this occurs is
  ::

    'UserWarning: WARNING: missing_value not used since it cannot be safely
    cast to variable data type'

"""

import datetime as dt
import os
import requests
import shutil
import sys
import tarfile

import numpy as np
import netCDF4
import pandas as pds
import pysat
from pysat import logger
from pysat.utils import files as futils
import xarray as xr


# ----------------------------------------------------------------------------
# Instrument attributes

platform = 'cosmic'
name = 'gps'
tags = {'ionprf': '',
        'sonprf': '',
        'wetprf': '',
        'atmprf': '',
        'scnlv1': ''}
inst_ids = {'': ['ionprf', 'sonprf', 'wetprf', 'atmprf', 'scnlv1']}

pandas_format = False

# ----------------------------------------------------------------------------
# Instrument test attributes

_test_dates = {'': {'ionprf': dt.datetime(2008, 1, 1),
                    'sonprf': dt.datetime(2008, 1, 1),
                    'wetprf': dt.datetime(2008, 1, 1),
                    'atmprf': dt.datetime(2008, 1, 1),
                    'scnlv1': dt.datetime(2008, 1, 1)}}
_test_download = {'': {kk: True for kk in tags.keys()}}
_password_req = {'': {kk: False for kk in tags.keys()}}


# ----------------------------------------------------------------------------
# Instrument methods

def init(self):
    """Initializes the Instrument object with instrument specific values.

    Runs once upon instantiation.

    """
    ack = ' '.join((''))
    refs = ' '.join(('Y. Liou et al., "FORMOSAT-3/COSMIC GPS',
                     'Radio Occultation Mission: Preliminary',
                     'Results," in IEEE Transactions on',
                     'Geoscience and Remote Sensing, vol. 45,',
                     'no. 11, pp. 3813-3826, Nov. 2007, doi:',
                     '10.1109/TGRS.2007.903365.\n',
                     'Additional information can be found at',
                     'https://cdaac-www.cosmic.ucar.edu/cdaac/doc/cosmic.html'))
    self.acknowledgements = ack
    self.references = refs
    logger.info(ack)

    return


def clean(self):
    """Return COSMIC GPS data cleaned to the specified level.

    Parameters
    ----------
    self : pysat.Instrument
        Instrument class object, whose attribute clean_level is used to return
        the desired level of data selectivity.

    Note
    ----
    Supports 'clean', 'dusty', 'dirty'

    """
    if self.tag == 'ionprf':
        # Ionosphere density profiles
        if self.clean_level == 'clean':
            # Filter out profiles where source provider processing doesn't work.
            self.data = self.data.where(self['edmaxalt'] != -999., drop=True)
            self.data = self.data.where(self['edmax'] != -999., drop=True)

            # Ensure 'edmaxalt' in "reasonable" range
            self.data = self.data.where(((self['edmaxalt'] >= 175.)
                                        & (self['edmaxalt'] <= 475.)),
                                        drop=True)

            # Filter densities when negative
            dens_copy = self['ELEC_dens'].values
            for i, profile in enumerate(self['time']):
                # Take out all densities below any altitude (< 325) with
                # a negative density.
                idx, = np.where((self[i, :, 'ELEC_dens'] < 0)
                                & (self[i, :, 'MSL_alt'] <= 325))
                if len(idx) > 0:
                    dens_copy[i, 0:idx[-1] + 1] = np.nan

                # Take out all densities above any altitude > 325 with a
                # negative density.
                idx, = np.where((self[i, :, 'ELEC_dens'] < 0)
                                & (self[i, :, 'MSL_alt'] > 325))
                if len(idx) > 0:
                    dens_copy[i, idx[0]:] = np.nan
            self[:, :, 'ELEC_dens'] = dens_copy

            # Do an altitude density gradient check to reduce number of
            # cycle slips.
            densDiff = self['ELEC_dens'].diff(dim='RO')
            altDiff = self['MSL_alt'].diff(dim='RO')
            normGrad = (densDiff / (altDiff * self[:, :-1, 'ELEC_dens']))

            # Calculate maximum gradient per profile
            normGrad = normGrad.max(dim='RO')

            # Remove profiles with high altitude gradients
            self.data = self.data.where(normGrad <= 1.)

    elif self.tag == 'scnlv1':
        # scintillation files
        if self.clean_level == 'clean':
            # Filter out profiles where source provider processing doesn't work.
            self.data = self.data.where(self['alttp_s4max'] != -999., drop=True)
            self.data = self.data.where(self['s4max9sec'] != -999., drop=True)

    return


# ----------------------------------------------------------------------------
# Instrument functions
#

# Set the list_files routine
def list_files(tag=None, inst_id=None, data_path=None, format_str=None):
    """Return a Pandas Series of every file for chosen satellite data.

    Parameters
    ----------
    tag : string or NoneType
        Denotes type of file to load.
        (default=None)
    inst_id : string or NoneType
        Specifies the satellite ID for a constellation.  Not used.
        (default=None)
    data_path : string or NoneType
        Path to data directory.  If None is specified, the value previously
        set in Instrument.files.data_path is used.  (default=None)
    format_str : NoneType
        User specified file format not supported here. (default=None)

    Return
    ------
    file_list : pysat.Files
        A class containing the verified available files

    """

    estr = 'Building a list of COSMIC files, which can possibly take time. '
    logger.info('{:s}~1s per 100K files'.format(estr))
    sys.stdout.flush()

    # Note that Files.from_os() could be used here except for the fact
    # that there are multiple COSMIC files per given time
    # here, we follow from_os() except a fictional microsecond
    # is added to file times to help ensure there are no file collisions

    # overloading revision keyword below
    if format_str is None:
        # COSMIC file format string
        format_str = ''.join(('*/*/*.{year:04d}.{day:03d}',
                              '.{hour:02d}.{minute:02d}.*_nc'))
    # process format string to get string to search for
    search_dict = futils.construct_searchstring_from_format(format_str)
    search_str = search_dict['search_string']
    # perform local file search
    files = futils.search_local_system_formatted_filename(data_path, search_str)
    # we have a list of files, now we need to extract the information
    # pull of data from the areas identified by format_str
    stored = futils.parse_delimited_filenames(files, format_str, delimiter='.')
    if len(stored['year']) > 0:
        year = np.array(stored['year'])
        day = np.array(stored['day'])
        hour = np.array(stored['hour'])
        minute = np.array(stored['minute'])
        try:
            uts = hour*3600.0 + minute*60.0
        except TypeError as err:
            raise TypeError(' '.join(('unable to construct time from',
                                      'filename\n{:}'.format(str(err)))))
        # do a pre-sort on uts to get files that may conflict with each other
        # due to multiple spacecraft and antennas
        # this ensures that we can make the times all unique for the file list
        idx = np.argsort(uts)
        # adding linearly increasing offsets less than 0.1 s
        shift_uts = np.mod(np.arange(len(year)), 9E4) * 1.E-5 + 1.E-5
        uts[idx] += shift_uts

        index = pysat.utils.time.create_datetime_index(year=year,
                                                       day=day,
                                                       uts=uts)
        if not index.is_unique:
            raise ValueError(' '.join(('Generated non-unique datetimes for',
                                       'COSMIC within list_files.')))
        # store sorted file names with unique times in index
        file_list = np.array(stored['files'])
        file_list = pds.Series(file_list, index=index)
        return file_list

    else:
        logger.info('Found no files, check your path or download them.')
        return pds.Series(None, dtype='object')


def load(fnames, tag=None, inst_id=None, altitude_bin=None):
    """Load COSMIC GPS files.

    Parameters
    ----------
    fnames : pandas.Series
        Series of filenames
    tag : str or NoneType
        tag or None (default=None)
    inst_id : str or NoneType
        satellite id or None (default=None)
    altitude_bin : integer
        Number of kilometers to bin altitude profiles by when loading.
        Currently only supported for tag='ionprf'.

    Returns
    -------
    output : pandas.DataFrame
        Object containing satellite data
    meta : pysat.Meta
        Object containing metadata such as column names and units

    """

    # Input check
    if altitude_bin is not None:
        if tag != 'ionprf':
            estr = 'altitude_bin keyword only supported for "tag=ionprf"'
            raise ValueError(estr)

    num = len(fnames)
    # Make sure there are files to read
    if num != 0:
        # Call generalized load_files routine
        output = load_files(fnames, tag=tag, inst_id=inst_id)

        # Create datetime index
        utsec = output.hour * 3600. + output.minute * 60. + output.second
        output['index'] = \
            pysat.utils.time.create_datetime_index(year=output.year.values,
                                                   month=output.month.values,
                                                   day=output.day.values,
                                                   uts=utsec.values)
        # Rename index
        output = output.rename(index='time')

        # Ensure time is increasing
        output = output.sortby('time')

        if tag == 'ionprf':
            # Set up coordinates
            coord_labels = ['MSL_alt', 'GEO_lat', 'GEO_lon', 'OCC_azi']
            var_labels = ['ELEC_dens', 'TEC_cal']

            # Apply coordinates to loaded data.
            output = output.set_coords(coord_labels)

            if altitude_bin is not None:
                # Deal with altitude binning, can't do it directly with
                # xarray since all dimensions get grouped.

                coord_labels.extend(['MSL_bin_alt'])
                all_labels = []
                all_labels.extend(coord_labels)
                all_labels.extend(var_labels)

                # Normalize and round actual altitude values by altitude_bin
                bin_alts = (output['MSL_alt'] / altitude_bin).round().values

                # Reconstruct altitude from bin_alts value
                alts = bin_alts * altitude_bin

                # Create array for bounds of each bin that data will be
                # grouped into.
                bin_arr = np.arange(np.nanmax(bin_alts))

                # Indexing information mapping which altitude goes to which bin
                dig_bins = np.digitize(bin_alts, bin_arr)

                # Create arrays to store results
                new_coords = {}
                for label in all_labels:
                    new_coords[label] = np.full(
                        (len(output['time']), len(bin_arr)), np.nan)

                # Go through each profile and mean values in each altitude bin.
                # Solution inspired by
                # (https://stackoverflow.com/questions/38013778/
                # is-there-any-numpy-group-by-function)
                # However, unique didn't work how I wanted for multi-dimensional
                # array, thus the for loop.
                for i in range(len(output['time'])):
                    ans = np.unique(dig_bins[i, :], return_index=True)

                    for label in all_labels:
                        if label == 'MSL_bin_alt':
                            temp_calc = np.split(alts[i, :], ans[1][1:])
                        else:
                            temp_calc = np.split(output[label].values[i, :],
                                                 ans[1][1:])
                        new_coords[label][i, :] = [np.mean(temp_vals) for
                                                   temp_vals in temp_calc]

                # Create new Dataset with binned data values.
                # First, prep coordinate data.
                coords = {}
                data_vars = {}
                for key in coord_labels:
                    coords[key] = (('time', 'RO'), new_coords[key])
                coords['time'] = output['time']

                # Create data_vars input dict
                for key in var_labels:
                    data_vars[key] = (('time', 'RO'), new_coords[key])

                # Create new Dataset
                new_set = xr.Dataset(data_vars=data_vars, coords=coords)

                # Copy over other variables
                for key in output.data_vars:
                    if key not in all_labels:
                        new_set[key] = output[key]

                # Replace initial Dataset
                output = new_set

        # Use the first available file to pick out meta information
        meta = pysat.Meta()
        ind = 0
        repeat = True
        while repeat:
            try:
                data = netCDF4.Dataset(fnames[ind])
                ncattrsList = data.ncattrs()
                for d in ncattrsList:
                    meta[d] = {meta.labels.units: '',
                               meta.labels.name: d}
                keys = data.variables.keys()
                for key in keys:
                    if 'units' in data.variables[key].ncattrs():
                        meta[key] = {
                            meta.labels.units: data.variables[key].units,
                            meta.labels.name: data.variables[key].long_name}
                repeat = False
            except RuntimeError:
                # file was empty, try the next one by incrementing ind
                ind += 1

        return output, meta
    else:
        # no data
        return xr.Dataset(None), pysat.Meta()


# separate routine for doing actual loading. This was broken off from main load
# because I was playing around with multiprocessor loading
# yielded about 20% improvement in execution time
def load_files(files, tag=None, inst_id=None, altitude_bin=None):
    """Load COSMIC data files directly from a given list.

    May be directly called by user, but in general is called by load.  This is
    separate from the main load function for future support of multiprocessor
    loading.

    Parameters
    ----------
    files : pandas.Series
        Series of filenames
    tag : str or NoneType
        tag or None (default=None)
    inst_id : str or NoneType
        satellite id or None (default=None)
    altitude_bin : integer
        Number of kilometers to bin altitude profiles by when loading.
        Currently only supported for tag='ionprf'.

    Returns
    -------
    output : list of dicts, one per file
        Object containing satellite data

    """
    output = [None] * len(files)
    drop_idx = []

    # Dict to store information about each data variable and data lengths
    # from each file loaded.
    main_dict = {}
    main_dict_len = {}

    # List of all data variables in the file
    data_var_keys = []

    # Iterate through files and load data
    for (i, fname) in enumerate(files):
        try:
            # Open file for access
            data = netCDF4.Dataset(fname)

            # Get list of file attributes, which includes information about
            # where the profile is observed, and store.
            ncattrsList = data.ncattrs()
            file_attrs = {}
            for d in ncattrsList:
                file_attrs[d] = data.getncattr(d)

            # Get a list of all data variables from the first file only
            if i == 0:
                for key in data.variables.keys():
                    data_var_keys.append(key)
                    main_dict[key] = []
                    main_dict_len[key] = []

            # Load all of the variables in the netCDF
            for key in data_var_keys:
                # Grab data
                t_list = data.variables[key][:]

                # Reverse byte order if needed and store
                if t_list.dtype.byteorder != '=':
                    main_dict[key].append(t_list.byteswap().newbyteorder())
                else:
                    main_dict[key].append(t_list)

                # Store length of data for the file
                main_dict_len[key].append(len(main_dict[key][-1]))

            output[i] = file_attrs
            data.close()

        except RuntimeError:
            # Some of the files have zero bytes, which causes a read error.
            # Store the index of these zero byte files so they can be dropped.
            drop_idx.append(i)

    # Drop anything that came from the zero byte files
    drop_idx.reverse()
    for i in drop_idx:
        del output[i]

    # Each GPS occultation has a different number of data points.
    # Generate numpy arrays based upon the largest size.
    for key in main_dict_len.keys():
        main_dict_len[key] = np.max(main_dict_len[key])

    for key in main_dict.keys():
        data_arr = np.full((len(main_dict[key]), main_dict_len[key]), np.nan)
        for i in range(len(main_dict[key])):
            data_arr[i, 0:len(main_dict[key][i])] = main_dict[key][i]

        main_dict[key] = data_arr

    # Collect all simple variable output into a Dataset
    output = pds.DataFrame(output).to_xarray()
    for key in main_dict:
        output[key] = (['index', 'RO'], main_dict[key])

    return output

    # if tag == 'atmprf':
    #     # this file has three groups of variable lengths
    #     # each goes into its own DataFrame
    #     # two are processed here, last is processed like other
    #     # file types
    #     # see code just after this if block for more
    #     # general explanation on lines just below
    #     p_keys = ['OL_vec2', 'OL_vec1', 'OL_vec3', 'OL_vec4']
    #     p_dict = {}
    #     # get indices needed to parse data
    #     plengths = main_dict_len['OL_vec1']
    #     max_p_length = np.max(plengths)
    #     plengths, plengths2 = _process_lengths(plengths)
    #     # collect data
    #     for key in p_keys:
    #         p_dict[key] = main_dict.pop(key)
    #         _ = main_dict_len.pop(key)
    #     psub_frame = pds.DataFrame(p_dict)
    #
    #     # change in variables in this file type
    #     # depending upon the processing applied at UCAR
    #     if 'ies' in main_dict.keys():
    #         q_keys = ['OL_ipar', 'OL_par', 'ies', 'hes', 'wes']
    #     else:
    #         q_keys = ['OL_ipar', 'OL_par']
    #     q_dict = {}
    #     # get indices needed to parse data
    #     qlengths = main_dict_len['OL_par']
    #     max_q_length = np.max(qlengths)
    #     qlengths, qlengths2 = _process_lengths(qlengths)
    #     # collect data
    #     for key in q_keys:
    #         q_dict[key] = main_dict.pop(key)
    #         _ = main_dict_len.pop(key)
    #     qsub_frame = pds.DataFrame(q_dict)
    #
    #     max_length = np.max([max_p_length, max_q_length])
    #     length_arr = np.arange(max_length)
    #     # small sub DataFrames
    #     for i in np.arange(len(output)):
    #         output[i]['OL_vecs'] = psub_frame.iloc[plengths[i]:plengths[i+1], :]
    #         output[i]['OL_vecs'].index = \
    #             length_arr[:plengths2[i+1]-plengths2[i]]
    #         output[i]['OL_pars'] = qsub_frame.iloc[qlengths[i]:qlengths[i+1], :]
    #         output[i]['OL_pars'].index = \
    #             length_arr[:qlengths2[i+1]-qlengths2[i]]
    #
    # if tag == 'ionprf':
    #     if altitude_bin is not None:
    #         for out in output:
    #             rval = (out['profiles']['MSL_alt']/altitude_bin).round().values
    #             out['profiles'].index = rval * altitude_bin
    #             out['profiles'] = \
    #                 out['profiles'].groupby(out['profiles'].index.values).mean()
    #     else:
    #         for out in output:
    #             out['profiles'].index = out['profiles']['MSL_alt']
    #
    # return output


def download(date_array, tag, inst_id, data_path=None,
             user=None, password=None):
    """Download COSMIC GPS data.

    Parameters
    ----------
    date_array : array-like
        list of datetimes to download data for. The sequence of dates need not
        be contiguous.
    tag : string
        Tag identifier used for particular dataset. This input is provided by
        pysat. (default='')
    inst_id : string
        Instrument ID string identifier used for particular dataset. This input
        is provided by pysat. (default='')
    data_path : string
        Path to directory to download data to. (default=None)
    user : string or NoneType
        User string input used for download. Provided by user and passed via
        pysat. If an account is required for downloads this routine here must
        error if user not supplied. (default=None)
    password : string or NoneType
        Password for data download. (default=None)

    Note
    ----
    This routine is invoked by pysat and is not intended for direct use by
    the end user.

    """

    if tag == 'ionprf':
        sub_dir = 'ionPrf'
    elif tag == 'sonprf':
        sub_dir = 'sonPrf'
    elif tag == 'wetprf':
        sub_dir = 'wetPrf'
    elif tag == 'atmprf':
        sub_dir = 'atmPrf'
    elif tag == 'scnlv1':
        sub_dir = 'scnLv1'
    else:
        raise ValueError('Unknown cosmic_gps tag')

    for date in date_array:
        logger.info('Downloading COSMIC data for ' + date.strftime('%D'))
        yr, doy = pysat.utils.time.getyrdoy(date)
        yrdoystr = '{year:04d}/{doy:03d}'.format(year=yr, doy=doy)

        # Try re-processed data (preferred)
        auth = requests.auth.HTTPBasicAuth(user, password)
        try:
            # Construct path string for online file
            dwnld = ''.join(("https://data.cosmic.ucar.edu/gnss-ro/cosmic1",
                             "/repro2013/level2/", yrdoystr, "/", sub_dir,
                             '_repro2013',
                             '_{year:04d}_{doy:03d}.tar.gz'.format(year=yr,
                                                                   doy=doy)))

            # Make online connection
            req = requests.get(dwnld, auth=auth)
            req.raise_for_status()
        except requests.exceptions.HTTPError:
            # If response is negative, try post-processed data
            try:
                # Construct path string for online file
                dwnld = ''.join(("https://data.cosmic.ucar.edu/gnss-ro/cosmic1",
                                 "/postProc/level2/", yrdoystr, "/", sub_dir,
                                 '_postProc',
                                 '_{year:04d}_{doy:03d}.tar.gz'))
                dwnld = dwnld.format(year=yr, doy=doy)

                # Make online connection
                req = requests.get(dwnld, auth=auth)
                req.raise_for_status()
            except requests.exceptions.HTTPError as err:
                estr = ''.join((str(err), '\n', 'Data not found'))
                logger.info(estr)

        # Copy request info to tarball file with generated name in `fname`.
        fname = os.path.join(data_path,
                             ''.join(('cosmic_', sub_dir,
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
