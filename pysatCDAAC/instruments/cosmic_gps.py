# -*- coding: utf-8 -*-
"""Loads data from the COSMIC satellite.

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
    Select profile type, or scintillation, using one of the following keys:
    {'ionprf': 'Ionospheric Profiles',
     'wetprf': 'Atmospheric profiles with moisture',
     'atmprf': 'Atmospheric profiles without moisture',
     'eraprf': 'ERA-40 Interim reanalysis data',
     'gfsprf': 'NCEP operational analysis data',
     'ionphs': 'Ionospheric excess phase',
     'podtec': 'Absolute Total Electron Content and auxiliary data',
     'scnlv1': 'S4 scintillation index and auxiliary data'}
inst_id
    None supported
altitude_bin
    Number of kilometers to bin altitude profiles by when loading.
    Works for all files except tag='scnlv1', 'podtec', or 'ionphs'.

Warnings
--------
- Routine was not produced by COSMIC team
- Files are labeled with times at minute resolution which can result in multiple
  COSMIC data profiles at the same time. pysat requires that instruments have
  monotonic and unique times, thus, to meet pysat requirements a time shift
  (based upon file/data parameters) is added to each profile to ensure
  all times are unique. This time shift within a minute
  is not considered significant given the released data structure.
  For level-1b data files time shifts are distributed throughout the minute,
  for level-2 files the time shifts are less than .0001 seconds.
  The difference in time distribution is related to the availability of a
  unique combination of parameters at the different file levels.

"""

import datetime as dt
import functools

import netCDF4
import numpy as np
import pandas as pds
import xarray as xr

import pysat
from pysat.utils import files as futils
from pysatCDAAC.instruments.methods import general as mm_cdaac


# ----------------------------------------------------------------------------
# Instrument attributes

platform = 'cosmic'
name = 'gps'

tags = {'ionprf': 'Ionospheric Profiles',
        'wetprf': 'Atmospheric profiles with moisture',
        'atmprf': 'Atmospheric profiles without moisture',
        'eraprf': 'ERA-40 Interim reanalysis data',
        'gfsprf': 'NCEP operational analysis data',
        'ionphs': 'Ionospheric excess phase',
        'podtec': 'Absolute Total Electron Content and auxiliary data',
        'scnlv1': 'S4 scintillation index and auxiliary data'}

inst_ids = {'': list(tags.keys())}

pandas_format = False

# Translate tags to levels and camelcase names
tag_translation = {'ionphs': {'level': 'level1b', 'substr': 'ionPhs'},
                   'podtec': {'level': 'level1b', 'substr': 'podTec'},
                   'scnlv1': {'level': 'level1b', 'substr': 'scnLv1'},
                   'ionprf': {'level': 'level2', 'substr': 'ionPrf'},
                   'wetprf': {'level': 'level2', 'substr': 'wetPrf'},
                   'atmprf': {'level': 'level2', 'substr': 'atmPrf'},
                   'eraprf': {'level': 'level2', 'substr': 'eraPrf'},
                   'gfsprf': {'level': 'level2', 'substr': 'gfsPrf'}}

coord_translation = {
    'ionprf': ['MSL_alt', 'GEO_lat', 'GEO_lon', 'OCC_azi'],
    'atmprf': ['MSL_alt', 'Lat', 'Lon', 'Azim'],
    'sonprf': ['MSL_alt', 'lat', 'lon'],
    'wetprf': ['MSL_alt', 'Lat', 'Lon'],
    'eraprf': ['MSL_alt', 'Lat', 'Lon', 'Pres', 'Temp', 'Vp', 'Ref'],
    'gfsprf': ['MSL_alt', 'Pres', 'Temp', 'Vp', 'Ref'],
    'ionphs': ['caL1Snr', 'pL1Snr', 'pL2Snr', 'xLeo', 'yLeo', 'zLeo', 'xdLeo',
               'ydLeo', 'zdLeo', 'xGps', 'yGps', 'zGps', 'xdGps', 'ydGps',
               'zdGps', 'exL1', 'exL2'],
    'podtec': ['x_GPS', 'y_GPS', 'z_GPS', 'x_LEO', 'y_LEO', 'z_LEO', 'TEC',
               'elevation', 'caL1_SNR', 'pL2_SNR', 'profile_time'],
    'scnlv1': ['alt_s4max', 'lat_s4max', 'lon_s4max', 'lct_s4max']}

# ----------------------------------------------------------------------------
# Instrument test attributes

_test_dates = {'': {}.fromkeys(list(tags.keys()), dt.datetime(2014, 5, 1))}

# ----------------------------------------------------------------------------
# Instrument methods


def init(self):
    """Initialize the Instrument object with instrument specific values.

    Note
    ----
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
    pysat.logger.info(ack)

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

            # Ensure 'edmaxalt' in "reasonable" range.
            self.data = self.data.where(((self['edmaxalt'] >= 175.)
                                        & (self['edmaxalt'] <= 475.)),
                                        drop=True)

            # Filter densities when negative.
            dens_copy = self['ELEC_dens'].values
            for i, profile in enumerate(self['time']):
                # Take out all densities below any altitude (< 325 km) with
                # a negative density.
                idx, = np.where((self[i, :, 'ELEC_dens'] < 0)
                                & (self[i, :, 'MSL_alt'] <= 325))
                if len(idx) > 0:
                    dens_copy[i, 0:idx[-1] + 1] = np.nan

                # Take out all densities above any altitude > 325 km with a
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

            # Calculate maximum gradient per profile.
            normGrad = normGrad.max(dim='RO')

            # Remove profiles with high altitude gradients
            self.data = self.data.where(normGrad <= 1.)

    elif self.tag == 'scnlv1':
        # Scintillation files
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
    tag : str or NoneType
        Denotes type of file to load.
        (default=None)
    inst_id : str or NoneType
        Specifies the satellite ID for a constellation.  Not used.
        (default=None)
    data_path : str or NoneType
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
    pysat.logger.info('{:s}~1s per 100K files'.format(estr))

    # Note that Files.from_os() could be used here except for the fact
    # that there are multiple COSMIC files per given time.
    # Instead, we follow from_os() except a fictional amount of time
    # is added to file times based upon the satellite and ground station
    # numbers to ensure there are no file collisions.

    # Overloading revision and cycle keyword below
    if format_str is None:
        # COSMIC file format string
        if tag_translation[tag]['level'] == 'level2' or (tag == 'ionphs'):
            format_str = ''.join(('*/*/*_C{revision:03d}.{year:04d}.',
                                  '{day:03d}.{hour:02d}.{minute:02d}.',
                                  'G{cycle:02d}_{version:04d}.????_nc'))
        elif tag_translation[tag]['level'] == 'level1b':
            format_str = ''.join(('*/*/*_C{revision:03d}.{year:04d}.',
                                  '{day:03d}.{hour:02d}.{minute:02d}.',
                                  '????.G{cycle:02d}.??_{version:04d}.????_nc'))

    # Process format string to get string to search for
    search_dict = futils.construct_searchstring_from_format(format_str)
    search_str = search_dict['search_string']

    # Perform local file search
    files = futils.search_local_system_formatted_filename(data_path, search_str)

    # We have a list of files, now we need to extract the information
    # from the areas identified by format_str.
    stored = futils.parse_fixed_width_filenames(files, format_str)

    # Process info
    if len(stored['year']) > 0:
        year = np.array(stored['year'])
        day = np.array(stored['day'])
        hour = np.array(stored['hour'])
        minute = np.array(stored['minute'])
        ver = np.array(stored['version'])

        # Satellite ID pulled out as revision
        rev = np.array(stored['revision'])

        # Ground Station pulled out as 'cycle'
        gs = np.array(stored['cycle'])

        # Create UTS time. Done with caution in case parsing above
        # wasn't correct for some reason.
        try:
            uts = hour * 3600.0 + minute * 60.0
        except TypeError as err:
            raise TypeError(' '.join(('unable to construct time from',
                                      'filename\n{:}'.format(str(err)))))

        # Add shift in time based upon ground station and satellite ID
        # to ensure files named by the minute are unique.
        uts += rev * 0.01 + gs * 0.0001

        index = pysat.utils.time.create_datetime_index(year=year,
                                                       day=day,
                                                       uts=uts)
        if not index.is_unique:
            # Look for duplicate times but different versions
            dups = index[index.duplicated()].unique()
            if len(dups) > 0:
                # Keep the highest version for duplicated times
                version = pds.Series(ver, index=index)
                frame = pds.DataFrame({'files': files, 'revive': version,
                                       'time': index}, index=index)
                frame = frame.sort_values(by=['time', 'revive'],
                                          ascending=[True, False])
                frame = frame.drop_duplicates(subset='time', keep='first')
                stored['files'] = frame['files'].values
                index = frame.index

            if not index.is_unique:
                raise ValueError(' '.join(('Generated non-unique datetimes for',
                                           'COSMIC within list_files.')))

        # Store file names with unique times in index
        file_list = np.array(stored['files'])
        file_list = pds.Series(file_list, index=index)
        return file_list
    else:
        return pds.Series(None, dtype='object')


def load(fnames, tag=None, inst_id=None, altitude_bin=None,
         altitude_bin_num=300):
    """Load COSMIC GPS files.

    Parameters
    ----------
    fnames : pandas.Series
        Series of filenames
    tag : str or NoneType
        tag or None (default=None)
    inst_id : str or NoneType
        satellite id or None (default=None)
    altitude_bin : int or NoneType
        Number of kilometers to bin altitude profiles by when loading.
        Works for all files except tag='scnlv1', 'podtec', or 'ionphs' as
        `MSL_alt` is required in the file. If None, no binning performed.
        (default=None)
    altitude_bin_num : int
        Number of bins to use when binning profile altitude if
        `altitude_bin` is not None. (default=300)

    Returns
    -------
    output : pandas.DataFrame
        Object containing satellite data
    meta : pysat.Meta
        Object containing metadata such as column names and units

    """

    # Input check.
    if altitude_bin is not None:
        if tag in ['podtec', 'scnlv1', 'ionphs']:
            estr = ' '.join(['altitude_bin keyword only supported if `MSL_alt`',
                             'present in the file.'])
            raise ValueError(estr)

    num = len(fnames)
    # Make sure there are files to read.
    if num != 0:
        coords = {}
        # Specify additional coords for datasets with nonstandard dimensions.
        if tag == 'atmprf':
            temp_keys = ['OL_vec2', 'OL_vec1', 'OL_vec3', 'OL_vec4']
            dim_label = 'dim1'
            for key in temp_keys:
                coords[key] = dim_label

            temp_keys = ['OL_ipar', 'OL_par', 'ies', 'hes', 'wes']
            dim_label = 'dim2'
            for key in temp_keys:
                coords[key] = dim_label

        # Call generalized load_files routine.
        output = load_files(fnames, tag=tag, inst_id=inst_id, coords=coords)

        # Create datetime index.
        utsec = output.hour * 3600. + output.minute * 60. + output.second

        # Not all profiles are unique within a minute sampling, thus
        # we add a small time offset to ensure unique times. A more consistent
        # offset time could be obtained by parsing the filenames as is done
        # in list files however load isn't passed `format_str`, thus this
        # solution wouldn't work in all cases.
        if tag_translation[tag]['level'] == 'level2' or (tag == 'ionphs'):
            # Add 1E-5 seconds to time based upon occulting_inst_id and an
            # additional 1E-6 seconds added based upon cosmic ID.
            # Get cosmic satellite ID.
            c_id = np.array([snip.values.tolist()[3]
                             for snip in output.fileStamp]).astype(int)
            # Time offset
            if tag != 'ionphs':
                utsec += output.occulting_sat_id * 1.e-5 + c_id * 1.e-6
            else:
                utsec += output.occsatId * 1.e-5 + c_id * 1.e-6
        else:
            # Construct time out of three different parameters:
            #   duration must be less than 100,000
            #   prn_id is allowed two characters
            #   antenna_id gets one character
            # prn_id and antenna_id alone are not sufficient for a unique time.
            if np.nanmax(output.duration) >= 1.e5:
                estr = ''.join(('Assumptions for the time shift calculation ',
                                'are not holding. Please contact pysatCDAAC ',
                                'developers.'))
                raise ValueError(estr)
            utsec += output.prn_id * 1.e-2 + output.duration.astype(int) * 1.E-6
            utsec += output.antenna_id * 1.E-7

        output['index'] = pysat.utils.time.create_datetime_index(
            year=output.year.values,
            month=output.month.values,
            day=output.day.values,
            uts=utsec.values)

        # Rename index to time.
        if tag_translation[tag]['level'] == 'level1b':
            # scnlv1 files already have a 2D time variable, it is a conflict.
            output = output.rename(time='profile_time')
        output = output.rename(index='time')

        # Ensure time is increasing.
        output = output.sortby('time')

        # Set up coordinates, depending upon file type.
        coord_labels = coord_translation[tag]

        # Apply coordinates to loaded data.
        output = output.set_coords(coord_labels)

        # Bin by altitude is requested by user
        if altitude_bin is not None and ('MSL_alt' in coord_labels):
            # Deal with altitude binning, can't do it directly with
            # xarray since all dimensions get grouped.

            # Technique depends upon altitude values being in sorted
            # ascending order. First, filter out negative altitudes.
            idx, idy, = np.where(output['MSL_alt'] < 0)
            output['MSL_alt'][idx, idy] = np.nan
            idx = output['MSL_alt'].argsort(axis=-1)

            # Get all variables/coordinates using same dimensions as MSL_alt
            # and sort.
            msl_dims = output['MSL_alt'].dims

            var_labels = [var for var in output.data_vars
                          if output[var].dims == msl_dims]
            coord_labels = [coord for coord in output.coords
                            if output[coord].dims == msl_dims]

            all_labels = []
            all_labels.extend(var_labels)
            all_labels.extend(coord_labels)

            for var in all_labels:
                if output[var].dims == msl_dims:
                    output[var] = (msl_dims,
                                   np.take_along_axis(output[var].values, idx,
                                                      axis=1))

            coord_labels.extend(['MSL_bin_alt'])
            all_labels.extend(['MSL_bin_alt'])

            # Normalize and floor actual altitude values by altitude_bin
            bin_alts = np.floor((output['MSL_alt'] / altitude_bin).values)

            # Reconstruct altitude from bin_alts value
            alts = bin_alts * altitude_bin

            # Create array for bounds of each bin that data will be
            # grouped into.
            bin_arr = np.arange(altitude_bin_num + 1)

            # Indexing information mapping which altitude goes to which bin
            dig_bins = np.digitize(bin_alts, bin_arr)

            # Create arrays to store results
            new_coords = {}
            for label in all_labels:
                new_coords[label] = np.full((len(output['time']),
                                             len(bin_arr) - 1), np.nan)

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
                    # Average all values in each bin. Guard against first
                    # realized bin being larger than first possible bin.
                    ir = dig_bins[i, 0] - 1
                    new_coords[label][i, ir:len(temp_calc) + ir] = [
                        np.mean(temp_vals)
                        for temp_vals in temp_calc][:len(bin_arr) - ir - 1]

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
                # File was empty, try the next one by incrementing ind
                ind += 1

        return output, meta
    else:
        # No data.
        return xr.Dataset(None), pysat.Meta()


# Separate routine for doing actual loading. This was broken off from main load
# because I was playing around with multiprocessor loading.
# Yielded about 20% improvement in execution time.
def load_files(files, tag=None, inst_id=None, coords=None):
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
    coords : dict or NoneType
        Dict keyed by data variable name that stores the coordinate name that
        should be assigned when loading data. If a variable name not provided
        will default to 'RO'. (default=None)

    Returns
    -------
    output : list of dicts, one per file
        Object containing satellite data

    """
    output = [None] * len(files)
    drop_idx = []

    if coords is None:
        coords = {}

    # Dict to store information about each data variable and data lengths
    # from each file loaded.
    main_dict = {}
    main_dict_len = {}

    # List of all data variables in the file.
    data_var_keys = []

    # Iterate through files and load data
    for (i, fname) in enumerate(files):
        try:
            # Open file for access.
            data = netCDF4.Dataset(fname)

            # Get list of file attributes, which includes information about
            # where the profile is observed, and store.
            ncattrsList = data.ncattrs()
            file_attrs = {}
            for d in ncattrsList:
                file_attrs[d] = data.getncattr(d)

            # Get a list of all data variables from the first file only.
            if i == 0:
                for key in data.variables.keys():
                    data_var_keys.append(key)
                    main_dict[key] = []
                    main_dict_len[key] = []

            # Load all of the variables in the netCDF.
            for key in data_var_keys:
                # Grab data.
                t_list = data.variables[key][:]

                # Reverse byte order if needed and store.
                if t_list.dtype.byteorder != '=':
                    main_dict[key].append(t_list.byteswap().newbyteorder())
                else:
                    main_dict[key].append(t_list)

                # Store length of data for the file.
                main_dict_len[key].append(len(main_dict[key][-1]))

            output[i] = file_attrs
            data.close()

        except RuntimeError:
            # Some of the files have zero bytes, which causes a read error.
            # Store the index of these zero byte files so they can be dropped.
            drop_idx.append(i)

    # Drop anything that came from the zero byte files.
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

    # Collect all simple variable output into a Dataset.
    output = pds.DataFrame(output).to_xarray()
    for key in main_dict:
        if key not in coords:
            coords[key] = 'RO'
        output[key] = (['index', coords[key]], main_dict[key])

    return output


# Use general CDAAC download routine.
download_tags = {
    '': {tag:
         {'remote_dir': ''.join(('gnss-ro/cosmic1/repro2013/',
                                 tag_translation[tag]['level'],
                                 '/{year:4d}/{day:03d}/')),
          'tar_name': ''.join((tag_translation[tag]['substr'],
                               '_repro2013_{year:4d}_{day:03d}.tar.gz')),
          'backup': ['repro2013', 'postProc']} for tag in tags.keys()}}
download = functools.partial(mm_cdaac.download, supported_tags=download_tags,
                             sub_path=True)
