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
- 'wetprf': 'wetPrf' files
- 'atmprf': 'atmPrf' files
- 'eraprf': 'eraPrf' files
- 'gfsprf': 'gfsPrf' files

Warnings
--------
- Routine was not produced by COSMIC team

"""

import datetime as dt
import os
import requests
import tarfile

import numpy as np
import netCDF4
import pandas as pds
import pysat
from pysat import logger
from pysat.utils import files as futils


# ----------------------------------------------------------------------------
# Instrument attributes

platform = 'cosmic'
name = 'gps'

l1_tags = ['ionPhs', 'podTec', 'scnLv1']
lower_l1_tags = [tag.lower() for tag in l1_tags]

l2_tags = ['ionPrf', 'wetPrf', 'atmPrf', 'eraPrf', 'gfsPrf']
lower_l2_tags = [tag.lower() for tag in l2_tags]

tags = {'ionprf': 'Ionospheric Profiles',
        'wetprf': '',
        'atmprf': 'Atmospheric Profiles',
        'eraprf': '',
        'gfsprf': ''}
for tag in l1_tags:
    tags[tag] = ''

inst_ids = {'': lower_l2_tags,
            'level_1b': lower_l1_tags}

# ----------------------------------------------------------------------------
# Instrument test attributes

_test_dates = {'': {'ionprf': dt.datetime(2008, 1, 1),
                    'wetprf': dt.datetime(2008, 1, 1),
                    'atmprf': dt.datetime(2008, 1, 1),
                    'eraprf': dt.datetime(2008, 1, 1),
                    'gfsprf': dt.datetime(2008, 1, 1)},
               'level_1b': {'ionphs': dt.datetime(2008, 1, 1),
                            'podtec': dt.datetime(2008, 1, 1),
                            'scnlv1': dt.datetime(2008, 1, 1)}
               }
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
        # ionosphere density profiles
        if self.clean_level == 'clean':
            # try and make sure all data is good
            # filter out profiles where source provider processing doesn't
            # get max dens and max dens alt
            self.data = self.data[((self['edmaxalt'] != -999.)
                                   & (self['edmax'] != -999.))]
            # make sure edmaxalt in "reasonable" range
            self.data = self.data[((self['edmaxalt'] >= 175.)
                                   & (self['edmaxalt'] <= 475.))]
            # filter densities when negative
            for i, profile in enumerate(self['profiles']):
                # take out all densities below the highest altitude negative
                # dens below 325
                idx, = np.where((profile.ELEC_dens < 0)
                                & (profile.index <= 325))
                if len(idx) > 0:
                    profile.iloc[0:(idx[-1] + 1)] = np.nan
                # take out all densities above the lowest altitude negative
                # dens above 325
                idx, = np.where((profile.ELEC_dens < 0)
                                & (profile.index > 325))
                if len(idx) > 0:
                    profile.iloc[idx[0]:] = np.nan

                # do an altitude density gradient check to reduce number of
                # cycle slips
                densDiff = profile.ELEC_dens.diff()
                altDiff = profile.MSL_alt.diff()
                normGrad = (densDiff / (altDiff * profile.ELEC_dens)).abs()
                idx, = np.where((normGrad > 1.) & normGrad.notnull())
                if len(idx) > 0:
                    self[i, 'edmaxalt'] = np.nan
                    self[i, 'edmax'] = np.nan
                    self[i, 'edmaxlat'] = np.nan
                    profile['ELEC_dens'] *= np.nan

        # filter out any measurements where things have been set to NaN
        self.data = self.data[self['edmaxalt'].notnull()]

    elif self.tag == 'scnlv1':
        # scintillation files
        if self.clean_level == 'clean':
            # try and make sure all data is good
            # filter out profiles where source provider processing doesn't
            # work
            self.data = self.data[((self['alttp_s4max'] != -999.)
                                   & (self['s4max9sec'] != -999.))]

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

    # Note that Files.from_os() could be used here except for the fact
    # that there are multiple COSMIC files per given time.
    # Instead, we follow from_os() except a fictional amount of time
    # is added to file times based upon the satellite and ground station
    # numbers to ensure there are no file collisions.

    # Overloading revision and cycle keyword below
    if format_str is None:
        # COSMIC file format string
        if inst_id == '' or (tag == 'ionphs'):
            format_str = ''.join(('*/*/*_C{revision:03d}.{year:04d}.',
                                  '{day:03d}.{hour:02d}.{minute:02d}.',
                                  'G{cycle:02d}_{version:04d}.????_nc'))
        elif inst_id == 'level_1b':
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
    global lower_l1_tags

    # input check
    if altitude_bin is not None:
        if tag != 'ionprf':
            estr = 'altitude_bin keyword only supported for "tag=ionprf"'
            raise ValueError(estr)

    num = len(fnames)
    # make sure there are files to read
    if num != 0:
        # call separate load_files routine, segmented for possible
        # multiprocessor load, not included and only benefits about 20%
        output = pds.DataFrame(load_files(fnames, tag=tag, inst_id=inst_id,
                                          altitude_bin=altitude_bin))
        utsec = output.hour * 3600. + output.minute * 60. + output.second
        # FIXME: need to switch to xarray so unique time stamps not needed
        # make times unique by adding a unique amount of time less than a second
        if tag not in lower_l1_tags or (tag == 'ionphs'):
            # add 1E-6 seconds to time based upon occulting_inst_id
            # additional 1E-7 seconds added based upon cosmic ID
            # get cosmic satellite ID
            c_id = np.array([snip[3] for snip in output.fileStamp]).astype(int)
            # time offset
            if tag != 'ionphs':
                utsec += output.occulting_sat_id * 1.e-5 + c_id * 1.e-6
            else:
                utsec += output.occsatId * 1.e-5 + c_id * 1.e-6
        else:
            # construct time out of three different parameters
            # duration must be less than 10,000
            # prn_id is allowed two characters
            # antenna_id gets one
            # prn_id and antenna_id are not sufficient for a unique time
            utsec += output.prn_id * 1.e-2 + output.duration.astype(int) * 1.E-6
            utsec += output.antenna_id * 1.E-7
        # move to Index
        output.index = \
            pysat.utils.time.create_datetime_index(year=output.year,
                                                   month=output.month,
                                                   day=output.day,
                                                   uts=utsec)
        if not output.index.is_unique:
            raise ValueError('Datetimes returned by load_files not unique.')
        # make sure UTS strictly increasing
        output.sort_index(inplace=True)
        # use the first available file to pick out meta information
        profile_meta = pysat.Meta()
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
                        profile_meta[key] = {
                            meta.labels.units: data.variables[key].units,
                            meta.labels.name: data.variables[key].long_name}
                repeat = False
            except RuntimeError:
                # file was empty, try the next one by incrementing ind
                ind += 1

        meta['profiles'] = profile_meta
        return output, meta
    else:
        # no data
        return pds.DataFrame(None), pysat.Meta()


def _process_lengths(lengths):
    """Prep lengths for parsing.

    Internal func used by load_files.
    """

    lengths = lengths.tolist()
    lengths.insert(0, 0)
    lengths = np.array(lengths)
    lengths2 = lengths.copy()
    lengths[-1] += 1
    return lengths, lengths2


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
    main_dict = {}
    main_dict_len = {}

    safe_keys = []
    for (i, fname) in enumerate(files):
        try:
            data = netCDF4.Dataset(fname)
            # build up dictionary will all ncattrs
            new = {}
            # get list of file attributes
            ncattrsList = data.ncattrs()
            # these include information about where the profile observed
            for d in ncattrsList:
                new[d] = data.getncattr(d)

            if i == 0:
                keys = data.variables.keys()
                for key in keys:
                    safe_keys.append(key)
                    main_dict[key] = []
                    main_dict_len[key] = []

            # load all of the variables in the netCDF
            for key in safe_keys:
                # grab data
                t_list = data.variables[key][:]
                # reverse byte order if needed
                if t_list.dtype.byteorder != '=':
                    main_dict[key].append(t_list.byteswap().newbyteorder())
                else:
                    main_dict[key].append(t_list)
                # store lengths
                main_dict_len[key].append(len(main_dict[key][-1]))

            output[i] = new
            data.close()
        except RuntimeError:
            # some of the files have zero bytes, which causes a read error
            # this stores the index of these zero byte files so I can drop
            # the Nones the gappy file leaves behind
            drop_idx.append(i)

    # drop anything that came from the zero byte files
    drop_idx.reverse()
    for i in drop_idx:
        del output[i]

    # combine different sub lists in main_dict into one
    for key in safe_keys:
        main_dict[key] = np.hstack(main_dict[key])
        main_dict_len[key] = np.cumsum(main_dict_len[key])

    if tag == 'atmprf':
        # this file has three groups of variable lengths
        # each goes into its own DataFrame
        # two are processed here, last is processed like other
        # file types
        # see code just after this if block for more
        # general explanation on lines just below
        p_keys = ['OL_vec2', 'OL_vec1', 'OL_vec3', 'OL_vec4']
        p_dict = {}
        # get indices needed to parse data
        plengths = main_dict_len['OL_vec1']
        max_p_length = np.max(plengths)
        plengths, plengths2 = _process_lengths(plengths)
        # collect data
        for key in p_keys:
            p_dict[key] = main_dict.pop(key)
            _ = main_dict_len.pop(key)
        psub_frame = pds.DataFrame(p_dict)

        # change in variables in this file type
        # depending upon the processing applied at UCAR
        if 'ies' in main_dict.keys():
            q_keys = ['OL_ipar', 'OL_par', 'ies', 'hes', 'wes']
        else:
            q_keys = ['OL_ipar', 'OL_par']
        q_dict = {}
        # get indices needed to parse data
        qlengths = main_dict_len['OL_par']
        max_q_length = np.max(qlengths)
        qlengths, qlengths2 = _process_lengths(qlengths)
        # collect data
        for key in q_keys:
            q_dict[key] = main_dict.pop(key)
            _ = main_dict_len.pop(key)
        qsub_frame = pds.DataFrame(q_dict)

        max_length = np.max([max_p_length, max_q_length])
        length_arr = np.arange(max_length)
        # small sub DataFrames
        for i in np.arange(len(output)):
            output[i]['OL_vecs'] = \
                psub_frame.iloc[plengths[i]:plengths[i + 1], :]
            output[i]['OL_vecs'].index = \
                length_arr[:plengths2[i + 1] - plengths2[i]]
            output[i]['OL_pars'] = \
                qsub_frame.iloc[qlengths[i]:qlengths[i + 1], :]
            output[i]['OL_pars'].index = \
                length_arr[:qlengths2[i + 1] - qlengths2[i]]

    # create a single data frame with all bits, then
    # break into smaller frames using views
    main_frame = pds.DataFrame(main_dict)
    # get indices needed to parse data
    lengths = main_dict_len[list(main_dict.keys())[0]]
    # get largest length and create numpy array with it
    # used to speed up reindexing below
    max_length = np.max(lengths)
    length_arr = np.arange(max_length)
    # process lengths for ease of parsing
    lengths, lengths2 = _process_lengths(lengths)
    # break main profile data into each individual profile
    for i in np.arange(len(output)):
        output[i]['profiles'] = main_frame.iloc[lengths[i]:lengths[i + 1], :]
        output[i]['profiles'].index = length_arr[:lengths2[i + 1] - lengths2[i]]

    if tag == 'ionprf':
        if altitude_bin is not None:
            for out in output:
                rval = (out['profiles']['MSL_alt']
                        / altitude_bin).round().values
                out['profiles'].index = rval * altitude_bin
                out['profiles'] = \
                    out['profiles'].groupby(out['profiles'].index.values).mean()
        else:
            for out in output:
                out['profiles'].index = out['profiles']['MSL_alt']

    return output


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
    global l1_tags, lower_l1_tags, l2_tags, lower_l2_tags

    if inst_id == '':
        level_str = 'level2'
        matches = [utag for ltag, utag in zip(lower_l2_tags, l2_tags)
                   if tag == ltag]
        sub_str = matches[0]
    elif inst_id == 'level_1b':
        level_str = 'level1b'
        matches = [utag for ltag, utag in zip(lower_l1_tags, l1_tags)
                   if tag == ltag]
        sub_str = matches[0]

    for date in date_array:
        logger.info('Downloading COSMIC data for ' + date.strftime('%D'))
        yr, doy = pysat.utils.time.getyrdoy(date)
        yrdoystr = '{year:04d}/{doy:03d}'.format(year=yr, doy=doy)

        # Try re-processed data (preferred)
        auth = requests.auth.HTTPBasicAuth(user, password)
        try:
            # Construct path string for online file
            dwnld = ''.join(("https://data.cosmic.ucar.edu/gnss-ro/cosmic1",
                             "/repro2013/", level_str, "/", yrdoystr, "/",
                             sub_str, '_repro2013',
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
                                 "/postProc/", level_str, "/", yrdoystr, "/",
                                 sub_str, '_postProc',
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
