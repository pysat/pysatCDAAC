Loading COSMIC2 IVM data
=====================

pysatCDAAC uses `pysat <https://github.com/pysat/pysat>`_ to load
space science instrument data.  As specified in the
`pysat tutorial <https://pysat.readthedocs.io/en/latest/tutorial.html>`_,
data may be loaded using the following commands.  Data from the Ion Velocity
Meter on board the Ionospheric CONnection Explorer `(ICON) <https://www.nasa.gov/icon>`_ is used as an example.

::


   import datetime as dt
   import pysat
   import pysatCDAAC as py_cdaac

   pysat.utils.registry.register_by_module(py_cdaac.instruments)

   date = dt.datetime(2021, 1, 1)
   ivm = pysat.Instrument(platform='cosmic2', name='ivm',
                          inst_id='e1', update_files=True)
   ivm.download(start=date)
   ivm.load(date=date)
   print(ivm)


The output shows some basic info about the instrument object, as well as
information about the data loaded.

::

  pysat Instrument object
  -----------------------
  Platform: 'cosmic2'
  Name: 'ivm'
  Tag: ''
  Instrument id: 'e1'

  Data Processing
  ---------------
  Cleaning Level: 'clean'
  Data Padding: None
  Custom Functions: 0 applied

  Local File Statistics
  ---------------------
  Number of files: 2
  Date Range: 01 January 2021 --- 02 January 2021

  Loaded Data Statistics
  ----------------------
  Date: 01 January 2021
  DOY: 001
  Time range: 01 January 2021 00:00:01 --- 02 January 2021 00:00:00
  Number of Times: 85276
  Number of variables: 46

  Variable Names:
  alt        ap_pot     ap_pot_var
                 ...
  rpa_flag   sc_flag    uts

  pysat Meta object
  -----------------
  Tracking 7 metadata values
  Metadata for 46 standard variables
  Metadata for 0 ND variables
  Metadata for 6 global attributes
