Installation
============

The following instructions will allow you to install pysatCDAAC.

Prerequisites
-------------

.. image:: figures/poweredbypysat.png
    :width: 150px
    :align: right
    :alt: powered by pysat Logo, blue planet with orbiting python


pysatCDAAC uses common Python modules, as well as modules developed by
and for the Space Physics community.  This module officially supports
Python 3.9+ and pysat 3.2.1+.

 ================== ====================
 Common modules     Community modules
 ================== ====================
  netCDF4            pysat>=3.2.1
  numpy
  pandas
  requests
  xarray
 ================== ====================

Installation Options
--------------------

1. Clone the git repository
::


   git clone https://github.com/pysat/pysatCDAAC.git


2. Install pysatCDAAC:
   Change directories into the repository folder and build the project.
   There are a few ways you can do this:

   A. Install on the system (root privileges required)::


        sudo pip install .

   B. Install at the user level::


        pip install --user .

   C. Install with the intent to change the code::


        pip install --user -e .


..
    .. extras-require:: test
        :pyproject:

..
    .. extras-require:: doc
        :pyproject:

.. _post-install:

Post Installation
-----------------

After installation, you may register the :py:mod:`pysatCDAAC`
:py:class:`Instrument` sub-modules with pysat.  If this is your first time using
pysat, check out the `quickstart guide
<https://pysat.readthedocs.io/en/latest/quickstart.html>`_ for pysat. Once pysat
is set up, you may choose to register the the :py:mod:`pysatCDAAC`
:py:class:`Instruments` sub-modules by:

.. code:: python


   import pysat
   import pysatCDAAC

   pysat.utils.registry.register_by_module(pysatCDAAC.instruments)

You may then use the pysat :py:attr:`platform` and :py:attr:`name` keywords to
initialize the model :py:class:`Instrument` instead of the
:py:attr:`inst_module` keyword argument.
