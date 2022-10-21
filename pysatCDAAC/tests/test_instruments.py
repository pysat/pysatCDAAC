"""Unit and Integration Tests for each instrument module.

Note
----
Imports test methods from pysat.tests.instrument_test_class

"""

import numpy as np
import pytest

import pysat
# Import the test classes from pysat
from pysat.tests.classes import cls_instrument_library as clslib

# Make sure to import your instrument library here
import pysatCDAAC


# Retrieve the lists of CDAAC instruments and testing methods
instruments = clslib.InstLibTests.initialize_test_package(
    clslib.InstLibTests, inst_loc=pysatCDAAC.instruments)


class TestInstruments(clslib.InstLibTests):
    """Main class for instrument tests.

    Note
    ----
    All standard tests, setup, and teardown inherited from the core pysat
    instrument test class.

    """

    @pytest.mark.second
    @pytest.mark.parametrize("inst_dict", [x for x in instruments['download']])
    @pytest.mark.parametrize("bin_num", [100, 200])
    def test_altitude_bin_keyword(self, inst_dict, bin_num):
        """Test altitude binning keywords.

        Parameters
        ----------
        inst_dict : dict
            Dictionary of instrument properties generated by pysat.
        bin_num : int
            Number of bins for altitude profiling.

        """

        if inst_dict['tag'] in ['scnlv1', 'podtec', 'ionphs']:
            pytest.skip("Binning not available for level-1 data")
            return

        self.test_inst = pysat.Instrument(inst_module=inst_dict['inst_module'],
                                          tag=inst_dict['tag'],
                                          inst_id=inst_dict['inst_id'],
                                          altitude_bin=5.,
                                          altitude_bin_num=bin_num)
        date = inst_dict['inst_module']._test_dates[inst_dict['inst_id']]
        date = date[inst_dict['tag']]
        self.test_inst.load(date=date)

        # Confirm presence of binned altitudes.
        assert 'MSL_bin_alt' in self.test_inst.data

        # Confirm binned altitudes are even factors of binning.
        rem = np.remainder(self.test_inst['MSL_bin_alt'].values, 5)
        idx, idy, = np.where(rem == 0)
        idx2, idy2, = np.where(np.isnan(rem))
        assert len(idx) + len(idx2) == np.product(rem.shape)

        # Confirm length of each profile corresponds to bin_num
        assert self.test_inst.data.dims['RO'] == bin_num

        return
