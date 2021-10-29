import numpy as np
import pytest
import tempfile

import pysat
# Import the test classes from pysat
from pysat.utils import generate_instrument_list
from pysat.tests.instrument_test_class import InstTestClass

# Make sure to import your instrument library here
import pysatCDAAC


# Retrieve the lists of CDAAC instruments and testing methods
instruments = generate_instrument_list(inst_loc=pysatCDAAC.instruments)

method_list = [func for func in dir(InstTestClass)
               if callable(getattr(InstTestClass, func))]

# Search tests for iteration via pytestmark, update instrument list
for method in method_list:
    if hasattr(getattr(InstTestClass, method), 'pytestmark'):
        # Get list of names of pytestmarks
        Nargs = len(getattr(InstTestClass, method).pytestmark)
        names = [getattr(InstTestClass, method).pytestmark[j].name
                 for j in range(0, Nargs)]
        # Add instruments from your library
        if 'all_inst' in names:
            mark = pytest.mark.parametrize("inst_name", instruments['names'])
            getattr(InstTestClass, method).pytestmark.append(mark)
        elif 'download' in names:
            mark = pytest.mark.parametrize("inst_dict", instruments['download'])
            getattr(InstTestClass, method).pytestmark.append(mark)
        elif 'no_download' in names:
            mark = pytest.mark.parametrize("inst_dict",
                                           instruments['no_download'])
            getattr(InstTestClass, method).pytestmark.append(mark)


class TestInstruments(InstTestClass):
    """Uses class level setup and teardown so that all tests use the same
    temporary directory. We do not want to geneate a new tempdir for each test,
    as the load tests need to be the same as the download tests.
    """

    def setup_class(self):
        """Runs once before the tests to initialize the testing setup."""
        # Make sure to use a temporary directory so that the user's setup is not
        # altered
        self.tempdir = tempfile.TemporaryDirectory()
        self.saved_path = pysat.params['data_dirs']
        pysat.params.data['data_dirs'] = [self.tempdir.name]

        # Assign the location of the Instrument sub-modules
        self.inst_loc = pysatCDAAC.instruments

    def teardown_class(self):
        """Runs after every method to clean up previous testing."""
        pysat.params.data['data_dirs'] = self.saved_path
        self.tempdir.cleanup()
        del self.inst_loc, self.saved_path, self.tempdir

    @pytest.mark.parametrize("inst_dict", [x for x in instruments['download']])
    def test_altitude_bin_keyword(self, inst_dict):
        """Test altitude bin keyword."""

        if inst_dict['tag'] in ['scnlv1', 'podtec', 'ionphs']:
            pytest.skip("Binning not available for scnlv1 data")
            return

        self.test_inst = pysat.Instrument(inst_module=inst_dict['inst_module'],
                                          tag=inst_dict['tag'],
                                          inst_id=inst_dict['inst_id'],
                                          altitude_bin=5.)
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

        return
