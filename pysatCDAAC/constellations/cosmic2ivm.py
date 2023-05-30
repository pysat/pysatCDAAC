"""Create a constellation with all COSMIC2 IVM instruments

Attributes
----------
instruments : list
    List of pysat.Instrument objects

"""
import pysat
import pysatCDAAC

inst_ids = pysatCDAAC.instrument.cosmic2_ivm.inst_ids.keys()
instruments = [pysat.Instrument(inst_module=pysatCDAAC.instruments.cosmic2_ivm,
                                inst_id=inst_id) for inst_id in inst_ids]
