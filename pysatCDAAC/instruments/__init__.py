__all__ = ['cosmic_gps']

for inst in __all__:
    exec("from pysatCDAAC.instruments import {x}".format(x=inst))
