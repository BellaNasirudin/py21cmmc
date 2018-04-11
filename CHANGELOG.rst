
Changelog
=========

feature/update_to_brads_version
-------------------------------

Things modified from Brad's March21st version:
* All C files copied over verbatim.
* Added USE_LIGHTCONE to ReturnData
* Added thin wrapper around the driver in _21cmfast/__init__.py, including class structure which unpacks the return.
* Added more transparent/automatic handling of inputs in _utils.py.
* Added input parameter structs for easier setting of arguments.
* Added CLI calls to do exactly what Brad does in his script (minus running init.c)
