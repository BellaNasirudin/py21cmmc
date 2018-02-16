"""
Provides thin wrappers around C functions.
"""
import os
import numpy as np
import glob
from ctypes import cdll
from ctypes import c_float, c_int, c_uint, c_double, Structure
from numpy.ctypeslib import ndpointer


__all__ = ['get_box_parameters', "drive_21cmMC", 'c_generatePS', 'c_likelihood_chisquare']

LOCATION = os.path.dirname(os.path.abspath(__file__))
SharedLibraryLocation = glob.glob("%s/drive_21cmMC_streamlined*.so" % LOCATION)[0]
lib21CMFAST = cdll.LoadLibrary(SharedLibraryLocation)

get_box_parameters = lib21CMFAST.GetBoxParameters
drive_21cmMC = lib21CMFAST.drive_21CMMC
c_generatePS = lib21CMFAST.generatePS
c_likelihood_chisquare = lib21CMFAST.likelihood_chi_square


# ======= get_box_parameters wrapper ===================================================================================
class BoxParameters(Structure):
    _fields_ = [
        ("box_len", c_float),
        ("N", c_int),
        ("vel_comp", c_int),
    ]


get_box_parameters.restype = BoxParameters
# ======================================================================================================================


# ======= drive_21cmMC wrapper =========================================================================================
class DriveResult(Structure):
    _fields_ = [
        ('ave', c_float),
        ('global_xH', c_float)
    ]


drive_21cmMC.restype = DriveResult
drive_21cmMC.argtypes = [
    ndpointer(np.complex64, ndim=1,flags='aligned, contiguous'),  # deltax_unfiltered
    ndpointer(np.float32, ndim=1,flags='aligned, contiguous'),    # v
    c_float,                                                      # REDSHIFT
    c_float,                                                      # ION_EFF_FACTOR
    c_float,                                                      # MFP
    c_float,                                                      # TVIR_MIN
    c_int,                                                        # PERFORM_PS
    ndpointer(np.float32, ndim=1, flags='aligned, contiguous')    # delta_T (input empty).
]
# ======================================================================================================================


# ======= c_generatePS wrapper =========================================================================================
class PowerSpectrumResult(Structure):
    _fields_ = [
        ("k_box", ndpointer(np.float64, ndim=1, flags='aligned, contiguous')),
        ("p_box", ndpointer(np.float64, ndim=1, flags='aligned, contiguous')),
        ("in_bin_ct", ndpointer(np.uint64, ndim=1, flags='aligned, contiguous')),
        ("nbins", c_int)
    ]

c_generatePS.restype = PowerSpectrumResult
c_generatePS.argtypes = [
    c_float,  # ave
    c_float,  # global_xH
    ndpointer(np.float32, ndim=1, flags='aligned, contiguous')   # delta_T
]
# ======================================================================================================================


# ======= c_likelihood_chisquare wrapper ===============================================================================
c_likelihood_chisquare.restype = c_double

c_likelihood_chisquare.argtypes = [
    PowerSpectrumResult,
    ndpointer(np.float64, ndim=1, flags='aligned, contiguous'),  # mock
    ndpointer(np.float64, ndim=1, flags='aligned, contiguous'),  # sensitivity
    c_float,  # ForegroundCut
    c_float,  # ShotNoiseCut
    c_float,  # ModellingUncertainty
]
# ======================================================================================================================



