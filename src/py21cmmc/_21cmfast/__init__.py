"""
Provides thin wrappers around C functions.
"""
import os
import numpy as np
import glob
from ctypes import cdll
from ctypes import c_float, c_int, c_uint, c_double, c_char, c_char_p, Structure, POINTER
from numpy.ctypeslib import ndpointer



LOCATION = os.path.dirname(os.path.abspath(__file__))
SharedLibraryLocation = glob.glob("%s/drive_21cmMC_streamlined*.so" % LOCATION)[0]
lib21CMFAST = cdll.LoadLibrary(SharedLibraryLocation)

_drive_21cmMC = lib21CMFAST.drive_21CMMC


class Point(Structure):
    _fields_ = [
        ('Sampled_z', POINTER(c_float)),
        ('GlobalNF', POINTER(c_float)),
        ('GlobalTb', POINTER(c_float)),
        ('PS_k', POINTER(c_float)),
        ('PS', POINTER(c_float)),
        ('LCBox', POINTER(c_float)),
        ('Variables', POINTER(c_float))
    ]

#TODO: this could be cleaned up just by requiring Point to contain structures itself.
class LightCone:
    def __init__(self, raw_return):
        """
        A data class representing the returned lightcone information from the MCMC driver.

        Parameters
        ----------
        raw_return : Point class
            The raw returned structure.
        """
        self.HII_DIM = int(raw_return.Variables[0])
        self.box_length = raw_return.Variables[1]
        self.nbins_ps = int(raw_return.Variables[2])
        self.ncells_los = int(raw_return.Variables[3])
        self.n_redshifts = int(raw_return.Variables[4])
        self.n_ps = int(raw_return.Variables[5])
        self.lightcone = bool(raw_return.Variables[6])

        self.redshifts = np.zeros(self.n_redshifts)
        self.average_nf = np.zeros(self.n_redshifts)
        self.average_Tb = np.zeros(self.n_redshifts)

        for i in range(self.n_redshifts):
            self.redshifts[i] = raw_return.Sampled_z[i]
            self.average_nf[i] = raw_return.GlobalNF[i]
            self.average_Tb[i] = raw_return.GlobalTb[i]

        if self.lightcone:
            self.redshifts = self.redshifts[::-1]
            self.average_nf = self.average_nf[::-1]
            self.average_Tb = self.average_Tb[::-1]

            k_offset = 2
        else:
            k_offset = 1

        self.power_spectrum = np.zeros((self.n_ps, self.nbins_ps - k_offset))
        self.k = np.zeros(self.nbins_ps - k_offset)

        for i in range(self.n_ps):
            for j in range(self.nbins_ps - k_offset):
                self.power_spectrum[i][j] = (raw_return.PS[(j + k_offset) + self.nbins_ps * i])
            if i==0:
                self.k[j] = (raw_return.PS_k[(j + k_offset) + self.nbins_ps * 0])

boxdir = os.path.expanduser(os.path.join("~",".py21cmmc","Boxes"))

_drive_21cmMC.restype = Point
_drive_21cmMC.argtypes = [c_char_p,c_char_p, POINTER(c_float), POINTER(c_float), POINTER(c_float)]


def drive_21cmMC(random_ids, flag_options, astro_params, cosmo_params):
    output = _drive_21cmMC(
        '%s' % (random_ids[0]),
        '%s' % (random_ids[1]),
        flag_options,
        astro_params,
        cosmo_params
    )

    return LightCone(output)


# def run_init(**init_params):
#     if init_params:
#         set_globals(**init_params)
#
#     print('-' * 80)
#     print("\t PERFORMING INIT.C")
#     print('-'*80)
#     _run_init(c_char_p(boxdir.encode("utf-8")))
#     print("-- DONE "+ '-'*72)
#     print()
#
#
# def run_perturb(redshift, **init_params):
#     if init_params:
#         set_globals(**init_params)
#
#     print('-' * 80)
#     print("\t PERFORMING PERTURB_FIELD.C")
#     print("-"*80)
#     _run_perturb(c_char_p(boxdir.encode("utf-8")), redshift)
#     print("-- DONE "+ '-'*72)
#
# # ====== Functions for getting/setting globals =========================================================================
# c_set_globals.argtypes = [
#     c_float, # BOX_LEN
#     c_int,   # DIM
#     c_int,   # HII_DIM
#     c_int,   # NUM_CORES
#     c_float, # RAM
# ]
#
#
# class InitGlobals(Structure):
#     _fields_ = [
#         ('BOX_LEN', c_float),
#         ("DIM", c_int),
#         ("HII_DIM", c_int),
#         ("NUMCORES", c_int),
#         ("RAM", c_float)
#     ]
#
#
# get_globals.restype = InitGlobals
#
#
# def set_globals(**kwargs):
#     globs = get_globals()
#     args = []
#     for fld in InitGlobals._fields_:
#         name = fld[0]
#         if name in kwargs:
#             args.append(kwargs.pop(name))
#         else:
#             args.append(getattr(globs,name))
#
#     if len(kwargs)>0:
#         print("WARNING: Some of the parameters passed are incorrect: %s. Possible arguments are %s"%(kwargs, [x[0] for x in Globals._fields_]))
#
#     c_set_globals(*args)
# # ======================================================================================================================
#
# # ======= get_box_parameters wrapper ===================================================================================
# class BoxParameters(Structure):
#     _fields_ = [
#         ("box_len", c_float),
#         ("N", c_int),
#         ("vel_comp", c_int),
#     ]
#
#
# get_box_parameters.restype = BoxParameters
# # ======================================================================================================================
#
#
# # ======= drive_21cmMC wrapper =========================================================================================
# class DriveResult(Structure):
#     _fields_ = [
#         ('ave', c_float),
#         ('global_xH', c_float)
#     ]
#
#
# drive_21cmMC.restype = DriveResult
# drive_21cmMC.argtypes = [
#     ndpointer(np.complex64, ndim=1,flags='aligned, contiguous'),  # deltax_unfiltered
#     ndpointer(np.float32, ndim=1,flags='aligned, contiguous'),    # v
#     c_float,                                                      # REDSHIFT
#     c_float,                                                      # ION_EFF_FACTOR
#     c_float,                                                      # MFP
#     c_float,                                                      # TVIR_MIN
#     c_int,                                                        # PERFORM_PS
#     ndpointer(np.float32, ndim=1, flags='aligned, contiguous')    # delta_T (input empty).
# ]
# # ======================================================================================================================
#
#
# # ======= c_generatePS wrapper =========================================================================================
# class PowerSpectrumResult(Structure):
#     _fields_ = [
#         ("k_box", ndpointer(np.float64, ndim=1, flags='aligned, contiguous')),
#         ("p_box", ndpointer(np.float64, ndim=1, flags='aligned, contiguous')),
#         ("in_bin_ct", ndpointer(np.uint64, ndim=1, flags='aligned, contiguous')),
#         ("nbins", c_int)
#     ]
#
# c_generatePS.restype = PowerSpectrumResult
# c_generatePS.argtypes = [
#     c_float,  # ave
#     c_float,  # global_xH
#     ndpointer(np.float32, ndim=1, flags='aligned, contiguous')   # delta_T
# ]
# # ======================================================================================================================
#
#
# # ======= c_likelihood_chisquare wrapper ===============================================================================
# c_likelihood_chisquare.restype = c_double
#
# c_likelihood_chisquare.argtypes = [
#     PowerSpectrumResult,
#     ndpointer(np.float64, ndim=1, flags='aligned, contiguous'),  # mock
#     ndpointer(np.float64, ndim=1, flags='aligned, contiguous'),  # sensitivity
#     c_float,  # ForegroundCut
#     c_float,  # ShotNoiseCut
#     c_float,  # ModellingUncertainty
# ]
# # ======================================================================================================================



