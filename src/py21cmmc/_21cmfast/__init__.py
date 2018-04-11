"""
Provides thin wrappers around C functions.
"""
import os
import numpy as np
import glob
from ctypes import cdll
from ctypes import c_float, c_int, c_uint, c_double, c_char, c_char_p, Structure, POINTER, c_ulonglong
from numpy.ctypeslib import ndpointer

LOCATION = os.path.dirname(os.path.abspath(__file__))
SharedLibraryLocation = glob.glob("%s/drive_21cmMC_streamlined*.so" % LOCATION)[0]
lib21CMFAST = cdll.LoadLibrary(SharedLibraryLocation)

_drive_21cmMC = lib21CMFAST.drive_21CMMC


class BaseStructure(Structure):
    def __init__(self, **kwargs):
        """
        Ctypes.Structure with integrated default values.

        :param kwargs: values different to defaults
        :type kwargs: dict
        """

        values = type(self)._defaults_.copy()
        for (key, val) in kwargs.items():
            values[key] = val

        super().__init__(**values)  # Python 3 syntax


class ReturnParams(Structure):
    """
    Descriptive parameters of the return data from drive_21cmMC.
    """
    _fields_ = [
        ("HII_DIM", c_int ),
        ("BOX_LEN", c_float),
        ("nbins_ps", c_int ),
        ("ncells_los", c_int ),
        ("n_redshifts", c_int ),
        ("n_ps", c_int ),
        ("lightcone", c_int ),
    ]


class DriveReturn(Structure):
    _fields_ = [
        ('Sampled_z', POINTER(c_float)),
        ('GlobalNF', POINTER(c_float)),
        ('GlobalTb', POINTER(c_float)),
        ('PS_k', POINTER(c_float)),
        ('PS', POINTER(c_float)),
        ('LCBox', POINTER(c_float)),
        ('Parameters', ReturnParams)
    ]


class FlagOptions(BaseStructure):
    """
    Structure defining flag-style options (with defaults) to :func:`drive_21cmMC`.

    Parameters
    ----------
    redshifts : list or float
        If a list, the redshifts at which to evaluate co-eval boxes. If a scalar, implies that a light cone is required,
        and its value corresponds to TS_calc_z, i.e. the lowest redshift of the light-cone.

    INCLUDE_ZETA_PL : bool, optional
        Should always be zero (have yet to include this option)

    READ_FROM_FILE : bool, optional
            Whether to read/write Walker file data to file, or pass directly between Python/C.

    GenerateNewICs : bool, optional
        Generate on the fly ICs (Only required for 21CMMC). If this flag is not set, it uses the seed in init.c

    SUBCELL_RSDS : bool, optional
        Add sub-cell RSDs (currently doesn't work if Ts is not used)

    IONISATION_FCOLL_TABLE : bool, optional
        An interpolation table for ionisation collapsed fraction (will be removing this at some point)

    USE_FCOLL_TABLE : bool, optional
        An interpolation table for ionisation collapsed fraction (for Ts.c; will be removing this at some point)

    PERFORM_TS_CALC : bool, optional
        Whether to compute the full spin temperature evolution

    INHOMO_RECO : bool, optional
        Whether to perform inhomogeneous recombinations

    OUTPUT_GLOBAL_AVE : bool, optional
        Output the global averaged data

    PRINT_DATA_FILES : bool, optional
        Output the neutral fraction and 21cm PS text files

    PRINT_COEVAL_BOXES : bool, optional
        Output the co-eval brightness temperature datacubes

    PRINT_LIGHTCONE_BOXES : bool, optional
        Output the light-cone brightness temperature datacubes

    """

    def __init__(self, redshifts, **kwargs):

        # Perform some logic on the inputs to get good defaults.
        if hasattr(redshifts, "__len__"):
            number_redshifts = len(redshifts)
            lightcone_flag = False
            TS_calc_z = 0  # not required
            redshifts = (c_double*number_redshifts)(*redshifts)
        else:
            number_redshifts = 0
            lightcone_flag = True
            TS_calc_z = 1.0 * redshifts  # mult by 1 to copy the value
            redshifts = (c_double*0)()

        super().__init__(
            USE_LIGHTCONE=lightcone_flag, REDSHIFT=TS_calc_z,
            N_USER_REDSHIFT=number_redshifts,redshifts=redshifts,
            **kwargs
        )

    _fields_ = [
        ("N_USER_REDSHIFT", c_int),
        ("USE_LIGHTCONE", c_int),
        ("INCLUDE_ZETA_PL", c_int),
        ("REDSHIFT", c_float),
        ("READ_FROM_FILE", c_int),
        ("GenerateNewICs", c_int),
        ("SUBCELL_RSD", c_int),
        ("USE_FCOLL_IONISATION_TABLE", c_int),
        ("SHORTEN_FCOLL", c_int),
        ("USE_TS_FLUCT", c_int),
        ("INHOMO_RECO", c_int),
        ("STORE_DATA", c_int),
        ("PRINT_FILES", c_int),
        ("PRINT_COEVAL_21cmBoxes", c_int),
        ("PRINT_LIGHTCONE_21cmBoxes", c_int),
        ("redshifts", POINTER(c_double)),
    ]

    _defaults_ = {
        "INCLUDE_ZETA_PL": False,
        "READ_FROM_FILE": False,
        "GenerateNewICs": True,
        "SUBCELL_RSD": False,
        "USE_FCOLL_IONISATION_TABLE": False,
        "SHORTEN_FCOLL": False,
        "USE_TS_FLUCT": False,
        "INHOMO_RECO": False,
        "STORE_DATA": False,
        "PRINT_FILES": False,
        "PRINT_COEVAL_21cmBoxes": False,
        "PRINT_LIGHTCONE_21cmBoxes": False
    }


# TODO: give each of the parameters a description!
class AstroParams(BaseStructure):
    """
    Ctypes structure containing astrophysical parameters (with defaults) for :func:`drive_21cmMC`.

    Parameters
    ----------
    INHOMO_RECO : float, optional
    EFF_FACTOR_PL_INDEX : float, optional
    HII_EFF_FACTOR : float, optional
    R_BUBBLE_MAX : float, optional
    ION_Tvir_MIN : float, optional
    L_X : float, optional
    NU_X_THRESH : float, optional
    NU_X_BAND_MAX : float, optional
    NU_X_MAX : float, optional
    X_RAY_SPEC_INDEX : float, optional
    X_RAY_Tvir_MIN : float, optional
    X_RAY_Tvir_LOWERBOUND : float, optional
    X_RAY_Tvir_UPPERBOUND : float, optional
    F_STAR : float, optional
    t_STAR : float, optional
    N_RSD_STEPS : float, optional
    LOS_direction : int < 3, optional
    Z_HEAT_MAX : float, optional
        Maximum redshift used in the Ts.c computation. Typically fixed at z = 35,
        but can be set to lower if the user wants light-cones in the saturated
        spin temperature limit (T_S >> T_CMB)

    ZPRIME_STEP_FACTOR : float, optional
        Used to control the redshift step-size for the Ts.c integrals. Typically fixed at 1.02,
        can consider increasing to make the code faster if the user wants
        light-cones in the saturated spin temperature limit (T_S >> T_CMB)
    """

    def __init__(self, INHOMO_RECO, **kwargs):
        super().__init__(**kwargs)

        if self.R_BUBBLE_MAX<0:
            self.R_BUBBLE_MAX = 15.0
            if INHOMO_RECO == 1:
                self.R_BUBBLE_MAX = 50.0

        if self.X_RAY_Tvir_MIN < 0:
            self.X_RAY_Tvir_MIN = self.ION_Tvir_MIN

    _fields_ = [
        ("EFF_FACTOR_PL_INDEX", c_float),
        ("HII_EFF_FACTOR", c_float),
        ("R_BUBBLE_MAX", c_float),
        ("ION_Tvir_MIN", c_float),
        ("L_X", c_float),
        ("NU_X_THRESH", c_float),
        ("NU_X_BAND_MAX", c_float),
        ("NU_X_MAX", c_float),
        ("X_RAY_SPEC_INDEX", c_float),
        ("X_RAY_Tvir_MIN", c_float),
        ("X_RAY_Tvir_LOWERBOUND", c_float),
        ("X_RAY_Tvir_UPPERBOUND", c_float),
        ("F_STAR", c_float),
        ("t_STAR", c_float),
        ("N_RSD_STEPS", c_int),
        ("LOS_direction", c_int),
        ("Z_HEAT_MAX", c_float),
        ("ZPRIME_STEP_FACTOR", c_float)
    ]

    _defaults_ = {
        "EFF_FACTOR_PL_INDEX": 0.0,
        "HII_EFF_FACTOR": 30.0,
        "R_BUBBLE_MAX": -1.0,
        "ION_Tvir_MIN": 4.69897,
        "L_X": 40.0,
        "NU_X_THRESH": 500.0,
        "NU_X_BAND_MAX": 2000.0,
        "NU_X_MAX": 10000.0,
        "X_RAY_SPEC_INDEX": 1.0,
        "X_RAY_Tvir_MIN": -1.0,
        "X_RAY_Tvir_LOWERBOUND": 4.0,
        "X_RAY_Tvir_UPPERBOUND": 6.0,
        "F_STAR": 0.05,
        "t_STAR": 0.5,
        "N_RSD_STEPS": 20,
        "LOS_direction": 2,
        "Z_HEAT_MAX": 35.0,
        "ZPRIME_STEP_FACTOR": 1.02,
    }


class CosmoParams(BaseStructure):
    """
    Ctypes Structure with cosmological parameters (with defaults) for :func:`drive_21cmMC`.

    Parameters
    ----------
    RANDOM_SEED : float, optional
        A seed to set the IC generator. If None, chosen from uniform distribution.

    SIGMA_8 : float, optional
        RMS mass variance (power spectrum normalisation).

    hlittle : float, optional
        H_0/100.

    OMEGA_M : float, optional
        Omega matter.

    OMEGA_b : float, optional
        Omega baryon, the baryon component.

    NS : float, optional
        Spectral index of the power spectrum.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.OMEGA_L = 1 - self.OMEGA_M
        if not self.RANDOM_SEED:
            self.RANDOM_SEED = int(np.random.randint(1, 1e12))

    _fields_ = [
        ("RANDOM_SEED", c_ulonglong),
        ("SIGMA8", c_float),
        ("hlittle", c_float),
        ("OMEGA_M", c_float),
        ("OMEGA_L", c_float),
        ("OMEGA_b", c_float),
        ("NS", c_float),
    ]

    _defaults_ = {
        "RANDOM_SEED": 12345,
        "SIGMA8": 0.82,
        'hlittle': 0.68,
        'OMEGA_M': 0.31,
        "OMEGA_b": 0.048,
        "NS": 0.97
    }


class BoxDim(BaseStructure):
    """
    Structure containin box size parameters (with defaults) for drive_21cmMC.

    Parameters
    ----------
    HII_DIM : int, optional
        Number of cells for the low-res box.

    DIM : int,optional
        Number of cells for the high-res box (sampling ICs) along a principal axis. To avoid
        sampling issues, DIM should be at least 3 or 4 times HII_DIM, and an integer multiple.
        By default, it is set to 4*HII_DIM.

    BOX_LEN : float, optional
        Length of the box, in Mpc.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.DIM < 0:
            self.DIM = 4*self.HII_DIM

    _fields_ = [
        ("HII_DIM", c_int),
        ("DIM", c_int),
        ("BOX_LEN", c_float)
    ]

    _defaults_ = {
        "BOX_LEN": 150.0,
        "DIM": -1,
        "HII_DIM": 100
    }

# TODO: this could be cleaned up just by requiring Point to contain structures itself.
class LightCone:
    def __init__(self, raw_return):
        """
        A data class representing the returned lightcone information from the MCMC driver.

        Parameters
        ----------
        raw_return : Point class
            The raw returned structure.
        """
        self.params = raw_return.Parameters
        self.is_lightcone = self.params.lightcone

        self.redshifts = np.ctypeslib.as_array(raw_return.Sampled_z, shape=(self.params.n_redshifts,))
        self.average_nf = np.ctypeslib.as_array(raw_return.GlobalNF, shape=(self.params.n_redshifts,))
        self.average_Tb = np.ctypeslib.as_array(raw_return.GlobalTb, shape=(self.params.n_redshifts,))

        # for i in range(self.n_redshifts):
        #     self.redshifts[i] = raw_return.Sampled_z[i]
        #     self.average_nf[i] = raw_return.GlobalNF[i]
        #     self.average_Tb[i] = raw_return.GlobalTb[i]

        if self.is_lightcone:
            self.redshifts = self.redshifts[::-1]
            self.average_nf = self.average_nf[::-1]
            self.average_Tb = self.average_Tb[::-1]

            k_offset = 2
        else:
            k_offset = 1

        self.power_spectrum = np.zeros((self.params.n_ps, self.params.nbins_ps - k_offset))
        self.k = np.zeros(self.params.nbins_ps - k_offset)

        for i in range(self.params.n_ps):
            for j in range(self.params.nbins_ps - k_offset):
                self.power_spectrum[i][j] = (raw_return.PS[(j + k_offset) + self.params.nbins_ps * i])
                if i == 0:
                    self.k[j] = (raw_return.PS_k[(j + k_offset) + self.params.nbins_ps * 0])

        self.lightcone_box = np.ctypeslib.as_array(
            raw_return.LCBox,
            shape=(self.params.HII_DIM, self.params.HII_DIM, self.params.ncells_los)
        )
        # # TODO: this must be much slower than just doing a numpy array
        # for k in range(self.ncells_los):
        #     for j in range(self.HII_DIM):
        #         for i in range(self.HII_DIM):
        #             self.lightcone_box[i, j, k] = raw_return.LCBox[k + self.ncells_los * (j + self.HII_DIM * i)]


boxdir = os.path.expanduser(os.path.join("~", ".py21cmmc", "Boxes"))

_drive_21cmMC.restype = DriveReturn
_drive_21cmMC.argtypes = [c_char_p, c_char_p, BoxDim, FlagOptions, AstroParams, CosmoParams]


def drive_21cmMC(random_ids, box_dim, flag_options, astro_params, cosmo_params):
    assert type(box_dim) == BoxDim
    assert type(flag_options) == FlagOptions
    assert type(astro_params) == AstroParams
    assert type(cosmo_params) == CosmoParams
    print(box_dim.DIM, box_dim.HII_DIM)
    # Force the types to adhere to Ctypes
    # flag_options = FlagOptions(**flag_options)
    # astro_params = AstroParams(**astro_params)
    # cosmo_params = CosmoParams(**cosmo_params)
    # print(cosmo_params.OMm)

    output = _drive_21cmMC(
        str.encode('%s' % (random_ids[0])),
        str.encode('%s' % (random_ids[1])),
        box_dim,
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
