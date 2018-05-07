"""
This module provides a set of higher-level wrappers to the underlying C functions of 21cmFAST/21CMMC.

Since the underlying C code is based on various global values being set in the driver (or sometimes in component
functions) and then filtering down into the component functions, to get a more "pythonic" interface, we attempt
to abstract that away from the user here, by providing functions that wrap each of the underlying functions, and
set the various globals within them, but expect *non-global* variables to be passed to them (with defaults etc.).
"""

from ._wrapped_21cmfast import ffi, lib
import numpy as np
from . import defaults
from ._utils import write_walker_cosmology_file, write_walker_file
from ._asarray import asarray
from os import path, mkdir
from tqdm import tqdm
from astropy.cosmology import Planck15, z_at_value
from astropy.units import Mpc


def _pyalloc(name, length, tp=None):
    "A function to allocate memory in python and set it to a C variable."
    v = getattr(lib, name)

    typemap = {
        'float *': np.float32,
        'double *': np.float64,
        'int *': np.int32
    }

    tp = tp or typemap[ffi.getctype(ffi.typeof(v))]

    if isinstance(length, int):
        x = np.zeros(length, dtype=tp)
    else:
        x = length
    setattr(lib, name, ffi.cast(ffi.typeof(v), ffi.from_buffer(x)))
    return x


class _StructWithDefaults:
    """
    A class which provides a convenient interface to create a C structure with defaults specified.

    Note: the actual C structure is accessed as the `cstruct` property of the class. This is not writeable, i.e. it is
    auto-generated when it is accessed, based on the parameters in the class. This provides a *fully initialised*
    structure, and will fail if not all fields are specified with defaults.

    It is provided for the purpose of *creating* C structures in Python to be passed to C functions, where sensible
    defaults are available. Structures which are created within C and passed back do not need to be wrapped.
    """
    _name = None

    def __init__(self, **kwargs):

        for k, v in self._defaults_.items():

            # Prefer arguments given to the constructor.
            if k in kwargs:
                v = kwargs[k]

            if isinstance(v, str):
                v = ffi.new('char[]', v.encode())

            try:
                setattr(self, k, v)
            except AttributeError:
                # The attribute has been defined as a property, save it as a hidden variable

                setattr(self, "_" + k, v)

        self._logic()

        # Set the name of this struct in the C code
        if self._name is None:
            self._name = self.__class__.__name__

    def _logic(self):
        pass

    @property
    def cstruct(self):
        "The C Structure corresponding to this instance."

        if hasattr(self, "_obj"):
            return self._obj

        self._obj = ffi.new("struct " + self._name + "*")
        self._logic()
        for fld in ffi.typeof(self._obj[0]).fields:
            try:
                setattr(self._obj, fld[0], getattr(self, fld[0]))
            except TypeError:
                print("For key %s, value %s:" % (fld[0], getattr(self, fld[0])))
                raise
        return self._obj

    @property
    def c_initializer(self):
        "A Python dictionary containing every field which needs to be initialized in the C struct."
        if hasattr(self, "_dct"):
            return self._dct

        obj = ffi.new("struct " + self._name + "*")
        self._logic()
        dct = {}
        for fld in ffi.typeof(obj[0]).fields:
            dct[fld[0]] = getattr(self, fld[0])
        return dct

    def set_globals(self):
        obj = ffi.new("struct " + self._name + "*")

        for k, v in self.c_initializer.items():
            try:
                setattr(lib, k, v)
            except TypeError:
                _pyalloc(k, v)


class FlagOptionStruct(_StructWithDefaults):
    """
    Structure defining flag-style options (with defaults) to :func:`drive_21cmMC`.

    Parameters
    ----------
    Z_HEAT_MAX : float
        The Astrophysical parameter, which helps to build the redshifts in the case that it's a lightcone.

    ZPRIME_STEP_FACTOR : float
        The Astrophysical parameter, which helps to build the redshifts in the case that it's a lightcone.

    redshifts : list or float, optional
        If a list, the redshifts at which to evaluate co-eval boxes. If a scalar (Default), implies that a light cone is required,
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
    _defaults_ = defaults.FlagOptions

    def __init__(self, Z_HEAT_MAX, ZPRIME_STEP_FACTOR, **kwargs):
        self.Z_HEAT_MAX = Z_HEAT_MAX
        self.ZPRIME_STEP_FACTOR = ZPRIME_STEP_FACTOR

        super().__init__(**kwargs)

    @property
    def redshifts(self):
        # First try returning cached value. Otherwise, new memory is assigned on each call!!
        if hasattr(self, "_redshifts_"):
            return self._redshifts_

        if self.USE_LIGHTCONE:
            z_prime = self.REDSHIFT * 1.0001
            while z_prime < self.Z_HEAT_MAX:
                z_prime = ((1. + z_prime) * self.ZPRIME_STEP_FACTOR - 1.)

            z_prime = ((1. + z_prime) / self.ZPRIME_STEP_FACTOR - 1.)
            redshifts = []
            while (z_prime > self.REDSHIFT):
                redshifts.append(z_prime)
                prev_z_prime = z_prime
                z_prime = ((1. + prev_z_prime) / self.ZPRIME_STEP_FACTOR - 1.)

            self._redshifts_ = np.array(redshifts)  # This keeps a persistent memory location for the data
            _pyalloc("redshifts", self._redshifts_)
            return self._redshifts_

        else:
            self._redshifts_ = np.array(self._redshifts)  # This keeps a persistent memory location for the data
            _pyalloc("redshifts", self._redshifts_)
            return self._redshifts_

    @property
    def USE_LIGHTCONE(self):
        return not hasattr(self._redshifts, "__len__")

    @property
    def N_USER_REDSHIFT(self):
        return len(self.redshifts)

    @property
    def REDSHIFT(self):
        return self._redshifts if self.USE_LIGHTCONE else 0

    def _logic(self):
        if self.GenerateNewICs and (self.USE_FCOLL_IONISATION_TABLE or self.SHORTEN_FCOLL):
            raise ValueError(
                """
                Cannot use interpolation tables when generating new initial conditions on the fly.
                (Interpolation tables are only valid for a single cosmology/initial condition)
                """
            )

        if self.USE_TS_FLUCT and self.INCLUDE_ZETA_PL:
            raise NotImplementedError(
                """
                Cannot use a non-constant ionising efficiency (zeta) in conjuction with the IGM spin temperature part of the code.
                (This will be changed in future)
                """
            )

        if self.USE_FCOLL_IONISATION_TABLE and self.INHOMO_RECO:
            raise ValueError(
                "Cannot use the f_coll interpolation table for find_hii_bubbles with inhomogeneous recombinations")

        if self.INHOMO_RECO and self.USE_TS_FLUCT:
            raise ValueError(
                """
                Inhomogeneous recombinations have been set, but the spin temperature is not being computed.
                "Inhomogeneous recombinations can only be used in combination with the spin temperature calculation (different from 21cmFAST).
                """
            )


class AstroParamStruct(_StructWithDefaults):
    """
    Ctypes structure containing astrophysical parameters (with defaults) for :func:`drive_21cmMC`.

    Parameters
    ----------
    INHOMO_RECO : bool, optional
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
    _defaults_ = defaults.AstroParams

    def __init__(self, INHOMO_RECO, **kwargs):
        self.INHOMO_RECO = INHOMO_RECO
        super().__init__(**kwargs)

    @property
    def R_BUBBLE_MAX(self):
        if not self._R_BUBBLE_MAX:
            return 50.0 if self.INHOMO_RECO else 15.0
        else:
            return self._R_BUBBLE_MAX

    @property
    def ION_Tvir_MIN(self):
        return 10 ** self._ION_Tvir_MIN

    @property
    def L_X(self):
        return 10 ** self._L_X

    @property
    def X_RAY_Tvir_MIN(self):
        return 10 ** self._X_RAY_Tvir_MIN if self._X_RAY_Tvir_MIN else self.ION_Tvir_MIN

    @property
    def NU_X_THRESH(self):
        return self._NU_X_THRESH * lib.NU_over_EV

    @property
    def NU_X_BAND_MAX(self):
        return self._NU_X_BAND_MAX * lib.NU_over_EV

    @property
    def NU_X_MAX(self):
        return self._NU_X_MAX * lib.NU_over_EV


class CosmoParamStruct(_StructWithDefaults):
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
    _defaults_ = defaults.CosmoParams

    @property
    def RANDOM_SEED(self):
        if not self._RANDOM_SEED:
            return int(np.random.randint(1, 1e12))
        else:
            return self._RANDOM_SEED

    @property
    def OMl(self):
        return 1 - self.OMm


class BoxDimStruct(_StructWithDefaults):
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
    _defaults_ = defaults.BoxDim

    @property
    def DIM(self):
        return self._DIM or 4 * self.HII_DIM


class LightCone:
    """
    This class takes numpy arrays, which are the raw input arrays defining the quantities, and just does some massaging
    to make them easier to handle in Python.
    """

    def __init__(self, redshifts, average_nf, average_Tb, power, k, lightcone_box, n_ps, num_bins, HII_DIM, box_len):

        self.redshifts_eval = redshifts[::-1]
        self.average_nf = average_nf[::-1]
        self.average_Tb = average_Tb[::-1]

        k_offset = 2
        self.power_spectrum = np.zeros((n_ps, num_bins - k_offset))
        self.k = np.zeros(num_bins - k_offset)


        for i in range(n_ps):
            for j in range(num_bins - k_offset):
                self.power_spectrum[i][j] = (power[(j + k_offset) + num_bins * i])
                if i == 0:
                    self.k[j] = (k[(j + k_offset) + num_bins * 0])

        self.lightcone_box = lightcone_box.reshape(
            (HII_DIM, HII_DIM, -1)
        )

        self.box_len = box_len
        self.box_len_lightcone = (self.box_len/HII_DIM)*self.lightcone_box.shape[-1]

        # Find the redshifts at each slice.
        init_d = Planck15.comoving_distance(self.redshifts_eval.min()).value
        final_d = init_d + self.box_len_lightcone
        d = np.linspace(init_d, final_d, self.lightcone_box.shape[-1]) * Mpc
        self.redshifts_slices = np.array([z_at_value(Planck15.comoving_distance, dd) for dd in d])

class CoEval:
    def __init__(self, redshifts, ave_nf, ave_Tb, delta_T, power_k=None, power_spectrum=None, ps_bins=None):
        """
        A data class representing the returned co-eval box information from the MCMC driver.
        """
        self.redshifts = redshifts
        self.n_redshifts = len(redshifts)

        self.average_nf = ave_nf
        self.average_Tb = ave_Tb

        k_offset = 1
        print(len(redshifts), ps_bins, k_offset)
        self.power_spectrum = np.zeros((self.n_redshifts, ps_bins - k_offset))
        self.k = np.zeros(ps_bins - k_offset)
        for i in range(len(redshifts)):
            for j in range(ps_bins - k_offset):
                self.power_spectrum[i][j] = (power_spectrum[(j + k_offset) + ps_bins * i])
                if i == 0:
                    self.k[j] = (power_k[(j + k_offset) + ps_bins * 0])

        self.delta_T = delta_T


# def drive_21cmmc(random_ids, box_dim, flag_options, astro_params, cosmo_params):
#     """
#     Wrapper of drive_21CMMC.
#
#     Returns
#     -------
#
#     """
#     box_dim = BoxDimStruct(**box_dim).c_initializer
#     flag_options = FlagOptionStruct(**flag_options).c_initializer
#     astro_params = AstroParamStruct(flag_options['INHOMO_RECO'], **astro_params).c_initializer
#     cosmo_params = CosmoParamStruct(**cosmo_params).c_initializer
#
#     if flag_options['READ_FROM_FILE']:
#         write_walker_file(flag_options, astro_params, random_ids)
#         write_walker_cosmology_file(flag_options, cosmo_params, random_ids)
#
#     out = lib.drive_21CMMC(
#         random_ids[0].encode(),
#         random_ids[1].encode(),
#         box_dim, flag_options, astro_params, cosmo_params
#     )
#
#     return LightCone(out)


def compute_initial_conditions(box_dim, cosmo_params, regenerate=False, write=True, dir='.'):
    box_dim.set_globals()
    cosmo_params.set_globals()
    # Initialize Power-Spectrum constants
    lib.init_ps()

    try:
        mkdir(path.join(dir, "Boxes"))
    except OSError:
        pass

    # The file where this data should be stored.
    # fl = path.join("Boxes", "init_boxes_%s_%.0fMpc.npz" % (box_dim.HII_DIM, box_dim.BOX_LEN)))

    suffix = "_{HII_DIM}_{BOX_LEN:.0f}Mpc".format(**box_dim.c_initializer)
    files = {
        'LOWRES_density': "smoothed_deltax_z0.00" + suffix,
        "LOWRES_vx": 'vxoverddot' + suffix,
        "LOWRES_vy": 'vyoverddot' + suffix,
        "LOWRES_vz": 'vzoverddot' + suffix,
        "LOWRES_vx_2LPT": 'vxoverddot_2LPT' + suffix,
        "LOWRES_vy_2LPT": 'vyoverddot_2LPT' + suffix,
        "LOWRES_vz_2LPT": 'vzoverddot_2LPT' + suffix,
        "HIRES_density": 'deltax_z0.00_{DIM}_{BOX_LEN}Mpc'.format(**box_dim.c_initializer)
    }

    if not regenerate:
        # Try to find the box
        # TODO: this would probably be better with a hash or something.
        try:
            init_boxes = {k: np.fromfile(path.join(dir, 'Boxes', v)) for k, v in files.items()}
            # give these to their C variables, so that the behaviour of this function is the same
            # whether or not the boxes already exist.
            for k, v in init_boxes.items():
                _pyalloc(k, v)

            return init_boxes
            # return np.load(fl)
        except FileNotFoundError:
            print("No init boxes found, computing...")

    # Also have to allocate some arrays which don't get allocated within the function
    hires_density = _pyalloc("HIRES_density", lib.TOT_FFT_NUM_PIXELS)
    lib.ComputeInitialConditions()

    dct = {k: asarray(ffi, getattr(lib, k), lib.HII_TOT_NUM_PIXELS)
           for k in ['LOWRES_density', "LOWRES_vx", "LOWRES_vy", "LOWRES_vz", "LOWRES_vx_2LPT",
                     "LOWRES_vy_2LPT", "LOWRES_vz_2LPT"]}
    dct["HIRES_density"] = hires_density

    if write:
        # Write each of the boxes to file, in the same way as 21cmFAST.
        for k, v in dct.items():
            v.tofile(path.join(dir, 'Boxes', files[k]))

    # Return some of the filled globals.
    # TODO: this is bad because the C program owns the data, so we can't easily delete it from Python.
    # Would be better to pass in arrays from Python.
    return dct


def compute_perturb_field(redshift, init_boxes, velocity_component=0, write=True, regenerate=False, read=False, dir='.'):
    """
    Note: the returned arrays are pointers to memory which C owns, and will overwrite if this function is called again.
          If the return values are used, they should be copied explicitly.

    You must have called compute_initial_conditions before calling this function -- it sets some of the globals in
    lib which are accessed here.

    """
    suffix = "_z{redshift:06.2f}_{HII_DIM}_{BOX_LEN:.0f}Mpc".format(redshift=redshift, HII_DIM=lib.HII_DIM,
                                                                    BOX_LEN=lib.BOX_LEN)
    files = {
        'LOWRES_density_REDSHIFT': "updated_smoothed_deltax" + suffix,
        'LOWRES_velocity_REDSHIFT': "updated_v{VEL_COMPONENT}".format(VEL_COMPONENT='xyz'[velocity_component]) + suffix
    }

    if not regenerate:
        # Try to find the boxes
        # TODO: this would probably be better with a hash or something.
        if all([path.exists(path.join(dir, 'Boxes', v)) for v in files.values()]):
            if read:
                return {k: np.fromfile(path.join(dir, 'Boxes', v)) for k, v in files.items()}
            else:
                return

    # Ensure we set the init_boxes into the global variables
    # If the init_boxes has been gotten from compute_initial_conditions, then it's already allocated, but
    # that should be ok (just resetting it to the same thing).
    init_boxes = dict(init_boxes)
    for k, v in init_boxes.items():
        _pyalloc(k, v)

    # Also have to allocate some arrays which don't get allocated within the function
    lowres_density = _pyalloc("LOWRES_density_REDSHIFT", lib.HII_TOT_NUM_PIXELS)
    lowres_vel = _pyalloc("LOWRES_velocity_REDSHIFT", lib.HII_TOT_NUM_PIXELS)

    lib.ComputePerturbField(redshift)
    dct = dict(
        LOWRES_density_REDSHIFT=lowres_density,
        LOWRES_velocity_REDSHIFT=lowres_vel
    )

    if write:
        # Write each of the boxes to file, in the same way as 21cmFAST.
        for k, v in dct.items():
            v.tofile(path.join(dir, 'Boxes', files[k]))
    return dct


def _setup_lightcone():
    lib.dR = (lib.BOX_LEN / lib.HII_DIM) * lib.CMperMPC
    lib.LC_BOX_PADDING = int(np.ceil(lib.LC_BOX_PADDING_IN_MPC / (lib.BOX_LEN * lib.HII_DIM)))
    lib.N_USER_REDSHIFT_LC = lib.N_USER_REDSHIFT - 1

    redshifts_lc = _pyalloc("redshifts_LC", lib.N_USER_REDSHIFT_LC)
    start_index_lc = _pyalloc("start_index_LC", lib.N_USER_REDSHIFT_LC)
    end_index_lc = _pyalloc("end_index_LC", lib.N_USER_REDSHIFT_LC)
    full_index_LC = _pyalloc("full_index_LC", 100000)
    slice_redshifts = _pyalloc('slice_redshifts', 100000)

    lib.z1_LC = lib.redshifts[lib.N_USER_REDSHIFT - 1]
    lib.z_LC = lib.z1_LC
    slice_ct = 0
    lib.total_slice_ct = 0
    lib.num_boxes_interp = 0

    i = 0
    while lib.z1_LC < lib.redshifts[0]:
        lib.z2_LC = lib.redshifts[lib.N_USER_REDSHIFT - 2 - i]
        while (lib.z_LC < lib.z2_LC):
            slice_redshifts[lib.total_slice_ct] = lib.z_LC
            full_index_LC[lib.total_slice_ct] = lib.total_slice_ct
            if slice_ct == lib.HII_DIM:
                lib.end_z = lib.z_LC
                lib.num_boxes_interp += 1
                slice_ct = 0

            slice_ct += 1
            lib.total_slice_ct += 1
            lib.z_LC -= lib.dR / lib.drdz(lib.z_LC)

        redshifts_lc[i] = lib.z1_LC
        if (i == 0):
            start_index_lc[i] = 0
            end_index_lc[i] = slice_ct
        else:
            start_index_lc[i] = end_index_lc[i - 1]
            end_index_lc[i] = slice_ct

        lib.z1_LC = lib.z2_LC
        i += 1

    lib.total_num_boxes = lib.num_boxes_interp
    lib.remainder_LC = lib.total_slice_ct - lib.num_boxes_interp * lib.HII_DIM;

    lib.final_z = lib.z_LC

    box_z1 = _pyalloc("box_z1", lib.HII_TOT_NUM_PIXELS)
    box_z2 = _pyalloc("box_z2", lib.HII_TOT_NUM_PIXELS)
    box_interpolate = _pyalloc("box_interpolate", lib.HII_TOT_NUM_PIXELS)
    box_interpolate_remainder = _pyalloc("box_interpolate_remainder", lib.HII_DIM * lib.HII_DIM * lib.remainder_LC)

    # Now, we need to remember to return every allocated memory space that still gets used elsewhere (just to keep it
    # alive)
    return (redshifts_lc, start_index_lc, end_index_lc, full_index_LC, slice_redshifts, box_z1, box_z2, box_interpolate,
            box_interpolate_remainder)


def _setup_erf():
    lib.erfc_arg_min = -15.0
    lib.erfc_arg_max = 15.0
    lib.ERFC_NUM_POINTS = 10000

    ERFC_VALS = _pyalloc("ERFC_VALS", lib.ERFC_NUM_POINTS)
    ERFC_VALS_DIFF = _pyalloc("ERFC_VALS_DIFF", lib.ERFC_NUM_POINTS)

    lib.ArgBinWidth = (lib.erfc_arg_max - lib.erfc_arg_min) / (lib.ERFC_NUM_POINTS - 1.)
    lib.InvArgBinWidth = 1. / lib.ArgBinWidth

    for i in range(lib.ERFC_NUM_POINTS):
        erfc_arg_val = lib.erfc_arg_min + lib.ArgBinWidth * i
        ERFC_VALS[i] = lib.splined_erfc(erfc_arg_val)

    ERFC_VALS_DIFF[:-1] = ERFC_VALS[1:] - ERFC_VALS[:-1]

    return (ERFC_VALS, ERFC_VALS_DIFF)


def compute_ionisation_boxes(flag_options={}, astro_params={}):
    """
    This *basically* runs  the whole of the driver, apart from the init stuff, and in the C code it frees all the
    arrays except for the "DataToBeReturned". This is *not* ideal, so we should change it. It means we need to use a
    few hacks here.
    """

    flag_options.set_globals()
    astro_params.set_globals()
    # ==================================================================================================================
    # Initialize some variables (found in the driver routine)
    # ==================================================================================================================
    if flag_options.GenerateNewICs:
        lowres_density_redshift = _pyalloc("LOWRES_density_REDSHIFT", lib.HII_TOT_NUM_PIXELS)
        lowres_velocity_redshift = _pyalloc("LOWRES_velocity_REDSHIFT", lib.HII_TOT_NUM_PIXELS)

    if flag_options.INHOMO_RECO:
        lib.INHOMO_RECO_R_BUBBLE_MAX = 50.0
        lib.R_BUBBLE_MAX = lib.INHOMO_RECO_R_BUBBLE_MAX

    if flag_options.USE_LIGHTCONE:
        # This conveniently removes most of the lightcone-specific logic. Several arrays are created in Python, but
        # set to C variables. These must be kept alive for the memory to not be de-allocated. These are kept alive
        # in the tmp_mem_buff object, which we shouldn't have to actually use. They will be killed at the end of this
        # function.
        tmp_mem_buff = _setup_lightcone()
    else:
        lib.Co_eval_box_counter = 0

    # Some weird stuff with LOS_direction
    lib.Original_LOS_direction = astro_params.LOS_direction
    lib.Default_LOS_direction = 2
    lib.Stored_LOS_direction_state_1 = astro_params.LOS_direction
    lib.Stored_LOS_direction_state_2 = lib.LOS_direction

    lib.NUM_BINS = int(np.ceil(np.log(lib.HII_DIM) / np.log(
        1.35)))  # This is gotten from the C code -- better than doing the whole init function.
    n_ps = lib.total_num_boxes if flag_options.USE_LIGHTCONE else lib.N_USER_REDSHIFT

    # Have to put these in the DataToBeReturned struct because they are filled there directly.
    PSdata_k = np.zeros(n_ps * lib.NUM_BINS, dtype=np.float32)
    PSdata = np.zeros(n_ps * lib.NUM_BINS, dtype=np.float32)
    lib.DataToBeReturned.PSData_k = ffi.cast("float *", ffi.from_buffer(PSdata_k))
    lib.DataToBeReturned.PSData = ffi.cast("float *", ffi.from_buffer(PSdata))
    lc_box = np.zeros((lib.HII_DIM * lib.HII_DIM * lib.total_slice_ct), dtype=np.float32)
    lib.DataToBeReturned.LCBox = ffi.cast("float*", ffi.from_buffer(lc_box))

    aveNF = _pyalloc("aveNF", lib.N_USER_REDSHIFT)
    aveTb = _pyalloc("aveTb", lib.N_USER_REDSHIFT)

    if flag_options.USE_FCOLL_IONISATION_TABLE and flag_options.USE_LIGHTCONE:
        lib.ReadFcollTable()

    if (lib.R_BUBBLE_MAX > lib.R_MFP_UB) and lib.USE_FCOLL_IONISATION_TABLE:
        raise ValueError(
            """
            The interpolation table for the ionisation box collapse fraction does not have the requisite sampling of the R_MFP.
            (The selected R_MFP exceeds the upper bound on R_MFP set in the interpolation table)
            (Either reduce R_MFP or re-run the creation of the interpolation table with a sufficiently large upper bound)
            """
        )

    if flag_options.USE_TS_FLUCT:
        Ts_z = _pyalloc("Ts_z", lib.HII_TOT_NUM_PIXELS)
        x_e_z = _pyalloc("x_e_z", lib.HII_TOT_NUM_PIXELS)

    # Setup error function arrays, and store them in temporary variable.
    tmp_erf_buff = _setup_erf()

    if flag_options.INHOMO_RECO:
        lib.init_MHR()
        lib.init_inhomo_reco()
    # ==================================================================================================================

    if flag_options.USE_TS_FLUCT:
        print("Computing Ts and Ionization Boxes...")
        lib.ComputeTsBoxes()

    else:
        if flag_options.USE_LIGHTCONE:
            tr = tqdm(range(lib.N_USER_REDSHIFT))
        else:
            tr = tqdm(range(lib.N_USER_REDSHIFT)[::-1])
            delta_T = []

        for i in tr:
            z = lib.redshifts[i]
            tr.set_description("Ionizing z=%.2f" % lib.redshifts[i])
            lib.ComputeIonisationBoxes(i, z, z + 0.2, 0)

            if not flag_options.USE_LIGHTCONE:
                delta_T.append(asarray(ffi, lib.delta_T, lib.HII_TOT_NUM_PIXELS).copy())

    if flag_options.INHOMO_RECO:
        lib.free_MHR()
        lib.destroy_inhomo_reco()

    if flag_options.USE_LIGHTCONE:
        return LightCone(flag_options.redshifts, aveNF, aveTb, PSdata, PSdata_k, lc_box,
                         n_ps, lib.NUM_BINS, lib.HII_DIM, lib.BOX_LEN)
    else:
        return CoEval(flag_options.redshifts, aveNF, aveTb, delta_T, PSdata_k, PSdata, lib.NUM_BINS)


def run_21cmfast(redshift, box_dim={}, flag_options={}, astro_params={}, cosmo_params={},
                 write=True, regenerate=False, run_perturb=True, run_ionize=True, init_boxes=None,
                 free_ps=True):
    """
    A high-level wrapper for 21cmFAST.

    This function runs 21cmFAST for the given redshifts and options. 21cmFAST is broadly broken into three sections:
    init, perturb_field, and ionize. By default all of these are run in sequence. However, the first section(s) can
    be run independently by specifying the run_ parameters.

    By default, all stages of the computation will be written to file to be read in later if required. If these files
    exist already, the default behaviour is to read them in, rather than recompute them. These behaviours can be
    modified with the `write` and `regenerate` parameters. These parameters over-rule the `GenerateNewICs` parameter
    in the `flag_options`, making this parameter redundant (but necessary for now).

    The return value at this point is None unless the ionize part is being run. In that case, it will return either a
    LightCone, or a CoEval box.

    Parameters
    ----------
    redshift : list or float
        The redshift(s) at which to perform 21cmFAST. If float, does a lightcone up to that redshift.

    box_dim : dict
        A dictionary of values to construct the BoxDim object with.

    flag_options : dict
        A dictionary of values to construct the FlagOptions object with.

    astro_params : dict
        A dictionary of values to construct the AstroParams object with.

    cosmo_params : dict
        A dictionary of values to construct the CosmoParams object with.

    write : bool, optional
        Whether to write out the init and perturb boxes.

    regenerate : bool, optional
        Whether to force regeneration of init and perturb boxes (even if they exist). Note that if they do not exist,
        they will be generated regardless of this parameter.

    run_perturb : bool, optional
        Whether to run the perturb_field part of the process.

    run_ionize : bool, optional
        Whether to run the ionization part of the process. Note that if `run_perturb` is False, this is automatically
        set to False.

    init_boxes : boxes, optional
        The boxes returned from :func:`compute_initial_conditions`. If this is passed in, the routine assumes
        that compute_initial_conditions has been called in the particular session, so that various initialising
        routines are *not* re-run. Do *not* just pass in boxes that you read in yourself.

    free_ps : bool, optional
        The initialization routine initialises the power spectrum, but if one wishes to repeatedly call this function
        with the same init_boxes, one should not free the power spectrum variables.

    Returns
    -------
    out : None, LightCone or CoEval
        If Ionization is not run, returns nothing. Otherwise, returns a :class:`LightCone` or :class:`CoEval` object
        depending on input options.
    """
    # ==================================================================================================================
    # Various setting up of parameters
    # ==================================================================================================================
    flag_options['redshifts'] = redshift

    # If we are regenerating initial conditions, and not writing them to file,
    # this indicates that we want to regenerate them *inside* the C code, without writing.
    # This is faster for MCMC where cosmology is being altered.
    # Otherwise, we handle all the generating and writing from Python, but don't pass in perturbed boxes in memory,
    # because that would be a huge amount of memory. Rather read it in within C on every iteration.
    # We *could* read in the files in Python on every iteration and pass by memory, but there's no real benefit to
    # doing this, and it *can't* be done for Ts.
    if regenerate and not write:
        flag_options['GenerateNewICs'] = True
    else:
        flag_options['GenerateNewICs'] = False

    # This *needs* to be done outside the functions, so that the Structure stick around in memory.
    if type(box_dim) != BoxDimStruct:
        box_dim = BoxDimStruct(**box_dim)

    if type(cosmo_params) != CosmoParamStruct:
        cosmo_params = CosmoParamStruct(**cosmo_params)

    if type(astro_params) != AstroParamStruct:
        astro_params = AstroParamStruct(flag_options.get("INHOMO_RECO", False), **astro_params)

    # Set up the various globals that we require.
    if type(flag_options) != FlagOptionStruct:
        flag_options = FlagOptionStruct(astro_params.Z_HEAT_MAX, astro_params.ZPRIME_STEP_FACTOR, **flag_options)

    dir = ffi.string(box_dim.DIREC).decode()
    if not run_perturb:
        run_ionize=False
    # ==================================================================================================================

    # ==================================================================================================================
    # The actual calls to the C code.
    # ==================================================================================================================
    if init_boxes is None:
        init_boxes = compute_initial_conditions(box_dim, cosmo_params, regenerate=regenerate, write=write, dir=dir)

    # We always run perturb if running ionize. If the boxes are already there, we waste minimal time, because all we do
    # is check for their existence. The only time we don't need to run perturb is if we are generating new ICs within
    # the C itself.
    # At this point, this function does *not* ever return the perturbed boxes -- they are only ever a side-effect.
    # Therefore, if it is run, it makes no sense *not* to write them out.
    if run_perturb and not flag_options.GenerateNewICs:
        tr = tqdm(range(flag_options.N_USER_REDSHIFT))
        for i in tr:
            z = flag_options.redshifts[i]
            tr.set_description("Perturbing Field, z=%.2f" % z)
            compute_perturb_field(
                z, init_boxes, velocity_component=astro_params.LOS_direction,
                write=True, regenerate=regenerate, dir=dir, read=False
            )

    if run_ionize:
        print(lib.HII_EFF_FACTOR)
        output = compute_ionisation_boxes(flag_options, astro_params)

        # TODO: This is probably a bad idea. We should probably just try to check if init_ps has been called
        # TODO: before ever calling it. Not sure if this then lends itself to a class setup instead.
        if free_ps:
            lib.free_ps()

        return output, init_boxes

    return init_boxes

def compute_tau(redshifts, xHI, cosmo_params={}):
    if type(cosmo_params) != CosmoParamStruct:
        cosmo_params = CosmoParamStruct(**cosmo_params)

    cosmo_params.set_globals()

    return lib.tau_e(0, redshifts[-1],
                     ffi.cast("float*", ffi.from_buffer(redshifts)),
                     ffi.cast("float*", ffi.from_buffer(xHI)),
                     len(redshifts))
