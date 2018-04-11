import numpy as np
from tqdm import tqdm
import os.path as pth
from collections import OrderedDict
import ctypes

from . import drive_21cmMC, AstroParams, CosmoParams, FlagOptions


def _set_flag_options(redshifts, mass_dependent_zeta = False, read_from_file=True, generate_new_ics = False,
                      subcell_rsds=False, ionisation_fcoll_table=False,
                      use_fcoll_table = False, perform_ts_calc=False, inhomo_reco = False, output_global_ave = False,
                      print_data_files = False, print_coeval_boxes=False, print_lightcone_cubes = False):
    """
    Create a dictionary of flag options for 21cmmc, with defaults.

    Parameters
    ----------
    redshifts : list or float
        If a list, the redshifts at which to evaluate co-eval boxes. If a scalar, implies that a light cone is required,
        and its value corresponds to TS_calc_z, i.e. the lowest redshift of the light-cone.

    mass_dependent_zeta : bool, optional
        Should always be zero (have yet to include this option)

    read_from_file : bool, optional
            Whether to read/write Walker file data to file, or pass directly between Python/C.

    generate_new_ics : bool, optional
        Generate on the fly ICs (Only required for 21CMMC). If this flag is not set, it uses the seed in init.c

    subcell_rsds : bool, optional
        Add sub-cell RSDs (currently doesn't work if Ts is not used)

    ionisation_fcoll_table : bool, optional
        An interpolation table for ionisation collapsed fraction (will be removing this at some point)

    use_fcoll_table : bool, optional
        An interpolation table for ionisation collapsed fraction (for Ts.c; will be removing this at some point)

    perform_ts_calc : bool, optional
        Whether to compute the full spin temperature evolution

    inhomo_reco : bool, optional
        Whether to perform inhomogeneous recombinations

    output_global_ave : bool, optional
        Output the global averaged data

    print_data_files : bool, optional
        Output the neutral fraction and 21cm PS text files

    print_coeval_boxes : bool, optional
        Output the co-eval brightness temperature datacubes

    print_lightcone_cubes : bool, optional
        Output the light-cone brightness temperature datacubes

    Returns
    -------
    flag_options : list
        Merely the input, but with defaults set.
    """

    if hasattr(redshifts, "__len__"):
        number_redshifts = len(redshifts)
        lightcone_flag = False
        TS_calc_z = 0 # not required
    else:
        number_redshifts = 0
        lightcone_flag = True
        TS_calc_z = redshifts

    flag_params = OrderedDict(
        number_redshifts = number_redshifts,
        lightcone_flag = int(lightcone_flag),
        mass_dependent_zeta = mass_dependent_zeta,
        TS_calc_z = int(TS_calc_z),
        read_from_file = int(read_from_file),
        generate_new_ics = int(generate_new_ics),
        subcell_rsds = int(subcell_rsds),
        ionisation_fcoll_table = int(ionisation_fcoll_table),
        use_fcoll_table = int(use_fcoll_table),
        perform_ts_calc = int(perform_ts_calc),
        inhomo_reco = int(inhomo_reco),
        output_global_ave = int(output_global_ave),
        print_data_files = int(print_data_files),
        print_coeval_boxes = int(print_coeval_boxes),
        print_lightcone_cubes = int(print_lightcone_cubes)
    )

    return flag_params

def _set_astro_params(INHOMO_RECO,
                      ALPHA = 0.0, ZETA = 30.0, MFP = None, TVIR_MIN = 4.69897, L_X = 40.0, NU_X_THRESH = 500.0,
                      NU_X_BAND_MAX=2000.0, NU_X_MAX = 10000.0, X_RAY_SPEC_INDEX = 1.0, X_RAY_TVIR_MIN = None,
                      X_RAY_TVIR_LB=4.0, X_RAY_TVIR_UB = 6.0, F_STAR = 0.05, t_STAR = 0.5,N_RSD_STEPS = 20,
                      LOS_direction=2, Z_HEAT_MAX=35.0, ZPRIME_STEP_FACTOR = 1.02):
    """

    Parameters
    ----------
    INHOMO_RECO : float, optional
    ALPHA : float, optional
    ZETA : float, optional
    MFP : float, optional
    TVIR_MIN : float, optional
    L_X : float, optional
    NU_X_THRESH : float, optional
    NU_X_BAND_MAX : float, optional
    NU_X_MAX : float, optional
    X_RAY_SPEC_INDEX : float, optional
    X_RAY_TVIR_MIN : float, optional
    X_RAY_TVIR_LB : float, optional
    X_RAY_TVIR_UB : float, optional
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

    Returns
    -------
    astro_params : list
        List which is required to be passed to drive_21cmmc.
    """
    if MFP is None:
        MFP = 15.0
        if INHOMO_RECO == 1:
            MFP = 50.0

    if X_RAY_TVIR_MIN is None:
        X_RAY_TVIR_MIN = TVIR_MIN


    AstrophysicalParameters = OrderedDict(
        ALPHA = ALPHA,
        ZETA = ZETA,
        MFP = MFP,
        TVIR_MIN = TVIR_MIN,
        L_X = L_X,
        NU_X_THRESH = NU_X_THRESH,
        NU_X_BAND_MAX = NU_X_BAND_MAX,
        NU_X_MAX = NU_X_MAX,
        X_RAY_SPEC_INDEX = X_RAY_SPEC_INDEX,
        X_RAY_TVIR_MIN = X_RAY_TVIR_MIN,
        X_RAY_TVIR_LB = X_RAY_TVIR_LB,
        X_RAY_TVIR_UB = X_RAY_TVIR_UB,
        F_STAR = F_STAR,
        t_STAR = t_STAR,
        N_RSD_STEPS = N_RSD_STEPS,
        LOS_direction = LOS_direction,
        Z_HEAT_MAX = Z_HEAT_MAX,
        ZPRIME_STEP_FACTOR = ZPRIME_STEP_FACTOR
    )

    return AstrophysicalParameters


def _set_cosmo_params(RANDOM_SEED=None, SIGMA_8 = 0.820000, littleh = 0.680000,
                      OMEGA_M=0.310000, OMEGA_b = 0.048, NS = 0.97000):
    """

    Parameters
    ----------
    RANDOM_SEED : float, optional
        A seed to set the IC generator. If None, chosen from uniform distribution.

    SIGMA_8 : float, optional
        RMS mass variance (power spectrum normalisation).

    littleh : float, optional
        H_0/100.

    OMEGA_M : float, optional
        Omega matter.

    OMEGA_b : float, optional
        Omega baryon, the baryon component.

    NS : float, optional
        Spectral index of the power spectrum.

    Returns
    -------
    cosmo_params : list
        The cosmological parameters as a list for feeding into drive_21cmmc.
    """

    OMEGA_L = 1.0 - OMEGA_M

    if RANDOM_SEED is None:
        RANDOM_SEED = int(np.random.uniform(1,1e12))

    CosmologicalParameters = OrderedDict(
        RANDOM_SEED = RANDOM_SEED,
        SIGMA_8 = SIGMA_8,
        littleh = littleh,
        OMEGA_M = OMEGA_M,
        OMEGA_L = OMEGA_L,
        OMEGA_b = OMEGA_b,
        NS = NS
    )

    return CosmologicalParameters


def _write_walker_file(redshifts, flag_options, astro_params, random_ids= ('3.000000', '3.000000')):
    separator = " "
    separator_other = "_"
    seq = []

    #TODO: this is probably better as a hash or something. Add the random thread ID
    seq.append("%s" % (random_ids[0]))
    # Add the second ID
    seq.append("%s" % (random_ids[1]))


    StringArgument_other = separator_other.join(seq)

    # Add number of redshifts
    # If using the light-cone version of the code, don't need to set a redshift
    seq.append("%d" % (flag_options.N_USER_REDSHIFTS))
    # Add light cone flag
    seq.append("%d" % (flag_options.USE_LIGHTCONE))
    # If power-law dependence on ionising efficiency is allowed. Add the flag here (support not yet included)
    seq.append("0")
    # Add redshift for Ts.c calculation
    seq.append("%g" % (flag_options.REDSHIFT))

    StringArgument = separator.join(seq)

    with open("Walker_%s.txt" % (StringArgument_other), "w") as f:
        f.write("FLAGS")
        for key in flag_options._fields_:
            if key[0] not in ['N_USER_REDSHIFTS', 'USE_LIGHTCONE', 'REDSHIFT']:
                f.write("    %s"%(getattr(flag_options, key[0])))
        f.write("\n")

        for key in astro_params._fields_:
            f.write("%s    %f\n"%(key[0],getattr(astro_params, key[0])))

        if hasattr(redshifts, "__len__"):
            for z in redshifts:
                f.write("CO-EVAL-Z    %f\n" % (float(z)))


def run_single_instance(redshifts, flag_options=None, astro_parameters = None, cosmo_parameters = None,
                        random_ids= ('3.000000', '3.000000') ):

    # We don't want to use empty dict as a default, as it is mutable. So set it here.
    flag_options = flag_options or {}
    astro_parameters = astro_parameters or {}
    cosmo_parameters = cosmo_parameters or {}

    # Turn each of the passed dictionaries into Ctypes structures.
    flag_options = FlagOptions(redshifts, **flag_options)
    astro_parameters = AstroParams(flag_options.INHOMO_RECO, **astro_parameters)
    cosmo_parameters = CosmoParams(**cosmo_parameters)

    if flag_options.READ_FROM_FILE:
        _write_walker_file(redshifts, flag_options, astro_parameters, random_ids)

    # flag_options_c = (ctypes.c_float * len(flag_options))(*flag_options.values())
    # astro_params_c = (ctypes.c_float * len(astro_parameters))(*astro_parameters.values())
    # cosmo_params_c = (ctypes.c_float * len(cosmo_parameters))(*cosmo_parameters.values())

    #TODO: should add the options to do the init or create_dens_veloc like in the script (actually better to do in
    #      another function)

    output = drive_21cmMC(
        random_ids,
        flag_options,
        astro_parameters,
        cosmo_parameters
    )

    return output


# def _read_perturbed(boxdir, redshifts, HII_DIM, BOX_LEN, VelocityComponent):
#     # Read in the relevant simulation information from the C code, and then read in all redshift boxes
#     # (velocity and density field) generated from 21cmFast.
#
#
#
#     DensityBoxes = np.zeros((len(redshifts), HII_DIM * HII_DIM * (HII_DIM // 2 + 1)), np.dtype('complex64'))
#     Velocity_Boxes = np.zeros((len(redshifts), HII_DIM * HII_DIM * 2 * (HII_DIM // 2 + 1)), np.dtype('float32'))
#
#     for i, redshift in enumerate(redshifts):
#
#         with open(fn_updated_smoothed(boxdir, redshift, HII_DIM, BOX_LEN), 'rb') as f:
#             IndividualBox = np.fromfile(f, dtype=np.dtype('complex64'), count=HII_DIM * HII_DIM * HII_DIM//2)
#
#         for ii in tqdm(range(HII_DIM), desc="Padding Densities "):
#             for jj in range(HII_DIM):
#                 for kk in range((HII_DIM // 2) + 1):
#                     if (kk * 2 + 1) < HII_DIM:
#
#                         DensityBoxes[i][kk + ((HII_DIM // 2) + 1) * (jj + HII_DIM * ii)] = IndividualBox[
#                                                                                               kk  + (HII_DIM//2) * (
#                                                                                                           jj + HII_DIM *
#                                                                                                       ii)]
#         with open(fn_updated_vel(boxdir, redshift, HII_DIM, BOX_LEN, VelocityComponent), 'rb') as f:
#             IndividualBox = np.fromfile(f, dtype=np.dtype('float32'), count=HII_DIM * HII_DIM * HII_DIM)
#
#         for ii in tqdm(range(HII_DIM),desc="Padding Velocities"):
#             for jj in range(HII_DIM):
#                 for kk in range(HII_DIM):
#                     Velocity_Boxes[i][kk + 2 * ((HII_DIM // 2) + 1) * (jj + HII_DIM * ii)] = IndividualBox[
#                         kk + HII_DIM * (jj + HII_DIM * ii)]
#
#     return DensityBoxes, Velocity_Boxes
#
#
# def fn_updated_smoothed(boxdir, redshift, HII_DIM, BOX_LEN):
#     return pth.expanduser(pth.join(boxdir,'updated_smoothed_deltax_z%06.2f_%i_%.0fMpc'% (redshift, HII_DIM, BOX_LEN)))
#
#
# def fn_updated_vel(boxdir, redshift, HII_DIM, BOX_LEN, VelocityComponent):
#     return pth.expanduser(pth.join(boxdir,'updated_v%s_z%06.2f_%i_%.0fMpc'% ('xyz'[VelocityComponent-1], redshift, HII_DIM, BOX_LEN)))
#
#
# def fn_deltak(boxdir, DIM, BOX_LEN):
#     return pth.expanduser(pth.join(boxdir,"deltak_z0.00_%i_%.0fMpc"%(DIM,BOX_LEN)))
#
# def check_init_existence(boxdir):
#     boxglobs = get_globals()
#     return pth.exists(fn_deltak(boxdir, boxglobs.DIM, boxglobs.BOX_LEN))
#
# def check_perturb_existence(boxdir, z):
#     box_params = get_box_parameters()
#     return (pth.exists(fn_updated_smoothed(boxdir, z, box_params.N, box_params.box_len)) and
#             pth.exists(fn_updated_vel(boxdir, z, box_params.N, box_params.box_len, box_params.vel_comp)))
#
# def run_perturb(boxdir, redshifts,init_params):
#     if init_params:
#         set_globals(**init_params)
#
#     boxglobs = get_globals()
#
#     box_params = get_box_parameters()
#
#     if not check_init_existence(boxdir):
#         run_init()  # No need to pass init params, since we have done that already.
#
#     for i, z in enumerate(redshifts):
#         if not check_perturb_existence(boxdir, z):
#             run_perturb(z)
#
#
#
# def get_single_box(boxdir, redshifts, zeta, mfp, log10_tvir, init_params=None, generate_ps=False):
#     # First, set the INIT global parameters.
#     if init_params:
#         set_globals(**init_params)
#
#     boxglobs = get_globals()
#
#     box_params = get_box_parameters()
#
#     if not pth.exists(fn_deltak(boxdir, boxglobs.DIM, boxglobs.BOX_LEN)):
#         run_init() # No need to pass init params, since we have done that already.
#
#     for i,z in enumerate(redshifts):
#         if not pth.exists(fn_updated_smoothed(boxdir, z, box_params.N, box_params.box_len)) or not\
#             pth.exists(fn_updated_vel(boxdir, z, box_params.N, box_params.box_len, box_params.vel_comp)):
#
#             run_perturb(z)
#
#     DensityBoxes, Velocity_Boxes = _read_perturbed(boxdir, redshifts, box_params.N, box_params.box_len, box_params.vel_comp)
#
#     # Create the empty delta_T array into which the brightness temperatures will be placed.
#     delta_T = np.empty((len(redshifts), box_params.N ** 3), dtype=np.float32)
#
#     output= [0]*len(redshifts)
#     PowerSpec = [0]*len(redshifts)
#     for i, z in enumerate(redshifts):
#         output[i] = drive_21cmMC(DensityBoxes[i], Velocity_Boxes[i],
#                               z, zeta, mfp, log10_tvir, 1, delta_T[i])
#
#         if generate_ps:
#             PowerSpec[i] = c_generatePS(output[i].ave, output[i].global_xH, delta_T[i])
#
#     if generate_ps:
#         return delta_T, output, PowerSpec, box_params
#     else:
#         return delta_T, output, box_params
