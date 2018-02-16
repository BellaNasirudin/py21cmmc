import numpy as np
from tqdm import tqdm
import time

import os.path as pth

from . import get_box_parameters, drive_21cmMC, c_generatePS, c_likelihood_chisquare


def read_noise(direc, tsc_name, redshifts, zeta, mfp, tvir):
    "Read NoiseData files"
    # Set for the desired telescope ['SKA_halveddipoles_compact', 'HERA331'] corresponding to the file structure in
    # "NoiseData"
    ErrorFileTelescope_Name = '%s' % (tsc_name)

    # Read in the length of the PS files for np arrays. The length of the files should be consistent
    # Note: This should be removed from final version.
    k_values = np.loadtxt(
        pth.expanduser(pth.join(direc, 'MockObs_%s_PS_250Mpc_z%s_%s_%s_%s.txt' % (tsc_name, redshifts[0],
                                                            zeta, mfp, tvir))),
        usecols=(0,)
    )
    multi_z_mockobs_k = np.zeros((len(redshifts), len(k_values)))
    multi_z_mockobs_PS = np.zeros((len(redshifts), len(k_values)))
    multi_z_PS_Error = np.zeros((len(redshifts), len(k_values)))

    # Read in the mock 21cm PS observation. Read in both k and the dimensionless PS
    # Important Note: The mock observations need to be binned "exactly" the same as the resolution of the sampled boxes
    # by 21cmFAST. If not, the chi^2 statistic is still going to be computed, but it will be wrong.
    for i, redshift in enumerate(redshifts):
        multi_z_mockobs_k[i], multi_z_mockobs_PS[i] = np.loadtxt(
            pth.expanduser(pth.join(direc,'MockObs_%s_PS_250Mpc_z%s_%s_%s_%s.txt' % (
            tsc_name, redshift, zeta, mfp, tvir))), usecols=(0, 1), unpack=True)

        # Total noise sensitivity as computed from 21cmSense. This total noise combines both the Poisson component of the
        # observationally measured PS (the high res 21cm FAST simulation) and the instrument thermal noise profile.
        multi_z_PS_Error[i] = np.loadtxt(pth.expanduser(pth.join(direc,'TotalError_%s_PS_250Mpc_z%s_%s_%s_%s.txt' % (
            ErrorFileTelescope_Name, redshift, zeta, mfp, tvir))), usecols=(1,))

    # Bad value flag returns a value far away from the highest likelihood value to ensure the MCMC never enters into
    # this region. This flag will be typically raised for a measured PS from a completely ionised simulation box
    # (neutral fraction zero). Hence, is just the chi^{2} for the PS itself (i.e. PS zero) multiplied by 100 to really
    # ensure this parameter combination can never be chosen
    bad_value_flag = -0.5 * (100. * sum(np.square((multi_z_mockobs_PS[0]) / multi_z_PS_Error[0])))

    return multi_z_mockobs_k, multi_z_mockobs_PS, multi_z_PS_Error


def read_boxes(boxdir, redshifts):
    # Read in the relevant simulation information from the C code, and then read in all redshift boxes
    # (velocity and density field) generated from 21cmFast.

    box_params = get_box_parameters()
    BOX_LEN = box_params.box_len
    HII_DIM = box_params.N
    VelocityComponent = box_params.vel_comp

    DensityBoxes = np.zeros((len(redshifts), HII_DIM * HII_DIM * (HII_DIM // 2 + 1)), np.dtype('complex64'))
    Velocity_Boxes = np.zeros((len(redshifts), HII_DIM * HII_DIM * 2 * (HII_DIM // 2 + 1)), np.dtype('float32'))

    for i, redshift in enumerate(redshifts):

        with open(pth.expanduser(pth.join(boxdir,'updated_smoothed_deltax_z%06.2f_%i_%.0fMpc'% (redshift, HII_DIM, BOX_LEN))), 'rb') as f:
            IndividualBox = np.fromfile(f, dtype=np.dtype('complex64'), count=HII_DIM * HII_DIM * HII_DIM//2)

        for ii in tqdm(range(HII_DIM), desc="Padding Densities "):
            for jj in range(HII_DIM):
                for kk in range((HII_DIM // 2) + 1):
                    if (kk * 2 + 1) < HII_DIM:

                        DensityBoxes[i][kk + ((HII_DIM // 2) + 1) * (jj + HII_DIM * ii)] = IndividualBox[
                                                                                              kk  + (HII_DIM//2) * (
                                                                                                          jj + HII_DIM *
                                                                                                      ii)]
        with open(pth.expanduser(pth.join(boxdir,'updated_v%s_z%06.2f_%i_%.0fMpc' % ('xyz'[VelocityComponent-1], redshift, HII_DIM, BOX_LEN))), 'rb') as f:
            IndividualBox = np.fromfile(f, dtype=np.dtype('float32'), count=HII_DIM * HII_DIM * HII_DIM)

        for ii in tqdm(range(HII_DIM),desc="Padding Velocities"):
            for jj in range(HII_DIM):
                for kk in range(HII_DIM):
                    Velocity_Boxes[i][kk + 2 * ((HII_DIM // 2) + 1) * (jj + HII_DIM * ii)] = IndividualBox[
                        kk + HII_DIM * (jj + HII_DIM * ii)]

    return DensityBoxes, Velocity_Boxes, HII_DIM


def get_inputs(boxdir, redshifts):

    # if len(redshifts) == 1:
    #     multiz_flag = 'singlez'
    # else:
    #     multiz_flag = 'multiz'

    # Set the parameter range for Zeta (ionising efficiency). CosmoHammer takes 4 values, 2nd and 3rd are min and max
    # allowed for the parameter, first and last
    # are not really important. The first is set to the middle of the range and the last is like a one-sigma val. I just
    # set this to be about 1-2 sigma of my parameter range. It is an MCMC, so will sample the entire parameter space.
    # It is most important for the initial locations of the "walkers", but I switched the initial position generator
    # params = np.array(
    #     [
    #         [50., config['zeta_range'][0], config['zeta_range'][1], 20.],
    #         [10., config['mfp_range'][0], config['mfp_range'][1], 3.]
    #         [np.log10(80000.), config['log10tvir_range'][0], config['log10tvir_range'][1], np.log10(50000.)]
    #     ]
    # )

    t0 = time.time()
    DensityBoxes, Velocity_Boxes, HII_DIM = read_boxes(boxdir, redshifts)
    t1 = time.time()

    print("Generating Inputs took: %s sec"%(t1-t0))

    return DensityBoxes, Velocity_Boxes, HII_DIM


def get_single_box(boxdir, redshifts, zeta, mfp, log10_tvir, generate_ps=False):

    DensityBoxes, Velocity_Boxes, HII_DIM = get_inputs(boxdir, redshifts)

    # Create the empty delta_T array into which the brightness temperatures will be placed.
    delta_T = np.empty((len(redshifts), HII_DIM ** 3), dtype=np.float32)

    output= [0]*len(redshifts)
    PowerSpec = [0]*len(redshifts)
    for i, z in enumerate(redshifts):
        output[i] = drive_21cmMC(DensityBoxes[i], Velocity_Boxes[i],
                              z, zeta, mfp, log10_tvir, 1, delta_T[i])

        if generate_ps:
            PowerSpec[i] = c_generatePS(output[i].ave, output[i].global_xH, delta_T[i])

    if generate_ps:
        return delta_T, output, PowerSpec
    else:
        return delta_T, output


def get_likelihood(direc, powerspec, tsc_name, redshifts, zeta, mfp, tvir, **ll_kwargs):
    multi_z_mockobs_k, multi_z_mockobs_PS, multi_z_PS_Error = read_noise(direc, tsc_name, redshifts, zeta, mfp, tvir)

    chi_squared = 0
    for i, (mock, sigma) in enumerate(zip(multi_z_mockobs_PS, multi_z_PS_Error)):
        chi_squared += c_likelihood_chisquare(
            powerspec[i], mock, sigma,
            ll_kwargs['foreground_cut'],
            ll_kwargs['shot_noise_cut'],
            ll_kwargs['ModUncert']
        )

    return chi_squared
