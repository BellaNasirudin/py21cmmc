import numpy as np
from tqdm import tqdm
import os.path as pth

from . import get_box_parameters, drive_21cmMC, c_generatePS, set_globals, run_init, run_perturb, get_globals


def _read_perturbed(boxdir, redshifts, HII_DIM, BOX_LEN, VelocityComponent):
    # Read in the relevant simulation information from the C code, and then read in all redshift boxes
    # (velocity and density field) generated from 21cmFast.



    DensityBoxes = np.zeros((len(redshifts), HII_DIM * HII_DIM * (HII_DIM // 2 + 1)), np.dtype('complex64'))
    Velocity_Boxes = np.zeros((len(redshifts), HII_DIM * HII_DIM * 2 * (HII_DIM // 2 + 1)), np.dtype('float32'))

    for i, redshift in enumerate(redshifts):

        with open(fn_updated_smoothed(boxdir, redshift, HII_DIM, BOX_LEN), 'rb') as f:
            IndividualBox = np.fromfile(f, dtype=np.dtype('complex64'), count=HII_DIM * HII_DIM * HII_DIM//2)

        for ii in tqdm(range(HII_DIM), desc="Padding Densities "):
            for jj in range(HII_DIM):
                for kk in range((HII_DIM // 2) + 1):
                    if (kk * 2 + 1) < HII_DIM:

                        DensityBoxes[i][kk + ((HII_DIM // 2) + 1) * (jj + HII_DIM * ii)] = IndividualBox[
                                                                                              kk  + (HII_DIM//2) * (
                                                                                                          jj + HII_DIM *
                                                                                                      ii)]
        with open(fn_updated_vel(boxdir, redshift, HII_DIM, BOX_LEN, VelocityComponent), 'rb') as f:
            IndividualBox = np.fromfile(f, dtype=np.dtype('float32'), count=HII_DIM * HII_DIM * HII_DIM)

        for ii in tqdm(range(HII_DIM),desc="Padding Velocities"):
            for jj in range(HII_DIM):
                for kk in range(HII_DIM):
                    Velocity_Boxes[i][kk + 2 * ((HII_DIM // 2) + 1) * (jj + HII_DIM * ii)] = IndividualBox[
                        kk + HII_DIM * (jj + HII_DIM * ii)]

    return DensityBoxes, Velocity_Boxes


def fn_updated_smoothed(boxdir, redshift, HII_DIM, BOX_LEN):
    return pth.expanduser(pth.join(boxdir,'updated_smoothed_deltax_z%06.2f_%i_%.0fMpc'% (redshift, HII_DIM, BOX_LEN)))


def fn_updated_vel(boxdir, redshift, HII_DIM, BOX_LEN, VelocityComponent):
    return pth.expanduser(pth.join(boxdir,'updated_v%s_z%06.2f_%i_%.0fMpc'% ('xyz'[VelocityComponent-1], redshift, HII_DIM, BOX_LEN)))


def fn_deltak(boxdir, DIM, BOX_LEN):
    return pth.expanduser(pth.join(boxdir,"deltak_z0.00_%i_%.0fMpc"%(DIM,BOX_LEN)))


def get_single_box(boxdir, redshifts, zeta, mfp, log10_tvir, init_params=None, generate_ps=False):
    # First, set the INIT global parameters.
    if init_params:
        set_globals(**init_params)

    boxglobs = get_globals()

    box_params = get_box_parameters()

    if not pth.exists(fn_deltak(boxdir, boxglobs.DIM, boxglobs.BOX_LEN)):
        run_init() # No need to pass init params, since we have done that already.

    for i,z in enumerate(redshifts):
        if not pth.exists(fn_updated_smoothed(boxdir, z, box_params.N, box_params.box_len)) or not\
            pth.exists(fn_updated_vel(boxdir, z, box_params.N, box_params.box_len, box_params.vel_comp)):

            run_perturb(z)

    DensityBoxes, Velocity_Boxes = _read_perturbed(boxdir, redshifts, box_params.N, box_params.box_len, box_params.vel_comp)

    # Create the empty delta_T array into which the brightness temperatures will be placed.
    delta_T = np.empty((len(redshifts), box_params.N ** 3), dtype=np.float32)

    output= [0]*len(redshifts)
    PowerSpec = [0]*len(redshifts)
    for i, z in enumerate(redshifts):
        output[i] = drive_21cmMC(DensityBoxes[i], Velocity_Boxes[i],
                              z, zeta, mfp, log10_tvir, 1, delta_T[i])

        if generate_ps:
            PowerSpec[i] = c_generatePS(output[i].ave, output[i].global_xH, delta_T[i])

    if generate_ps:
        return delta_T, output, PowerSpec, box_params
    else:
        return delta_T, output, box_params
