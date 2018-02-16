from py21cmmc import c_likelihood_chisquare, c_generatePS
import numpy as np
import os.path as pth


def get_likelihood(delta_T, box_properties, physical_properties, redshifts, **ll_kwargs):
    """
    Log-Likelihood generator using a simple power-spectrum estimate, based on the original C code in 21CMMC.

    Note, each pluggable likelihood function must have the same signature as this function.

    Parameters
    ----------
    delta_T : nd-array
        The differential brightness temperature box. It is a single-dimensional array of shape N**3, and N can
        be found in `box_properties`.
    box_properties : dict
        A dictionary of properties of the output box.
    physical_properties : list of dict-like
        A list of dictionary-like objects of physical properties of the box(one for each redshift). Each entry of the
        list has the following entries:
            * ave : average temp
            * global_xH : global ionized fraction.
    redshifts : list
        A list of redshifts that have been determined.
    ll_kwargs :
        These are any other options to the likelihood function.

    Returns
    -------
    ll : float
        The log-likelihood.
    """
    # Unpack the kwargs for this likelihood
    direc = ll_kwargs['noisedir']
    tsc_name = ll_kwargs['telescope_name']
    zeta = ll_kwargs['zeta']
    mfp = ll_kwargs['mfp']
    tvir = ll_kwargs['tvir']

    multi_z_mockobs_k, multi_z_mockobs_PS, multi_z_PS_Error = read_noise(direc, tsc_name, redshifts, zeta, mfp, tvir)

    # Generate power spectra in C
    powerspec = [0] * len(redshifts)
    for i in range(len(redshifts)):
        powerspec[i] = c_generatePS(physical_properties[i].ave, physical_properties[i].global_xH, delta_T[i])

    chi_squared = 0
    for i, (ps, mock, sigma) in enumerate(zip(powerspec, multi_z_mockobs_PS, multi_z_PS_Error)):
        chi_squared += c_likelihood_chisquare(
            ps, mock, sigma,
            ll_kwargs['foreground_cut'],
            ll_kwargs['shot_noise_cut'],
            ll_kwargs['ModUncert']
        )

    return -chi_squared


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
    # bad_value_flag = -0.5 * (100. * sum(np.square((multi_z_mockobs_PS[0]) / multi_z_PS_Error[0])))

    return multi_z_mockobs_k, multi_z_mockobs_PS, multi_z_PS_Error



