"""
The idea here is to create some data with py21cmmc, and then fit it with the same. We use no foregrounds or
anything like that here.
"""
import os
from py21cmmc.likelihood import Likelihood1DPowerLightconeNoErrors as Likelihood
from py21cmmc.mcmc import run_mcmc

# ====== Manually set parameters for the run =================================
force_redraw = True

parameters = {
    "HII_EFF_FACTOR": ['alpha', 30.0, 10.0, 50.0, 3.0],
    "R_BUBBLE_MAX": ['mfp', 50.0, 30, 70.0, 3.0],
    "ION_Tvir_MIN": ['tvir', 4.69897, 4., 6.0, 0.1]
}

storage_options = {
    "DATADIR": os.path.expanduser("~/Documents/MCMCDataSimplestTest"),
    "KEEP_ALL_DATA": False,
    "KEEP_GLOBAL_DATA": False,
}

box_dim = {
    "HII_DIM": 75,
    "BOX_LEN": 150.0
}

flag_options = {
    'redshifts': 7.0
}

astro_params = {
    "HII_EFF_FACTOR": 30.0,
    "R_BUBBLE_MAX": 50.0,
    "ION_Tvir_MIN": 4.69897
}

cosmo_params = {}
# ============================================================================


filename_of_data = os.path.join(storage_options['DATADIR'], "data.txt")
box_dim['DIREC'] = storage_options['DATADIR']

try:
    os.mkdir(box_dim['DIREC'])
except:
    pass


lk = Likelihood(filename_of_data, box_dim=box_dim, flag_options=flag_options, astro_params=astro_params,
                cosmo_params=cosmo_params, n_psbins=20)

if not os.path.exists(filename_of_data) or force_redraw:
    lk.simulate_data(save_lightcone=True)

run_mcmc(
    redshift=flag_options['redshifts'],
    parameters=parameters,
    storage_options = storage_options,
    box_dim=box_dim,
    flag_options=flag_options,
    astro_params=astro_params,
    cosmo_params=cosmo_params,
    likelihood_modules=[lk],
    walkersRatio=10,
    burninIterations=0,
    sampleIterations=100,
    threadCount=6
)
