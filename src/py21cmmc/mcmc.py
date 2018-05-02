from .likelihood import Core21cmFastModule
import sys
from os import path, mkdir
from cosmoHammer import LikelihoodComputationChain, CosmoHammerSampler
from cosmoHammer.util import Params


def run_mcmc(redshift, parameters, storage_options, box_dim = {}, flag_options={}, astro_params={}, cosmo_params={},
             extra_core_modules=[], likelihood_modules = [], **mcmc_options):

    flag_options['redshifts'] = redshift

    try:
        mkdir(storage_options['DATADIR'])
    except:
        pass

    # Setup parameters.
    params = Params(*[(k, v[1:]) for k, v in parameters.items()])

    # Setup the Core Module
    core_module = Core21cmFastModule(
        parameters, storage_options=storage_options, box_dim = box_dim, flag_options=flag_options,
        astro_params=astro_params, cosmo_params=cosmo_params
    )

    # Get all the likelihood modules.
    likelihoods = []
    for lk in likelihood_modules:
        if isinstance(lk, tuple):
            if isinstance(lk[0], str):
                likelihoods += [getattr(sys.modules['py21cmmc.likelihood'],
                                        "Likelihood%s" % lk[0])(box_dim = box_dim,
                                                                astro_params=astro_params,
                                                                cosmo_params=cosmo_params,
                                                                flag_options=flag_options,
                                                                **lk[1])]
            else:
                likelihoods += [lk[0](**lk[1])]

        else:
            likelihoods += [lk]

    chain = LikelihoodComputationChain(min=params[:, 1], max=params[:, 2])

    # Add default plus extra core modules.
    chain.addCoreModule(core_module)
    for cm in extra_core_modules:
        chain.addCoreModule(cm)

    for lk in likelihoods:
        chain.addLikelihoodModule(lk)
    chain.setup()

    sampler = CosmoHammerSampler(
        params=params,
        likelihoodComputationChain=chain,
        filePrefix=path.join(storage_options['DATADIR'], "mcmc_"),
        **mcmc_options
    )

    # The sampler writes to file, so no need to save anything ourselves.
    sampler.startSampling()
