"""
Module that contains the command line app.

Why does this file exist, and why not put this in __main__?

  You might be tempted to import things from __main__ later, but that will cause
  problems: the code will get executed twice:

  - When you run `python -mpy21cmmc` python will execute
    ``__main__.py`` as a script. That means there won't be any
    ``py21cmmc.__main__`` in ``sys.modules``.
  - When you import __main__ it will get executed again (as a module) because
    there's no ``py21cmmc.__main__`` in ``sys.modules``.

  Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""
import click
#from ._utils import run_single_instance#get_single_box, run_perturb, check_init_existence, check_perturb_existence
from . import plotting as plts
import yaml
import pickle
import py21cmmc as p21c

from .likelihood import Core21cmFastModule
import sys
from os import path
from cosmoHammer import LikelihoodComputationChain, CosmoHammerSampler
from cosmoHammer.util import Params
import os

def _get_config(config=None):
    if config is None:
        config = path.expanduser(path.join("~", ".py21cmmc","example_config.yml"))

    with open(config,"r") as f:
        cfg = yaml.load(f)

    return cfg


main = click.Group()


def _single(config=None, write=True, regen=False, outdir="", datafile=None, plot=[],
            run_perturb=True, run_ionize=True):
    cfg = _get_config(config)

    if outdir is None: outdir = "."
    cfg['box_dim']['DIREC'] = outdir

    outputs = p21c.run_21cmfast(
        cfg['redshifts'], cfg['box_dim'], cfg['flag_options'], cfg['astro_parameters'],
        cfg['cosmo_parameters'], write=write, regenerate=regen, run_perturb=run_perturb, run_ionize=run_ionize
    )

    if outputs is not None: #outputs is None if ionization is not run.
        if datafile is None:
            datafile = "" #TODO: would be better if we made an automatic filename system.

        if datafile:
            with open(path.join(outdir, datafile), 'wb') as f:
                pickle.dump(outputs, f)

        if "global" in plot:
            fig, ax = plts.plot_global_data(outputs)
            fig.savefig(path.join(outdir, datafile+"_global_plot.pdf"))

        if "power" in plot:
            fig, ax = plts.plot_power_spec_1D(outputs)
            fig.savefig(path.join(outdir, datafile + "_power_plot.pdf"))

        if "slice" in plot:
            fig, ax = plts.plot_lightcone_slice(outputs)
            fig.savefig(path.join(outdir, datafile + "_slice_plot.pdf"))


@main.command()
@click.option("--config", type=click.Path(exists=True, dir_okay=False), default=None,
              help="Path to the configuration file (default ~/.py21cmmc/example_config.yml)")
@click.option('--write/--no-write', default=True,
              help="Whether to write out intermediate files (from init and perturb_field) for later use")
@click.option("--regen/--no-regen", default=False,
              help="Whether to force regeneration of init/perturb files if they already exist.")
@click.option("--outdir", type=click.Path(exists=True, dir_okay=True), default=None,
              help="directory to write data and plots to -- must exist.")
@click.option("--datafile", type=str, default=None, help="name of outputted datafile (default empty -- no writing)")
@click.option("--plot", multiple=True, help="types of pdf plots to save. Valid values are [global, power, slice]")
@click.option("--perturb/--no-perturb", default = True,
              help="Whether to run the perturbed field calculation")
@click.option("--ionize/--no-ionize", default = True,
              help="Whether to run the ionization calculation")
def single(config, write, regen, outdir, datafile, plot, perturb, ionize):
    return _single(config, write, regen, outdir, datafile, plot, perturb, ionize)


@main.command()
@click.option("--config", type=click.Path(exists=True, dir_okay=False), default=None,
              help="Path to the configuration file (default ~/.py21cmmc/example_config.yml)")
@click.option("--regen/--no-regen", default=False,
              help="Whether to force regeneration of init/perturb files if they already exist.")
@click.option("--outdir", type=click.Path(exists=True, dir_okay=True), default=None,
              help="directory to write data and plots to -- must exist.")
def init(config, regen, outdir):
    """
    Run a single iteration of 21cmFAST init, saving results to file.
    The same operation can be done with ``py21cmmc single --no-perturb``.
    """
    _single(config, write=True, regen=regen, outdir=outdir, datafile=None, plot=[], run_perturb=False, run_ionize=False)


@main.command()
@click.option("--config", type=click.Path(exists=True, dir_okay=False), default=None,
              help="Path to the configuration file (default ~/.py21cmmc/example_config.yml)")
@click.option("--regen/--no-regen", default=False,
              help="Whether to force regeneration of init/perturb files if they already exist.")
@click.option("--outdir", type=click.Path(exists=True, dir_okay=True), default=None,
              help="directory to write data to.")
def perturb_field(config, regen, outdir):
    "Run a single iteration of 21cmFAST init, saving results to file."
    _single(config, write=True, regen=regen, outdir=outdir, datafile=None, plot=[], run_perturb=True, run_ionize=False)



@main.command()
@click.option("--config", type=click.Path(exists=True, dir_okay=False), default=None,
              help="Path to the configuration file (default ~/.py21cmmc/example_config_mcmc.yml)")
def mcmc(config):
    cfg = _get_config(config or path.expanduser(path.join("~", ".py21cmmc","example_config_mcmc.yml")))
    cfg['flag_options'].update(redshifts=cfg['redshifts'])

    try:
        os.mkdir(cfg['storage_options']['DATADIR'])
    except:
        pass

    # Setup parameters.
    params = Params(*[(k, v[1:]) for k, v in cfg['parameters'].items()])

    # Setup the Core Module
    core_module = Core21cmFastModule(
        cfg['parameters'], cfg['storage_options'],cfg['box_dim'], cfg['flag_options'],
        cfg['astro_parameters'], cfg['cosmo_parameters']
    )

    # Get all the likelihood modules.
    likelihoods = [getattr(sys.modules['py21cmmc.likelihood'], "Likelihood%s"%k)(**v) for k,v in cfg['likelihoods'].items()]

    chain = LikelihoodComputationChain(min=params[:,1], max = params[:,2])
    chain.addCoreModule(core_module)
    for lk in likelihoods:
        chain.addLikelihoodModule(lk)
    chain.setup()

    sampler = CosmoHammerSampler(
        params= params,
        likelihoodComputationChain=chain,
        filePrefix=path.join(cfg['storage_options']['DATADIR'], "ReionModel_21cmFAST"),
        walkersRatio = cfg['walkersRatio'],
        burninIterations=cfg['burninIterations'],
        sampleIterations=cfg['sampleIterations'],
        threadCount = cfg['threadCount']
    )

    # The sampler writes to file, so no need to save anything ourselves.
    sampler.startSampling()


@main.command()
def defaults():
    for nm, inst in [("Flag Options", p21c.FlagOptionStruct(redshifts=[9.0])),
                     ("Astrophysical Parameters", p21c.AstroParamStruct(INHOMO_RECO=False)),
                     ("Cosmological Parameters",p21c.CosmoParamStruct()),
                     ("Box Dimensions", p21c.BoxDimStruct())]:
        print(nm+": ")
        print("-"*26)
        for k in p21c.ffi.typeof(inst.cstruct[0]).fields:
            print("%-26s: %s"%(k[0], getattr(inst, k[0])))
        print()
