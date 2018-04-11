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
from ._utils import run_single_instance#get_single_box, run_perturb, check_init_existence, check_perturb_existence
from . import plotting as plts
# from ._21cmfast import run_init
import yaml
# import numpy as np
# import importlib
from os.path import join, expanduser
import pickle
from ._21cmfast import AstroParams, CosmoParams, FlagOptions

def _get_config(config=None):
    if config is None:
        config = expanduser(join("~", ".py21cmmc","example_config.yml"))

    with open(config,"r") as f:
        cfg = yaml.load(f)

    return cfg

main = click.Group()

@main.command()
@click.option("--config", type=click.Path(exists=True, dir_okay=False), default=None,
              help="Path to the configuration file (default ~/.py21cmmc/example_config.yml)")
@click.option("--outdir", type=click.Path(exists=True, dir_okay=True), default=None,
              help="directory to write data and plots to -- must exist.")
@click.option("--datafile", type=str, default=None, help="name of outputted datafile (default empty -- no writing)")
@click.option("--plot", multiple=True, help="types of pdf plots to save. Valid values are [global, power, slice]")
def single(config, outdir, datafile, plot):
    "Run a single iteration of 21cmFAST, outputting a series of co-eval boxes, or a lightcone."
    cfg = _get_config(config)

    outputs = run_single_instance(cfg['redshifts'], cfg['flag_options'],
                                  cfg['astro_parameters'], cfg['cosmo_parameters'],
                                  cfg['random_ids'])

    if outdir is None:
        outdir = ""

    if datafile is None:
        datafile = "" #TODO: would be better if we made an automatic filename system.

    if datafile:
        with open(join(outdir, datafile), 'wb') as f:
            pickle.dump(outputs, f)

    if "global" in plot:
        fig, ax = plts.plot_global_data(outputs)
        fig.savefig(join(outdir, datafile+"_global_plot.pdf"))

    if "power" in plot:
        fig, ax = plts.plot_power_spec_1D(outputs)
        fig.savefig(join(outdir, datafile + "_power_plot.pdf"))

    if "slice" in plot:
        fig, ax = plts.plot_lightcone_slice(outputs)
        fig.savefig(join(outdir, datafile + "_slice_plot.pdf"))

@main.command()
def defaults():
    for nm, inst in [("Flag Options", FlagOptions(redshifts=[9.0])),
                   ("Astrophysical Parameters", AstroParams(INHOMO_RECO=False)),
                   ("Cosmological Parameters",CosmoParams())]:
        print(nm+": ")
        print("-"*26)
        for k in inst._fields_:
            print("%-26s: %s"%(k[0], getattr(inst, k[0])))
        print()
#
# @main.command()
# @click.option("--config", type=click.Path(exists=True, dir_okay=False), default=None)
# @click.option("--ps/--no-ps", default=False, help="generate the power spectrum of the box")
# @click.option("--likelihood", default="traditional_c_based", help='the kind of likelihood to use')
# @click.option("--output", type=click.Path(), default=None)
# @click.option("--no-ll/--ll", default=False, help="don't determine the likelihood")
# def single(config, ps, likelihood, output, no_ll):
#     cfg = _get_config(config)
#
#     outputs = get_single_box(cfg['boxdir'], cfg['redshifts'], cfg['zeta_val'],
#                              cfg['MFP_val'], np.log10(cfg['TVir_val']),
#                              init_params=cfg['21cmfast_init'],
#                              generate_ps=ps)
#
#     if ps:
#         delta_T, properties, PS, box_params = outputs
#     else:
#         delta_T, properties, box_params = outputs
#
#     if not no_ll:
#         # Import the correct function
#         module = importlib.import_module("py21cmmc.likelihoods.%s"%likelihood)
#
#         ll = module.get_likelihood(delta_T, box_params, properties, cfg['redshifts'],
#                                    **cfg['likelihood_kwargs']
#                                    )
#
#         print("Log-Likelihood: ", ll)
#
#     if output is not None:
#         with open(output, 'wb') as f:
#             np.save(f, delta_T)

# @main.command()
# @click.option("--config", type=click.Path(exists=True, dir_okay=False), default=None)
# def init(config):
#     cfg = _get_config(config)
#
#     if check_init_existence(cfg['boxdir']):
#         cont = input("A box with these specs already exists. Continue? (y/N)")
#
#     if cont.lower() == 'y':
#         run_init(**cfg['21cmfast_init'])
#     else:
#         print("Exiting without running.")
#
# @main.command()
# @click.option("--config", type=click.Path(exists=True, dir_okay=False), default=None)
# def perturb_field(config):
#     cfg = _get_config(config)
#
#     run_perturb(cfg['boxdir'], cfg['redshifts'], cfg['21cmfast_init'])


