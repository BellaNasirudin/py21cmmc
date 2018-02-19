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
from ._utils import get_single_box
from ._21cmfast import set_globals
import yaml
import numpy as np
import importlib

main = click.Group()

@main.command()
@click.argument("config", type=click.Path(exists=True, dir_okay=False))
@click.option("--ps/--no-ps", default=False, help="generate the power spectrum of the box")
@click.option("--likelihood", default="traditional_c_based", help='the kind of likelihood to use')
@click.option("--output", type=click.Path(), default=None)
@click.option("--no-ll/--ll", default=False, help="don't determine the likelihood")
def single(config, ps, likelihood, output, no_ll):
    with open(config,"r") as f:
        cfg = yaml.load(f)

    outputs = get_single_box(cfg['boxdir'], cfg['redshifts'], cfg['zeta_val'],
                             cfg['MFP_val'], np.log10(cfg['TVir_val']),
                             init_params=cfg['21cmfast_init'],
                             generate_ps=ps)

    if ps:
        delta_T, properties, PS, box_params = outputs
    else:
        delta_T, properties, box_params = outputs

    if not no_ll:
        # Import the correct function
        module = importlib.import_module("py21cmmc.likelihoods.%s"%likelihood)

        ll = module.get_likelihood(delta_T, box_params, properties, cfg['redshifts'],
                                   **cfg['likelihood_kwargs']
                                   )

        print("Log-Likelihood: ", ll)

    if output is not None:
        with open(output, 'wb') as f:
            np.save(f, delta_T)
