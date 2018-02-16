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
from ._utils import get_single_box, get_likelihood
import os.path as pth
import yaml
import numpy as np

main = click.Group()

@main.command()
@click.argument("config", type=click.Path(exists=True, dir_okay=False))
@click.option("--ps/--no-ps", default=False, help="generate the power spectrum of the box")
@click.option("--ll/--no-ll", default=False, help='determine the likelihood (eg. chi-square) against data')
@click.option("--output", type=click.Path(), default=None)
def single(config, ps, ll, output):
    with open(config,"r") as f:
        cfg = yaml.load(f)

    outputs = get_single_box(cfg['boxdir'], cfg['redshifts'], cfg['zeta_val'],
                             cfg['MFP_val'], np.log10(cfg['TVir_val']),
                             generate_ps=ps)

    if ps:
        delta_T, properties, PS = outputs
    else:
        delta_T, properties = outputs

    if ll:
        likelihood = get_likelihood(cfg['noisedir'], PS, cfg['telescope_name'], cfg['redshifts'],
                                    cfg['zeta_val'], cfg['MFP_val'], cfg['TVir_val'],
                                    foreground_cut = cfg['foreground_cut'],
                                    shot_noise_cut= cfg['shot_noise_cut'],
                                    ModUncert = cfg['ModUncert'])

        print("Log-Likelihood: ", likelihood)

    if output is not None:
        with open(output, 'wb') as f:
            np.save(f, delta_T)
