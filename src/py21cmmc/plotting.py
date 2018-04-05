import pylab as P
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

EoR_color = matplotlib.colors.LinearSegmentedColormap.from_list('mycmap', [(0, 'white'), (0.21, 'yellow'),
                                                                                (0.42, 'orange'), (0.63, 'red'),
                                                                                (0.86, 'black'), (0.9, 'blue'),
                                                                                (1, 'cyan')])
plt.register_cmap(cmap=EoR_color)


def plot_global_data(lightcone, cmap=None, fig=None, ax = None,
                	 min_val = -250, max_val = 50.):
    if cmap is None:
        cmap = EoR_color

    z = lightcone.redshifts

    if fig is None:
        fig, ax = plt.subplots(1,2, figsize=(12, 1.9))

    ax[0].axis([z[0], z[-1], 0, 1])

    if len(z) < 10:
        ax[0].plot(z, lightcone.average_nf, 'ko')
    else:
        ax[0].plot(z, lightcone.average_nf)

    ax[1].axis([z[0], z[-1], min_val, max_val])

    if len(z) < 10:
        ax.plot(z, lightcone.average_Tb, 'ko')
    else:
        ax.plot(z, lightcone.average_Tb)

    return fig, ax


def plot_power_spec_1D(lightcone,fig=None, ax = None, PSmin = 0.01,PSmax = 100):
    if fig is None:
        fig, ax = plt.subplots(np.ceil(np.sqrt(lightcone.n_ps)),np.ceil(np.sqrt(lightcone.n_ps)), figsize=(12, 8.9),
                               subplot_kw={"xscale":'log', "yscale":'log'})

    for i,axi in enumerate(ax.flatten()):
        axi.axis([lightcone.k[0], lightcone.k[-1], PSmin, PSmax])
        axi.plot(lightcone.k, lightcone.power_spectrum[i])

    return fig, ax


def plot_lightcone_slice(lightcone, slice_index=None, fig=None, ax=None, min_val=-250., max_val=50.0, cmap=None):
    if not lightcone.lightcone:
        raise ValueError("This 21cmFAST output is not a lightcone. Aborting.")

    if cmap is None:
        cmap = EoR_color

    # The random voxel used for creating the light-cone slice
    RAND_DIR = lightcone.HII_DIM//2

    # Array to hold the light-cone slice
    SliceData = lightcone.lightcone_box[RAND_DIR]

    X_Vals = np.zeros(lightcone.HII_DIM)
    LC_Vals = np.zeros(lightcone.ncells_los)

    #TODO: should be able to use linspace
    for ii in range(lightcone.HII_DIM):
        X_Vals[ii] = (0.0 + lightcone.box_length * (float(ii) + 0.5) / ((lightcone.HII_DIM - 1)))


    # Indexing for the voxels along the line-of-sight (should be redshift, could return that data too...)
    for ii in range(lightcone.ncells_los):
        LC_Vals[ii] = ii

    if fig is None:
        fig, ax = plt.subplots(1,1, figsize=(12, 1.9))


    ax.axis([LC_Vals[0], LC_Vals[-1], X_Vals[0], X_Vals[-1]])

    CS = ax.pcolormesh(LC_Vals, X_Vals, SliceData, vmin=min_val, vmax=max_val, cmap=cmap, shading='gouraud',
                       zorder=-2)


    #TODO: really unsure why we're doing this, not just calling colorbar()
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size="1%", pad=0.05, pack_start=False)
    fig.add_axes(cax)

    m = plt.cm.ScalarMappable(cmap=cmap)
    m.set_array(SliceData)
    m.set_clim(vmin=min_val, vmax=max_val)

    cbar = plt.colorbar(m, cax=cax)

    cbar.set_label(r"$\boldsymbol{\delta T_{\rm b}\,[{\rm mK}]}$", rotation=270, labelpad=2, fontsize='large')

    cbar.set_ticks(np.arange(-250., 51., 50.))
    cbar.set_ticklabels([r"$\,$", r"$\,$", r"$-150$", r"$\,$", r"$-50$", r"$\,$", r"$50$"])

    cbar.ax.tick_params(axis='y', direction='out', pad=2)

    return fig, ax
