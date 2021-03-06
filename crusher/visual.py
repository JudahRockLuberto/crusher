#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Visualization of the results.""" 

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Ellipse
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import NullFormatter, FormatStrFormatter
from matplotlib import rcParams
from matplotlib.colors import Normalize
import matplotlib.patches as mpatches

from kungpao.display import display_single

from crusher.data import Z_SUN

plt.rc('text', usetex=True)
rcParams.update({'axes.linewidth': 1.5})
rcParams.update({'xtick.direction': 'in'})
rcParams.update({'ytick.direction': 'in'})
rcParams.update({'xtick.minor.visible': 'True'})
rcParams.update({'ytick.minor.visible': 'True'})
rcParams.update({'xtick.major.pad': '7.0'})
rcParams.update({'xtick.major.size': '8.0'})
rcParams.update({'xtick.major.width': '1.5'})
rcParams.update({'xtick.minor.pad': '7.0'})
rcParams.update({'xtick.minor.size': '4.0'})
rcParams.update({'xtick.minor.width': '1.5'})
rcParams.update({'ytick.major.pad': '7.0'})
rcParams.update({'ytick.major.size': '8.0'})
rcParams.update({'ytick.major.width': '1.5'})
rcParams.update({'ytick.minor.pad': '7.0'})
rcParams.update({'ytick.minor.size': '4.0'})
rcParams.update({'ytick.minor.width': '1.5'})
rcParams.update({'axes.titlepad': '10.0'})
rcParams.update({'font.size': 25})

__all__ = [
    'show_maps',
    'show_aper',
    'prepare_show_ellipse',
    'overplot_ellipse',
    'plot_ell_prof',
    'plot_ell_fourier',
]

# Color maps
IMG_CMAP = plt.get_cmap('Greys')
IMG_CMAP.set_bad(color='w')


def show_maps(maps, aper, age_info=True, met_info=True, cid=None, logms=None, figsize=(15, 15)):
    """Visualize the stellar mass, age, metallicity, and velocity dispersion maps.

    Parameters
    ----------
    maps : dict
        Dictionary that contains all stellar mass, age, and metallicity maps.
    aper : dict
        Dictionary that contains basic shape information of the galaxy.
    cid : int, optional
        `catsh_id`, sub-halo ID in the simulation. Used to identify galaxy.
        Default: None
    logms : float, optional
        Stellar mass in log10 unit. Default: None.
    figsize : tuple, optional
        Size of the 3x3 figure. Default: (15, 15)

    """
    # Setup the figure and grid of axes
    fig_sum = plt.figure(figsize=figsize, constrained_layout=False)
    # add rows based on if age or met is True (hence 1+age+met)
    grid_sum = fig_sum.add_gridspec(1+age_info+met_info, 3, wspace=0.0, hspace=0.0)
    fig_sum.subplots_adjust(
        left=0.005, right=0.995, bottom=0.005, top=0.995,
        wspace=0.00, hspace=0.00)

    # List of the maps need to be plot
    list_maps = ['mass_gal', 'mass_ins', 'mass_exs']
    
    if age_info == True:
        list_maps += ['age_gal', 'age_ins', 'age_exs']
        
    if met_info == True:
        list_maps += ['met_gal', 'met_ins', 'met_exs']

    for ii, name in enumerate(list_maps):
        ax = fig_sum.add_subplot(grid_sum[ii])
        if 'mass' in name:
            if ii % 3 == 0:
                _ = display_single(
                    maps[name], ax=ax, stretch='log10', zmin=6.0, zmax=10.5,
                    color_bar=True, scale_bar=False, no_negative=True,
                    color_bar_height='5%', color_bar_width='85%', color_bar_fontsize=20,
                    cmap=IMG_CMAP, color_bar_color='k')
                _ = ax.text(0.05, 0.06, r'$\log[M_{\star}/M_{\odot}]$', fontsize=25,
                            transform=ax.transAxes,
                            bbox=dict(facecolor='w', edgecolor='none', alpha=0.8))
            else:
                _ = display_single(
                    maps[name], ax=ax, stretch='log10', zmin=6.0, zmax=10.5,
                    color_bar=False, scale_bar=False, no_negative=True,
                    cmap=IMG_CMAP, color_bar_color='k')
            if ii == 0:
                # Label the center of the galaxy
                ax.scatter(aper['x'], aper['y'], marker='+', s=200,
                           c='orangered', linewidth=2.0, alpha=0.6)
                # Show the isophote shape
                e = Ellipse(xy=(aper['x'], aper['y']),
                            height=80.0 * aper['ba'], width=80.0, angle=aper['pa'])
                e.set_facecolor('none')
                e.set_edgecolor('orangered')
                e.set_alpha(0.5)
                e.set_linewidth(2.0)
                ax.add_artist(e)
                # Central
                ax.text(0.75, 0.06, r'$\rm Total$', fontsize=25,
                        transform=ax.transAxes)
            # Put the ID
            if ii == 1 and cid is not None and logms is not None:
                ax.text(
                    0.5, 0.88, r'$\mathrm{ID}: %d\ \ \log M_{\star}: %5.2f$' % (cid, logms),
                    fontsize=25, transform=ax.transAxes, horizontalalignment='center',
                    bbox=dict(facecolor='w', edgecolor='none', alpha=0.7))
                ax.text(0.75, 0.06, r'$\rm In\ situ$', fontsize=25,
                        transform=ax.transAxes)
            if ii == 2:
                ax.text(0.75, 0.06, r'$\rm Ex\ situ$', fontsize=25,
                        transform=ax.transAxes)
        # skip over if not age's turn yet or if age not here
        if 'age' in name and age_info is True:
            if ii % 3 == 0:
                _ = display_single(
                    maps[name], ax=ax, stretch='linear', zmin=1.0, zmax=8.5,
                    color_bar=True, scale_bar=False, no_negative=True,
                    color_bar_height='5%', color_bar_width='85%', color_bar_fontsize=20,
                    cmap=IMG_CMAP, color_bar_color='k')
                _ = ax.text(0.06, 0.06, r'$\rm Age/Gyr$', fontsize=25,
                            transform=ax.transAxes,
                            bbox=dict(facecolor='w', edgecolor='none', alpha=0.8))
            else:
                _ = display_single(
                    maps[name], ax=ax, stretch='linear', zmin=1.0, zmax=8.5,
                    color_bar=False, scale_bar=False, no_negative=True,
                    cmap=IMG_CMAP, color_bar_color='k')
        # skip over if not met's turn yet or if met not here
        if 'met' in name and met_info is True:
            if ii % 3 == 0:
                _ = display_single(
                    maps[name] / Z_SUN, ax=ax, stretch='log10', zmin=-0.6, zmax=0.9,
                    color_bar=True, scale_bar=False, no_negative=True,
                    color_bar_height='5%', color_bar_width='85%', color_bar_fontsize=20,
                    cmap=IMG_CMAP, color_bar_color='k')
                _ = ax.text(0.06, 0.06, r'$\log[Z_{\star}/Z_{\odot}]$', fontsize=25,
                            transform=ax.transAxes,
                            bbox=dict(facecolor='w', edgecolor='none', alpha=0.8))
            else:
                _ = display_single(
                    maps[name] / Z_SUN, ax=ax, stretch='log10', zmin=-0.6, zmax=0.9,
                    color_bar=False, scale_bar=False, no_negative=True,
                    cmap=IMG_CMAP, color_bar_color='k')

    return fig_sum


def show_aper(info, aper, age_info=True, met_info=True, figsize=(8, 18), rad_min=5.5, rad_max=170.):
    """Make a summary plot of the aperture measurements.

    Parameters
    ----------
    info : dict
        A dictionary that contains basic information of the galaxy
    aper : dict
        A dictionary that contains the aperture measurements of stellar mass, age,
        and metallicity.
    rad_min : float, optional
        Minimum radius to plot, in unit of kpc. Default: 5.5.
    rad_max : float, optional
        Maximum radius to plot, in unit of kpc. Default: 170.
    figsize : tuple, optional
        Size of the 3x3 figure. Default: (15, 15)

    """
    # Integrated properties of the galaxy
    logms = info['logms']
    if age_info is True and met_info is True:
        age = info['age']
        logz = np.log10(info['metallicity'] / Z_SUN)

    # Radial mask
    rad_mask = aper['rad_mid'] >= rad_min

    # Setup the figure
    fig_prof = plt.figure(figsize=figsize, constrained_layout=False)
    # add rows based on if age or met are True, or, 1 (why did 2+age+met)
    grid_prof = fig_prof.add_gridspec(2+age_info+met_info, 1, wspace=0.0, hspace=0.0)
    fig_prof.subplots_adjust(
        left=0.175, right=0.93, bottom=0.055, top=0.995,
        wspace=0.00, hspace=0.00)

    # Integrated mass profile
    ax0 = fig_prof.add_subplot(grid_prof[0])
    ax0.scatter(
        aper['rad_mid'] ** 0.25, np.log10(aper['maper_gal']),
        c='darkgrey', marker='s', s=60, label=r'$\rm Total$')
    ax0.scatter(
        aper['rad_mid'] ** 0.25, np.log10(aper['maper_ins']),
        c='orangered', marker='o', alpha=0.8, s=70, label=r'$\rm In\ situ$')
    ax0.scatter(
        aper['rad_mid'] ** 0.25, np.log10(aper['maper_exs']),
        c='steelblue', marker='h', alpha=0.8, s=80, label=r'$\rm Ex\ situ$')
    ax0.axhline(logms, linewidth=2.5, linestyle='--', alpha=0.8, c='k',
                label='__no_label__')

    ax0.legend(fontsize=22, loc='best')
    ax0.grid(linestyle='--', alpha=0.5)

    _ = ax0.set_xlim(rad_min ** 0.25, rad_max ** 0.25)

    mass_arr = np.stack(
        [np.log10(aper['maper_gal'][rad_mask]),
         np.log10(aper['maper_ins'][rad_mask]),
         np.log10(aper['maper_exs'][rad_mask])])

    _ = ax0.set_ylim(np.nanmin(mass_arr) - 0.09, logms + 0.15)
    _ = ax0.set_ylabel(r'$\rm Curve\ of\ Growth$', fontsize=28)

    # Radial mass bin profile
    ax1 = fig_prof.add_subplot(grid_prof[1])
    ax1.scatter(
        aper['rad_mid'] ** 0.25, np.log10(aper['mprof_ins']),
        c='orangered', marker='o', alpha=0.8, s=70, label=r'$\rm In\ situ$')
    ax1.scatter(
        aper['rad_mid'] ** 0.25, np.log10(aper['mprof_exs']),
        c='steelblue', marker='h', alpha=0.8, s=80, label=r'$\rm Ex\ situ$')

    ax1.grid(linestyle='--', alpha=0.5)

    _ = ax1.set_xlim(rad_min ** 0.25, rad_max ** 0.25)

    # `mprof` could be 0.0 in the outskirt for some compact galaxies
    # Setup a threshold to filter out those data points.
    ins_mask = aper['mprof_ins'] >= 1e7
    exs_mask = aper['mprof_exs'] >= 1e7

    mbins_arr = np.concatenate(
        [np.log10(aper['mprof_ins'][rad_mask & ins_mask]),
         np.log10(aper['mprof_exs'][rad_mask & exs_mask])], axis=None)
    _ = ax1.set_ylim(np.nanmin(mbins_arr) * 0.95, np.nanmax(mbins_arr) * 1.05)
    _ = ax1.set_ylabel(r'$\log [M_{\star}/M_{\odot}]$', fontsize=28)

    # Ex-situ fraction
    fexs = aper['mprof_exs'] / aper['mprof_gal']

    ax1_b = fig_prof.add_axes(ax1.get_position())
    ax1_b.patch.set_visible(False)
    ax1_b.xaxis.set_visible(False)
    ax1_b.spines['right'].set_color('maroon')
    ax1_b.tick_params(axis='y', colors='maroon')
    ax1_b.yaxis.set_label_position('right')
    ax1_b.yaxis.set_ticks_position('right')
    ax1_b.plot(aper['rad_mid'][rad_mask] ** 0.25, fexs[rad_mask], linestyle='--',
               c='maroon', linewidth=3.5, alpha=0.8, label=r'$\rm Ex\ situ\ fraction$')
    ax1_b.set_ylim(0.02, 0.98)
    ax1_b.legend(fontsize=22, loc='best')

    # Metallicity profiles
    if met_info is True: 
        ax2 = fig_prof.add_subplot(grid_prof[2])
        ax2.scatter(
            aper['rad_mid'] ** 0.25, np.log10(aper['met_gal_w'] / Z_SUN),
            c='darkgrey', marker='s', s=60, label='__no_label__')
        ax2.scatter(
            aper['rad_mid'] ** 0.25, np.log10(aper['met_ins_w'] / Z_SUN),
            c='orangered', marker='o', s=70, alpha=0.8, label='__no_label__')
        ax2.scatter(
            aper['rad_mid'] ** 0.25, np.log10(aper['met_exs_w'] / Z_SUN),
            c='steelblue', marker='h', s=80, alpha=0.8, label='__no_label__')
        ax2.scatter(
            aper['rad_mid'] ** 0.25, np.log10(aper['met_ins'] / Z_SUN),
            edgecolor='orangered', marker='o', s=80, alpha=0.8, label='__no_label__',
            facecolor='none', linewidth=2)
        ax2.scatter(
            aper['rad_mid'] ** 0.25, np.log10(aper['met_exs'] / Z_SUN),
            edgecolor='steelblue', marker='h', s=90, alpha=0.8, label='__no_label__',
            facecolor='none', linewidth=2)
        ax2.axhline(logz, linewidth=2.5, linestyle='--', alpha=0.8, c='k',
                    label=r'$\rm Catalog\ value$')

        ax2.grid(linestyle='--', alpha=0.5)
        ax2.legend(fontsize=22, loc='best')

        _ = ax2.set_xlim(rad_min ** 0.25, rad_max ** 0.25)
        met_arr = np.stack(
            [np.log10(aper['met_ins_w'][rad_mask] / Z_SUN),
             np.log10(aper['met_exs_w'][rad_mask] / Z_SUN)])
        _ = ax2.set_ylim(np.nanmin(met_arr) - 0.15, np.nanmax(met_arr) + 0.09)
        _ = ax2.set_ylabel(r'$\log [Z_{\star}/Z_{\odot}]$', fontsize=28)

    # Age profiles
    if age_info is True:
        ax3 = fig_prof.add_subplot(grid_prof[3])
        ax3.scatter(
            aper['rad_mid'] ** 0.25, aper['age_gal_w'],
            c='darkgrey', marker='s', s=60, label='__no_label__')
        ax3.scatter(
            aper['rad_mid'] ** 0.25, aper['age_ins_w'],
            c='orangered', marker='o', s=70, alpha=0.8, label=r'$\rm Weighted$')
        ax3.scatter(
            aper['rad_mid'] ** 0.25, aper['age_exs_w'],
            c='steelblue', marker='h', s=80, alpha=0.8, label='__no_label__')
        ax3.scatter(
            aper['rad_mid'] ** 0.25, aper['age_ins'],
            edgecolor='orangered', marker='o', s=80, alpha=0.8, label=r'$\rm Not\ Weighted$',
            facecolor='none', linewidth=2)
        ax3.scatter(
            aper['rad_mid'] ** 0.25, aper['age_exs'],
            edgecolor='steelblue', marker='h', s=90, alpha=0.8, label='__no_label__',
            facecolor='none', linewidth=2)
        ax3.axhline(age, linewidth=2.5, linestyle='--', alpha=0.8, c='k',
                    label=r'__no_label__')

        ax3.grid(linestyle='--', alpha=0.5)
        ax3.legend(fontsize=20, loc='best')

        _ = ax3.set_xlim(rad_min ** 0.25, rad_max ** 0.25)
        age_arr = np.stack(
            [aper['age_ins_w'][rad_mask], aper['age_exs_w'][rad_mask]])
        _ = ax3.set_ylim(np.nanmin(age_arr) * 0.8, np.nanmax(age_arr) + 1.5)

        _ = ax3.set_xlabel(r'$[R/{\rm kpc}]^{1/4}$', fontsize=28)
        _ = ax3.set_ylabel(r'$[\rm Age/Gyrs]$', fontsize=28)

    return fig_prof


def prepare_show_ellipse(info, maps, ell_sum):
    """Prepare the data for visualizing the 1-D profiles.

    Parameters
    ----------
    info : dict
        A dictionary that contains basic information of the galaxy
    maps : dict
        A dictionary that contains all stellar mass, age, and metallicity maps.
    ell_sum : dict
        A dictionary that contains all the 1-D Ellipse profiles.

    Returns
    -------
    ell_plot : dict
        A dictionary that summarizes all necessary data for visualizing Ellipse profiles.

    """
    return {'catsh_id': info['catsh_id'],
            'logms': info['logms'],
            'pix': info['pix'],
            'mass_gal': maps['mass_gal'],
            'mass_ins': maps['mass_ins'],
            'mass_exs': maps['mass_exs'],
            'ell_gal_2': ell_sum['gal_shape'],
            'ell_gal_3': ell_sum['gal_mprof'],
            'ell_ins_2': ell_sum['ins_shape'],
            'ell_ins_3': ell_sum['ins_mprof'],
            'ell_exs_2': ell_sum['exs_mprof'],
            'ell_exs_3': ell_sum['exs_mprof']
           }


def overplot_ellipse(ell_plot, zmin=3.5, zmax=10.5):
    """Overplot the elliptical isophotes on the stellar mass maps.

    Parameters
    ----------
    ell_plot : dict
        A dictionary that summarizes all necessary data for visualizing Ellipse profiles.
    zmin : float, optional
        Minimum log10(Mass) value used to show the stellar mass map. Default: 3.5
    zmax : float, optional
        Maximum log10(Mass) value used to show the stellar mass map. Default: 10.5

    """
    # Setup the figure
    fig = plt.figure(figsize=(10, 10))
    fig.subplots_adjust(
        left=0.005, right=0.995, bottom=0.005, top=0.995,
        wspace=0.00, hspace=0.00)

    # Build the grid
    gs = GridSpec(2, 2)
    gs.update(wspace=0.0, hspace=0.00)

    # Central galaxy: step 2
    ax1 = fig.add_subplot(gs[0])
    ax1.yaxis.set_major_formatter(NullFormatter())
    ax1.xaxis.set_major_formatter(NullFormatter())
    ax1 = display_single(
        ell_plot['mass_gal'], ax=ax1, stretch='log10', zmin=zmin, zmax=zmax,
        cmap=IMG_CMAP, no_negative=True, color_bar=True, scale_bar=False,
        color_bar_color='k')

    if ell_plot['ell_gal_2'] is not None:
        for k, iso in enumerate(ell_plot['ell_gal_2']):
            if k % 3 == 0 and iso['sma'] >= 6.0:
                e = Ellipse(xy=(iso['x0'], iso['y0']), height=iso['sma'] * 2.0,
                            width=iso['sma'] * 2.0 * (1.0 - iso['ell']),
                            angle=iso['pa'])
                e.set_facecolor('none')
                e.set_edgecolor('k')
                e.set_alpha(0.6)
                e.set_linewidth(2.0)
                ax1.add_artist(e)
        ax1.set_aspect('equal')

    _ = ax1.text(0.05, 0.06, r'$\rm Total$', fontsize=25,
                 transform=ax1.transAxes,
                 bbox=dict(facecolor='w', edgecolor='none', alpha=0.8))

    # Central galaxy: step 3
    ax2 = fig.add_subplot(gs[1])
    ax2.yaxis.set_major_formatter(NullFormatter())
    ax2.xaxis.set_major_formatter(NullFormatter())
    ax2 = display_single(
        ell_plot['mass_gal'], ax=ax2, stretch='log10', zmin=zmin, zmax=zmax,
        cmap=IMG_CMAP, no_negative=True, color_bar=False, scale_bar=True,
        pixel_scale=1., physical_scale=ell_plot['pix'], scale_bar_loc='right',
        scale_bar_length=50., scale_bar_color='k', scale_bar_y_offset=1.3)

    ax2.text(
        0.5, 0.92,
        r'$\mathrm{ID}: %d\ \ \log M_{\star}: %5.2f$' % (
            ell_plot['catsh_id'], ell_plot['logms']),
        fontsize=21, transform=ax2.transAxes,
        horizontalalignment='center', verticalalignment='center',
        bbox=dict(facecolor='w', edgecolor='none', alpha=0.5))

    # Show the average isophotal shape
    if ell_plot['ell_ins_3'] is not None:
        n_iso = len(ell_plot['ell_ins_3'])
        if n_iso > 15:
            idx_use = n_iso - 6
        else:
            idx_use = n_iso - 1
        for k, iso in enumerate(ell_plot['ell_ins_3']):
            if k == idx_use:
                e = Ellipse(xy=(iso['x0'], iso['y0']), height=iso['sma'] * 2.0,
                            width=iso['sma'] * 2.0 * (1.0 - iso['ell']),
                            angle=iso['pa'])
                e.set_facecolor('none')
                e.set_edgecolor('k')
                e.set_linestyle('--')
                e.set_alpha(0.8)
                e.set_linewidth(2.5)
                ax2.add_artist(e)
        ax2.set_aspect('equal')

    _ = ax2.text(0.05, 0.06, r'$\rm Total$', fontsize=25,
                 transform=ax2.transAxes,
                 bbox=dict(facecolor='w', edgecolor='none', alpha=0.8))

    # In situ component: step 2
    ax3 = fig.add_subplot(gs[2])
    ax3.yaxis.set_major_formatter(NullFormatter())
    ax3.xaxis.set_major_formatter(NullFormatter())
    ax3 = display_single(
        ell_plot['mass_ins'], ax=ax3, stretch='log10', zmin=zmin, zmax=zmax,
        cmap=IMG_CMAP, no_negative=True, color_bar=False, scale_bar=False)

    if ell_plot['ell_ins_2'] is not None:
        for k, iso in enumerate(ell_plot['ell_ins_2']):
            if k % 3 == 0 and iso['sma'] >= 6.0:
                e = Ellipse(xy=(iso['x0'], iso['y0']), height=iso['sma'] * 2.0,
                            width=iso['sma'] * 2.0 * (1.0 - iso['ell']),
                            angle=iso['pa'])
                e.set_facecolor('none')
                e.set_edgecolor('orangered')
                e.set_alpha(0.9)
                e.set_linewidth(2.0)
                ax3.add_artist(e)
        ax3.set_aspect('equal')

    _ = ax3.text(0.05, 0.06, r'$\rm In\ situ$', fontsize=25,
                 transform=ax3.transAxes,
                 bbox=dict(facecolor='w', edgecolor='none', alpha=0.8))

    # Ex situ component: step 2
    ax4 = fig.add_subplot(gs[3])
    ax4.yaxis.set_major_formatter(NullFormatter())
    ax4.xaxis.set_major_formatter(NullFormatter())
    ax4 = display_single(
        ell_plot['mass_exs'], ax=ax4, stretch='log10', zmin=zmin, zmax=zmax,
        cmap=IMG_CMAP, no_negative=True, color_bar=False, scale_bar=False)

    if ell_plot['ell_exs_2'] is not None:
        for k, iso in enumerate(ell_plot['ell_exs_2']):
            if k % 3 == 0 and iso['sma'] >= 6.0:
                e = Ellipse(xy=(iso['x0'], iso['y0']), height=iso['sma'] * 2.0,
                            width=iso['sma'] * 2.0 * (1.0 - iso['ell']),
                            angle=iso['pa'])
                e.set_facecolor('none')
                e.set_edgecolor('steelblue')
                e.set_alpha(0.9)
                e.set_linewidth(2.0)
                ax4.add_artist(e)
        ax4.set_aspect('equal')

    _ = ax4.text(0.05, 0.06, r'$\rm Ex\ situ$', fontsize=25,
                 transform=ax4.transAxes,
                 bbox=dict(facecolor='w', edgecolor='none', alpha=0.8))

    return fig


def plot_ell_prof(ell_plot, r_min=3.0, r_max=190.0, insitu=False, exsitu=False):
    """Plot a summary plot for the ellipse result.

    Parameters
    ----------
    ell_plot : dict
        A dictionary that summarizes all necessary data for visualizing Ellipse profiles.
    r_min : float, optional
        Minimum radius to plot, in unit of kpc. Default: 3.0.
    r_max : float, optional
        Maximum radius to plot, in unit of kpc. Default: 190.

    """
    # Setup the figure and axes
    fig = plt.figure(figsize=(10, 10))
    fig.subplots_adjust(
        left=0.0, right=1.0, bottom=0.00, top=1.0,
        wspace=0.00, hspace=0.00)

    ax1 = fig.add_axes([0.09, 0.10, 0.90, 0.48])
    ax2 = fig.add_axes([0.09, 0.58, 0.90, 0.21])
    ax3 = fig.add_axes([0.09, 0.79, 0.90, 0.20])

    # 1-D profile
    ax1.grid(linestyle='--', alpha=0.4, linewidth=2)

    if ell_plot['ell_gal_3'] is not None:
        ax1.errorbar(
            ell_plot['ell_gal_3']['r_kpc'] ** 0.25,
            np.log10(ell_plot['ell_gal_3']['intens']),
            yerr=ell_plot['ell_gal_3']['sbp_err'], markersize=8,
            color='darkgrey', alpha=0.8, fmt='s', capsize=3,
            capthick=1, elinewidth=1, label=r'$\mathrm{Total}$')

    if insitu is not False and ell_plot['ell_ins_3'] is not None:
        ax1.errorbar(
            ell_plot['ell_ins_3']['r_kpc'] ** 0.25 + 0.05,
            np.log10(ell_plot['ell_ins_3']['intens']),
            yerr=ell_plot['ell_ins_3']['sbp_err'], markersize=9,
            color='orangered', alpha=0.8, fmt='o', capsize=3,
            capthick=1, elinewidth=1, label=r'$\mathrm{In\ Situ}$')

    if exsitu is not False and ell_plot['ell_exs_3'] is not None:
        ax1.errorbar(
            ell_plot['ell_exs_3']['r_kpc'] ** 0.25 - 0.05,
            np.log10(ell_plot['ell_exs_3']['intens']),
            yerr=ell_plot['ell_exs_3']['sbp_err'], markersize=9,
            color='steelblue', alpha=0.8, fmt='h', capsize=3,
            capthick=1, elinewidth=1, label=r'$\mathrm{Ex\ Situ}$')

    ax1.legend(loc='best', fontsize=23)

    if (ell_plot['ell_gal_3'] is not None and ell_plot['ell_ins_3'] is not None and
            ell_plot['ell_exs_3'] is not None):
        mass_arr = (
            list(ell_plot['ell_exs_3']['intens'][ell_plot['ell_exs_3']['r_kpc'] > r_min]) +
            list(ell_plot['ell_ins_3']['intens'][ell_plot['ell_ins_3']['r_kpc'] > r_min]) +
            list(ell_plot['ell_gal_3']['intens'][ell_plot['ell_gal_3']['r_kpc'] > r_min]))
        min_mass = np.nanmin(mass_arr)
        if min_mass > 100.0:
            ax1.set_ylim(np.log10(min_mass) - 0.3, np.log10(np.nanmax(mass_arr)) + 0.5)
        else:
            ax1.set_ylim(2.01, np.log10(np.nanmax(mass_arr)) + 0.5)

    ax1.set_xlim(r_min ** 0.25, r_max ** 0.25)

    _ = ax1.set_xlabel(r'$R/\mathrm{kpc}^{1/4}$', fontsize=28)
    _ = ax1.set_ylabel(r'$\log\ (\mu_{\star}/[M_{\odot}\ \mathrm{kpc}^{-2}])$', fontsize=28)

    # Ellipticity profile
    ax2.grid(linestyle='--', alpha=0.4, linewidth=2)

    if ell_plot['ell_gal_3'] is not None:
        ax2.axhline(ell_plot['ell_gal_3']['ell'][1], c='k', linestyle='--',
                    linewidth=3, alpha=0.5)

    if ell_plot['ell_gal_2'] is not None:
        ax2.errorbar(
            ell_plot['ell_gal_2']['r_kpc'] ** 0.25, ell_plot['ell_gal_2']['ell'],
            yerr=ell_plot['ell_gal_2']['ell_err'], color='darkgrey', alpha=0.7, fmt='s',
            capsize=3, capthick=1, elinewidth=2, markersize=8)

    if ell_plot['ell_ins_2'] is not None:
        ax2.errorbar(
            ell_plot['ell_ins_2']['r_kpc'] ** 0.25 + 0.05, ell_plot['ell_ins_2']['ell'],
            yerr=ell_plot['ell_ins_2']['ell_err'], color='orangered', alpha=0.7, fmt='o',
            capsize=3, capthick=1, elinewidth=1, markersize=9)

    if ell_plot['ell_exs_2'] is not None:
        ax2.errorbar(
            ell_plot['ell_exs_2']['r_kpc'] ** 0.25 - 0.05, ell_plot['ell_exs_2']['ell'],
            yerr=ell_plot['ell_exs_2']['ell_err'], color='steelblue', alpha=0.6, fmt='h',
            capsize=3, capthick=1, elinewidth=1, markersize=9)

    if ell_plot['ell_exs_2'] is not None and ell_plot['ell_ins_2'] is not None:
        ell_arr = (
            list(ell_plot['ell_exs_2']['ell'][ell_plot['ell_exs_2']['r_kpc'] > r_min]) +
            list(ell_plot['ell_ins_2']['ell'][ell_plot['ell_ins_2']['r_kpc'] > r_min]))
        ax2.set_ylim(np.nanmin(ell_arr) - 0.05, np.nanmax(ell_arr) + 0.05)

    ax2.xaxis.set_major_formatter(NullFormatter())
    ax2.set_xlim(r_min ** 0.25, r_max ** 0.25)
    _ = ax2.set_ylabel(r'$e$', fontsize=25)

    # Position Angle profile
    ax3.grid(linestyle='--', alpha=0.4, linewidth=2)

    if ell_plot['ell_gal_3'] is not None:
        ax3.axhline(ell_plot['ell_gal_3']['pa'][1], c='k', linestyle='--',
                    linewidth=3, alpha=0.5)

    if ell_plot['ell_gal_2'] is not None:
        ax3.errorbar(
            ell_plot['ell_gal_2']['r_kpc'] ** 0.25, ell_plot['ell_gal_2']['pa'],
            yerr=ell_plot['ell_gal_2']['pa_err'], color='darkgrey', alpha=0.7, fmt='s',
            capsize=3, capthick=1, elinewidth=2, markersize=8)

    if ell_plot['ell_ins_2'] is not None:
        ax3.errorbar(
            ell_plot['ell_ins_2']['r_kpc'] ** 0.25 + 0.05, ell_plot['ell_ins_2']['pa'],
            yerr=ell_plot['ell_ins_2']['pa_err'], color='orangered', alpha=0.7, fmt='o',
            capsize=3, capthick=1, elinewidth=1, markersize=9)

    if ell_plot['ell_exs_2'] is not None:
        ax3.errorbar(
            ell_plot['ell_exs_2']['r_kpc'] ** 0.25 - 0.05, ell_plot['ell_exs_2']['pa'],
            yerr=ell_plot['ell_exs_2']['pa_err'], color='steelblue', alpha=0.6, fmt='h',
            capsize=3, capthick=1, elinewidth=1, markersize=9)

    if ell_plot['ell_exs_2'] is not None and ell_plot['ell_ins_2'] is not None:
        pa_arr = (
            list(ell_plot['ell_exs_2']['pa'][ell_plot['ell_exs_2']['r_kpc'] > r_min]) +
            list(ell_plot['ell_ins_2']['pa'][ell_plot['ell_ins_2']['r_kpc'] > r_min]))
        ax3.set_ylim(np.nanmin(pa_arr) - 10.0, np.nanmax(pa_arr) + 10.0)

    ax3.xaxis.set_major_formatter(NullFormatter())
    ax3.set_xlim(r_min ** 0.25, r_max ** 0.25)
    _ = ax3.set_ylabel(r'$\mathrm{PA\ [deg]}$', fontsize=20)

    return fig


def plot_ell_fourier(fourier, pix=1.0, r_min=6.0, r_max=190.0, show_both=False):
    """Plot a summary plot for the ellipse result."""
    # Setup the figure and axes
    fig = plt.figure(figsize=(10, 10))
    fig.subplots_adjust(
        left=0.0, right=1.0, bottom=0.00, top=1.0,
        wspace=0.00, hspace=0.00)

    ax1 = fig.add_axes([0.13, 0.091, 0.865, 0.224])
    ax2 = fig.add_axes([0.13, 0.315, 0.865, 0.225])
    ax3 = fig.add_axes([0.13, 0.540, 0.865, 0.225])
    ax4 = fig.add_axes([0.13, 0.765, 0.865, 0.224])

    # A1 and/or B1
    ax1.grid(linestyle='--', alpha=0.4, linewidth=2)
    ax1.axhline(0.0, linestyle='--', alpha=0.7, linewidth=2, color='k')

    ax1.errorbar(
        (fourier['r_pix'] * pix) ** 0.25, fourier['a1'],
        yerr=fourier['a1_err'], markersize=8,
        color='coral', alpha=0.9, fmt='o', capsize=3, capthick=2, elinewidth=2,
        label=r'$\mathrm{a}$')

    if show_both:
        ax1.errorbar(
            (fourier['r_pix'] * pix) ** 0.25, fourier['b1'],
            yerr=fourier['b1_err'], markersize=8,
            color='dodgerblue', alpha=0.6, fmt='s', capsize=3, capthick=2, elinewidth=2,
            label=r'$\mathrm{b}$')

        ax1.legend(loc='best', fontsize=15)

    ax1.set_xlim(r_min ** 0.25, r_max ** 0.25)
    ax1.yaxis.set_major_formatter(FormatStrFormatter(r'$%4.1f$'))

    _ = ax1.set_xlabel(r'$R/\mathrm{kpc}^{1/4}$', fontsize=28)
    if not show_both:
        _ = ax1.set_ylabel(r'$\rm a_{1}$', fontsize=30)
    else:
        _ = ax1.set_ylabel(r'$\rm a_{1}\ or\ b_{1}$', fontsize=30)

    # A2 and/or B2
    ax2.grid(linestyle='--', alpha=0.4, linewidth=2)
    ax2.axhline(0.0, linestyle='--', alpha=0.7, linewidth=2)

    ax2.errorbar(
        (fourier['r_pix'] * pix) ** 0.25, fourier['a2'],
        yerr=fourier['a2_err'], markersize=8,
        color='coral', alpha=0.9, fmt='o', capsize=3, capthick=2, elinewidth=2,
        label=r'$\mathrm{a2}$')

    if show_both:
        ax2.errorbar(
            (fourier['r_pix'] * pix) ** 0.25, fourier['b2'],
            yerr=fourier['b2_err'], markersize=8,
            color='dodgerblue', alpha=0.6, fmt='s', capsize=3, capthick=2, elinewidth=2,
            label=r'$\mathrm{b2}$')

    ax2.xaxis.set_major_formatter(NullFormatter())
    ax2.yaxis.set_major_formatter(FormatStrFormatter(r'$%4.1f$'))
    ax2.set_xlim(r_min ** 0.25, r_max ** 0.25)

    if not show_both:
        _ = ax2.set_ylabel(r'$\rm a_{2}$', fontsize=30)
    else:
        _ = ax2.set_ylabel(r'$\rm a_{2}\ or\ b_{2}$', fontsize=30)

    # A3 and/or B3
    ax3.grid(linestyle='--', alpha=0.4, linewidth=2)
    ax3.axhline(0.0, linestyle='--', alpha=0.7, linewidth=2)

    ax3.errorbar(
        (fourier['r_pix'] * pix) ** 0.25, fourier['a3'],
        yerr=fourier['a3_err'], markersize=8,
        color='coral', alpha=0.9, fmt='o', capsize=3, capthick=2, elinewidth=2,
        label=r'$\mathrm{a}$')

    if show_both:
        ax3.errorbar(
            (fourier['r_pix'] * pix) ** 0.25, fourier['b3'],
            yerr=fourier['b3_err'], markersize=8,
            color='dodgerblue', alpha=0.6, fmt='s', capsize=3, capthick=2, elinewidth=2,
            label=r'$\mathrm{b}$')

    ax3.xaxis.set_major_formatter(NullFormatter())
    ax3.yaxis.set_major_formatter(FormatStrFormatter(r'$%4.1f$'))
    ax3.set_xlim(r_min ** 0.25, r_max ** 0.25)

    if not show_both:
        _ = ax3.set_ylabel(r'$\rm a_{3}$', fontsize=30)
    else:
        _ = ax3.set_ylabel(r'$\rm a_{3}\ or\ b_{3}$', fontsize=30)

    # A4 and/or B4
    ax4.grid(linestyle='--', alpha=0.4, linewidth=2)
    ax4.axhline(0.0, linestyle='--', alpha=0.7, linewidth=2)

    ax4.errorbar(
        (fourier['r_pix'] * pix) ** 0.25, fourier['a4'],
        yerr=fourier['a4_err'], markersize=8,
        color='coral', alpha=0.9, fmt='o', capsize=3, capthick=2, elinewidth=2,
        label=r'$\mathrm{a}$')

    if show_both:
        ax4.errorbar(
            (fourier['r_pix'] * pix) ** 0.25, fourier['b4'],
            yerr=fourier['b4_err'], markersize=8,
            color='dodgerblue', alpha=0.6, fmt='s', capsize=3, capthick=2, elinewidth=2,
            label=r'$\mathrm{b}$')

    ax4.xaxis.set_major_formatter(NullFormatter())
    ax4.yaxis.set_major_formatter(FormatStrFormatter(r'$%4.1f$'))
    ax4.set_xlim(r_min ** 0.25, r_max ** 0.25)

    if not show_both:
        _ = ax4.set_ylabel(r'$\rm a_{4}$', fontsize=30)
    else:
        _ = ax4.set_ylabel(r'$\rm a_{4}\ or\ b_{4}$', fontsize=30)

    return fig


def plot_combined_maps(ell_array, distances=None, is_primaries=None, r_min=3.0, r_max=190.0, save_to=None):
    """Plot combined maps of all centals, satellites, all. If is_primaries == None,
    then plot color depending on distance to central. If is_primaries is given, 
    plot centrals and black and satellites as red

    Parameters
    ----------
    ell_plot_1: array
        An array of ell_plots
    r_min : float, optional
        Minimum radius to plot, in unit of kpc. Default: 3.0.
    r_max : float, optional
        Maximum radius to plot, in unit of kpc. Default: 190.

    """

    # before we start, let's normalize the distance array
    if distances is not None:
        distance_norm = list()
        for i in distances:
            i_norm = (i - min(distances)) / (max(distances) - min(distances))
            distance_norm.append(i_norm)


    # Setup the figure and axes
    fig = plt.figure(figsize=(10, 10))
    fig.subplots_adjust(
        left=0.0, right=1.0, bottom=0.00, top=1.0,
        wspace=0.00, hspace=0.00)   


    ax1 = fig.add_axes([0.09, 0.09, 0.90, 0.90])

    ax1.set_ylabel(r'$\log\ (\mu_{\star}/[M_{\odot}\ \mathrm{kpc}^{-2}])$', fontsize=28)
    ax1.set_xlabel(r'$R/\mathrm{kpc}^{1/4}$', fontsize=28)

    ax1.set_xlim(r_min ** 0.25, r_max ** 0.25)

    ax1.legend(['All'])

    ax1.grid(linestyle='--', alpha=0.4, linewidth=2)

    # 1-D profiles
    if is_primaries is None:
        ax1.set_title(r'Satellite Profiles ($11.2 < \log\ (M_{\odot})$)', fontsize=28)
        
        for i in range(len(ell_array)):
            if ell_array[i]['ell_gal_3'] is not None:
                ax1.errorbar(
                    ell_array[i]['ell_gal_3']['r_kpc'] ** 0.25,
                    np.log10(ell_array[i]['ell_gal_3']['intens']),
                    yerr=ell_array[i]['ell_gal_3']['sbp_err'], markersize=5,
                    color='darkgrey', alpha=0.8, fmt='s', capsize=3,
                    capthick=1, elinewidth=1)

                #also making a scatterplot inside the errorbar for distances
                if distances is not None:
                    ax = ax1.plot(ell_array[i]['ell_gal_3']['r_kpc'] ** 0.25,
                        np.log10(ell_array[i]['ell_gal_3']['intens']), c=cm.viridis_r(distance_norm[i]))

        if distances is not None:
            distance_colorbar = fig.colorbar(cm.ScalarMappable(norm=Normalize(min(distances), max(distances)), cmap='viridis_r'), ax=ax1)
            distance_colorbar.set_label(r'Distance to Central (ckpc/h)')

    else:
        ax1.set_title(r'Satellite Profiles ($11.2 < \log\ (M_{\odot}) < 11.4$)', fontsize=28)
        
        # add legend to plot
        centrals = mpatches.Patch(color='darkgrey', label='Centrals')
        sats = mpatches.Patch(color='darkorange', label='Satellites')
        ax1.legend(handles=[centrals, sats])
        
        for i in range(len(ell_array)):
            if ell_array[i]['ell_gal_3'] is not None:
                if is_primaries[i] == True:
                    ax1.errorbar(
                        ell_array[i]['ell_gal_3']['r_kpc'] ** 0.25,
                        np.log10(ell_array[i]['ell_gal_3']['intens']),
                        yerr=ell_array[i]['ell_gal_3']['sbp_err'], markersize=5,
                        color='darkgrey', alpha=0.8, fmt='s', capsize=3,
                        capthick=1, elinewidth=1)

                    ax = ax1.plot(ell_array[i]['ell_gal_3']['r_kpc'] ** 0.25,
                        np.log10(ell_array[i]['ell_gal_3']['intens']), c='dimgrey')

                else:
                    ax1.errorbar(
                        ell_array[i]['ell_gal_3']['r_kpc'] ** 0.25,
                        np.log10(ell_array[i]['ell_gal_3']['intens']),
                        yerr=ell_array[i]['ell_gal_3']['sbp_err'], markersize=5,
                        color='darkgrey', alpha=0.8, fmt='s', capsize=3,
                        capthick=1, elinewidth=1)

                    ax = ax1.plot(ell_array[i]['ell_gal_3']['r_kpc'] ** 0.25,
                        np.log10(ell_array[i]['ell_gal_3']['intens']), c='darkorange')

    if save_to is not None:
        plt.savefig(save_to, bbox_inches='tight')
                    
    return fig

def plot_m10_vs_m100(dict_1, dict_2, dict_3, dict_4, dict_5, dict_6, dict_7, save_to=None):
    """
    Plot mass of r<10kpc vs r<100kpc of many galaxies

    Parameters
    ----------
    dict_1 ... dict_8: each dictionary is a separate redshift where keys are some differentiator (like galaxy number)
    and value is tuple of form (M10, M100, redshift)

    save_to: filename to save to in directory where run script

    Notes:
        Make sure each index of mass_array, redshift_array, galaxy_array aligns with each other. So mass_array[i]
        corresponds to redshift_array[i] corresponds to galaxy_array[i]
    """

    # setup the figure and axes
    fig = plt.figure(figsize=(10, 10))
    fig.subplots_adjust(
        left=0.0, right=1.0, bottom=0.00, top=1.0,
        wspace=0.00, hspace=0.00)   


    ax1 = fig.add_axes([0.09, 0.09, 0.90, 0.90])

    ax1.set_title(r'M10 vs M100 (11.2 < $\log\ M_{\odot}$)', fontsize=28)
    ax1.set_ylabel(r'logM10', fontsize=28)
    ax1.set_xlabel(r'logM100', fontsize=28)

    ax1.grid(linestyle='--', alpha=0.4, linewidth=2)

    # making colors work for legends
    # first find the redshift color dict, so know what color to make each point and make sure redshift_list is not a tuple
    color_list = list()
    for i in range(8):
        color_list.append(i / 8)
        
    # converting the dictionaries into scatterplot available arrays
    M100_one, M10_one, red_one = dict_convert(dict_1)
    M100_two, M10_two, red_two = dict_convert(dict_2)
    M100_three, M10_three, red_three = dict_convert(dict_3)
    M100_four, M10_four, red_four = dict_convert(dict_4)
    M100_five, M10_five, red_five = dict_convert(dict_5)
    M100_six, M10_six, red_six = dict_convert(dict_6)
    M100_seven, M10_seven, red_seven = dict_convert(dict_7)

    # now for the plotting itself
    ax1.scatter(M100_one, M10_one, c = cm.viridis_r(color_list[0]), label = 'z = ' + str(red_one))
    ax1.scatter(M100_two, M10_two, c = cm.viridis_r(color_list[1]), label = 'z = ' + str(red_two))
    ax1.scatter(M100_three, M10_three, c = cm.viridis_r(color_list[2]), label = 'z = ' + str(red_three))
    ax1.scatter(M100_four, M10_four, c = cm.viridis_r(color_list[3]), label = 'z = ' + str(red_four))
    ax1.scatter(M100_five, M10_five, c = cm.viridis_r(color_list[4]), label = 'z = ' + str(red_five))
    ax1.scatter(M100_six, M10_six, c = cm.viridis_r(color_list[5]), label = 'z = ' + str(red_six))
    ax1.scatter(M100_seven, M10_seven, c = cm.viridis_r(color_list[6]), label = 'z = ' + str(red_seven))

    # create legends
    legend = ax1.legend(fontsize='xx-small')

    if save_to is not None:
        plt.savefig(save_to, bbox_inches='tight')

def dict_convert(dictionary):
    """
    Convert a dictionary to scatterplot available arrays for plot_m10_vs_m100()

    Parameters
    ----------
    dict: dictionary with keys as galaxy number and value is tuple where first is M10 and second is M100

    """

    x_array = list()
    y_array = list()
    for i in dictionary:
        y_array.append(dictionary[i][0])
        x_array.append(dictionary[i][1])

    redshift = dictionary[list(dictionary.keys())[0]][2]

    return x_array, y_array, redshift

def sigma_plot(info, aper, figsize=(8,18), rad_min=5.5, rad_max=170.0):
    """
    Create scatterplot of sigma profiles  
    """
    print(aper)
    fig = plt.figure()
    
    fig.subplots_adjust(
        left=0.0, right=1.0, bottom=0.00, top=1.0,
        wspace=0.00, hspace=0.00)   

    ax1 = fig.add_axes([0.09, 0.09, 0.90, 0.90])
    ax1.grid(linestyle='--', alpha=0.4, linewidth=2)

    _ = ax1.set_xlim(rad_min ** 0.25, rad_max ** 0.25)
    # ylim: 0 - 220 km/s
    # _ = ax1.set_ylim(0, 220) 

    _ = ax1.set_xlabel(r'$[R/{\rm kpc}]^{1/4}$', fontsize=28)
    _ = ax1.set_ylabel(r'$V_{disp}[\rm km/s]$', fontsize=28)
    
    ax1.scatter(
        aper['rad_mid'] ** 0.25, aper['sigma_gal_w'],
        c='black', s=60, label='__no_label__')
    
    """
    ax1.scatter(
        aper['rad_mid'] ** 0.25, aper['sigma_ins_w'],
        c='orangered', marker='o', s=70, alpha=0.8, label=r'$\rm Weighted$')
    ax1.scatter(
        aper['rad_mid'] ** 0.25, aper['sigma_exs_w'],
        c='steelblue', marker='h', s=80, alpha=0.8, label='__no_label__')
    ax1.scatter(
        aper['rad_mid'] ** 0.25, aper['sigma_ins'],
        edgecolor='orangered', marker='o', s=80, alpha=0.8, label=r'$\rm Not\ Weighted$',
        facecolor='none', linewidth=2)
    ax1.scatter(
        aper['rad_mid'] ** 0.25, aper['sigma_exs'],
        edgecolor='steelblue', marker='h', s=90, alpha=0.8, label='__no_label__',
        facecolor='none', linewidth=2)
        """
    
    return fig
