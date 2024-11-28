"""
Author(s)
---------
    - Pietro Ferraiuolo : Written in 2024

Description
-----------
Module which contains useful plots for visualizing globular cluster's kinetic data.
Refer to the individual functions documentation for more information about their use.

How to Use
----------
Just import the module

    >>> from ggcas import plots as gplt
    >>> gplt.scatter_2hist(...) # your data

"""
import os
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, Union
from ggcas._utility import osutils as osu
from ggcas.statistics import kde_estimator as _kde_estimator

label_font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 15,
        }
title_font = {'family': 'cursive',
        'style': 'italic',
        'color':  'black',
        'weight': 'semibold',
        'size': 20,
        }

default_figure_size = (6.4, 5.2)

def scatter_2hist(x, y, kde=False, kde_kind:str='gaussian', **kwargs):
    """
    Make a 2D scatter plot of two data arrays, with the respective histogram distributions
    projected on each axis. The kde option allows for Gaussian regression on the plotted data.

    Parameters
    ----------
    x : ArrayLike
        Data to be scattered on x-axis.
    y : ArrayLike
        Data to be scattered on y-axis.
    kde : bool, optional
        Option to show the Gaussian regression on the data. Default is False.

    Other Parameters
    ----------------
    **kwargs : Additional parameters for customizing the plot.
        xlabel : str
            Label of the x-axis.
        ylabel : str
            Label of the y-axis.
        title : str
            Title of the figure.
        alpha : float
            Transparency of the data points of the scatter.
        colorx : str
            Color of the histogram on x-axis.
        colory : str
            Color of the histogram on y-axis.
        scatter_color : str
            Color of the scattered dots.
        size : int or float
            Size of the scattered data points.
        figsize : tuple
            Size of the figure.

    """
    xlabel=kwargs.get('xlabel','')
    ylabel=kwargs.get('ylabel','')
    xlim=kwargs.get('xlim', None)
    ylim=kwargs.get('ylim', None)
    title=kwargs.get('title','')
    alpha=kwargs.get('alpha', 0.7)
    colorx=kwargs.get('colorx', 'green')
    colory=kwargs.get('colory', 'blue')
    sc=kwargs.get('scatter_color', 'black')
    s=osu.get_kwargs(('size', 's'), 5, kwargs)
    fsize=kwargs.get('figsize', default_figure_size)
    fig = plt.figure(figsize=fsize)
    gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                            left=0.1, right=0.9, bottom=0.1, top=0.9,
                            wspace=0.025, hspace=0.025)
    ax = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)
    # the scatter plot:
    ax.scatter(x, y, color=sc, alpha=alpha, s=s)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax_histx.set_ylabel('Counts')
    ax_histy.set_xlabel('Counts')
    ax.set_xlabel(xlabel, fontdict=label_font)
    ax.set_ylabel(ylabel, fontdict=label_font)
    bins = int(1.5*np.sqrt(len(x)))
    ax_histx.hist(x, bins=bins, color=colorx, alpha=0.6)
    ax_histy.hist(y, bins=bins, orientation='horizontal', color=colory, alpha=0.6)
    plt.suptitle(title, size=20, style='italic', family='cursive')
    if kde:
        x_kdex,x_kdey,x_c = _kde_estimator(x, kde_kind)
        y_kdex,y_kdey,y_c = _kde_estimator(y, kde_kind)
        ax_histx.plot(x_kdex, x_kdey, color=colorx, label=f"$\mu$={x_c[1]:.3f}\n$\sigma^2$={x_c[2]:.3f}")
        ax_histy.plot(y_kdey, y_kdex, color=colory, label=f"$\mu$={y_c[1]:.3f}\n$\sigma^2$={y_c[2]:.3f}")
        ax_histx.legend(loc='best', fontsize='small')
        ax_histy.legend(loc='best', fontsize='small')
    ax_histx.set_xlim(xlim)
    ax_histy.set_ylim(ylim)
    plt.show()

def colorMagnitude(g, b_r, teff_gspphot, **kwargs):
    """
    Make a scatter plot to create a color-magnitude diagram of the sample, using BP and RP photometry and temperature information.

    Parameters
    ----------
    g : float | ArrayLike
        Gaia mean magnitude in the G band.
    b_r : float | ArrayLike
        Gaia color, defined as the BP mean magnitude minus the RP mean magnitude.
    teff_gspphot : float | ArrayLike
        Gaia computed effective surface temperature.

    Other Parameters
    ----------------
    **kwargs : dict
        bgc : tuple
            A tuple with three float values, indicating the RGB gradient which define a color (placeholder for the ax.set_facecolor function).
            Aliases: 'bgc', 'bgcolor', 'background_color'.
        alpha : float
            Transparency of the scatter points.
        cmap : str
            All accepted matplotlib colormap strings.
        figsize : tuple
            Size of the figure.

    """
    bgc = osu.get_kwargs(('bgc', 'bgcolor', 'background_color'), (0.9,0.9,0.9), kwargs)
    a = kwargs.get('alpha', 0.8)
    cmap = kwargs.get('cmap', 'rainbow_r')
    fsize = kwargs.get('figsize', default_figure_size)
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=fsize)
    ax.set_facecolor(bgc)
    plt.scatter(b_r, g, c=teff_gspphot, alpha=a, cmap=cmap)
    plt.colorbar(label=r'$T_{eff}$')
    plt.ylim(max(g)+0.51, min(g)-0.51)
    plt.xlabel(r"$G_{BP} - G_{RP}$", fontdict=label_font)
    plt.ylabel(r"$G$", fontdict=label_font)
    plt.title("Color-Magnitude Diagram", fontsize=17)
    plt.show()

def properMotion(sample, **kwargs):
    """
    Make a scatter plot in the proper motion space of the sample.

    Parameters
    ----------
    sample : _Sample or dict
        The sample data containing 'pmra' and 'pmdec' fields.

    Other Parameters
    ----------------
    **kwargs : additional parameters for customizing the plot.
        color : str
            Color of the scattered data points.
        s : int or float
            Size of the scattered data points.
        alpha : float
            Transparency of the scattered data points.
        figsize : tuple
            Size of the figure.

    """
    col  = osu.get_kwargs(('color','c'), 'black', kwargs)
    size = osu.get_kwargs(('s','size'), 3, kwargs)
    alpha= kwargs.get('alpha', 0.5)
    fsize = kwargs.get('figsize', default_figure_size)
    from ._query import _Sample
    if isinstance(sample, _Sample):
        data = sample.sample
        pmra = data['pmra']
        pmdec = data['pmdec']
    else:
        pmra = sample['pmra']
        pmdec = sample['pmdec']
    fig, ax = plt.subplots(figsize=fsize)
    plt.xlabel(r'$\mu_{\alpha*}$ [deg]', fontdict=label_font)
    plt.ylabel(r'$\mu_\delta$ [deg]', fontdict=label_font)
    plt.title('Proper Motion Distribution', fontdict=title_font)
    ax.axis('equal')
    plt.scatter(pmra, pmdec, c=col, alpha=alpha, s=size)
    plt.show()

def spatial(sample, **kwargs):
    """
    Make a scatter plot in the spatial plot, that is in the Ra-Dec plane.

    Parameters
    ----------
    sample : _Sample or dict
        The sample data containing 'ra' and 'dec' fields.

    Other Parameters
    ----------------
    **kwargs : additional parameters for customizing the plot.
        color : str
            Color of the scattered data points.
        s : int or float
            Size of the scattered data points.
        alpha : float
            Transparency of the scattered data points.
        figsize : tuple
            Size of the figure.

    """
    col  = osu.get_kwargs(('color','c'), 'black', kwargs)
    fsize= kwargs.get('figsize', default_figure_size)
    size = osu.get_kwargs(('s','size'), 5, kwargs)
    alpha= kwargs.get('alpha', 0.5)
    from ._query import _Sample
    if isinstance(sample, _Sample):
        data = sample.sample
        ra = data['ra']
        dec = data['dec']
    else:
        ra = sample['ra']
        dec = sample['dec']
    fig, ax = plt.subplots(figsize=fsize)
    plt.xlabel(r'$DEC$ [deg]', fontdict=label_font)
    plt.ylabel(r'$RA$ [deg]', fontdict=label_font)
    plt.title('Spatial Distribution', fontdict=title_font)
    ax.axis('equal')
    plt.scatter(ra, dec, c=col, alpha=alpha, s=size)
    plt.show()

def histogram(data, kde=False, kde_kind:str='gaussian', out:bool=False, **kwargs):
    """
    Plots the data distribution with a histogram. The number of bins is defined as 1.5*sqrt(N).
    If kde is True, the kernel density estimation will be computed and plotted over the histogram.

    Parameters
    ----------
    data : ArrayLike
        Input dataset for the histogram.
    kde : bool, optional
        Option for the computation of the Gaussian Kernel density estimation of the histogram. The default is False.
    kde_kind : str, optional
        Kind of kernel density estimation to be computed. The default is 'gaussian'.
        Options:
            'gaussian'
            'exponential'
            'boltzmann'
            'king'

    Other Parameters
    ----------------
    **kwargs : Additional parameters for customizing the plot.
        kde_verbose : bool
            If True, the kde iteration will be printed.
        xlabel : str
            Label of the plot's x-axis.
        alpha : float
            Transparency of the bins.
        hcolor : str
            Color of the histogram's bins.
        kcolor : str
            Color of the kde curve.
        figsize : tuple
            Size of the figure.
        xlim : tuple
            Limits for the x-axis.

    Returns
    -------
    hist_result : dict
        Dictionary containing the histogram results and the distribution mean and standard deviation. The keys are 'h' and 'kde'.

    Example
    -------
    The output can be used to make other plots and computations. For example

    >>> import numpy as np
    >>> from ggcas import plots as gplt
    >>> x = np.random.randn(1000)
    >>> y = np.random.randn(1000)
    >>> b, n = gplt.histogram(x, y)

    Then, with the result, you can make for example a scatter plot

    >>> plt.scatter(x, y)

    and so on.
    """
    xlabel = kwargs.get('xlabel','')
    alpha  = kwargs.get('alpha', 1)
    hcolor = osu.get_kwargs(('hist_color','hcolor', 'hc'),'gray', kwargs)
    kcolor = osu.get_kwargs(('kde_color','kcolor', 'kc'),'gray', kwargs)
    title  = kwargs.get('title', xlabel+' Distribution')
    fsize  = kwargs.get('figsize', default_figure_size)
    verbose= osu.get_kwargs(('kde_verbose', 'verbose', 'v'), False, kwargs)
    if 'xlim' in kwargs :
        if isinstance(kwargs['xlim'], tuple):
            xlim = kwargs['xlim']
        else:
            raise TypeError("'xlim' arg must be a tuple")
    else:
        xlim = None
    n_bin = int(1.5*np.sqrt(len(data)))
    plt.figure(figsize=fsize)
    h = plt.hist(data, bins=n_bin, color=hcolor, alpha=alpha)
    plt.ylabel('counts')
    plt.xlabel(xlabel, fontdict=label_font)
    plt.title(title, fontdict=label_font)
    bins = h[1][:len(h[0])]
    counts = h[0]
    res={'h': [bins, counts]}
    if kde:
        x_kde,y_kde,coeffs = _kde_estimator(data, kde_kind, verbose=verbose)
        res['kde'] = coeffs
        label=_kde_labels(kde_kind, *coeffs)
        plt.plot(x_kde, y_kde, c=kcolor, label=label)
        plt.legend(loc='best', fontsize='medium')
    if xlim is not None:
        plt.xlim(xlim)
    plt.show()
    if out:
        return res

def scat_xhist(x, y, xerr: Optional[Union[float, np.ndarray]] = None, **kwargs):
    """
    Make a scatter plot of a quantity 'x', with its projected histogram, relative to a quantity 'y'.

    Parameters
    ----------
    x : ndarray
        Input data.
    y : ndarray
        Related data.
    xerr : int or ndarray, optional
        Error on the input data points. Can be either a single value or an array with the same size as the input data. The default is None.

    Other Parameters
    ----------------
    **kwargs : additional arguments for customizing the plot.
        xlabel : str
            Label on x axis. The default is 'x'.
        ylabel : str
            Label on y axis. The default is 'y'.
        color : str
            Color of the scattered data points.
        size : int or float
            Size of the scattered data points.
        figsize : tuple
            Size of the figure.

    Returns
    -------
    list
        List containing the mean and error of x.
    """
    # to do: add title, remove unit (if present) to mean text and automatic title
    # change the text in a legend, so that size is fixed
    xlabel = kwargs.get('xlabel', 'x')
    ylabel = kwargs.get('ylabel', 'y')
    color  = osu.get_kwargs(('c','color'), 'gray', kwargs)
    s      = osu.get_kwargs(('s','size'), 7.5, kwargs)
    fsize  = kwargs.get('figsize', default_figure_size)
    nb2    = int(1.5*np.sqrt(len(x)))
    mean_x = np.mean(x)
    fig, (ax0, ax1) = plt.subplots(nrows=2, ncols=1, height_ratios=[1,3.5], \
                                    figsize=fsize, sharex=True)
    fig.subplots_adjust(hspace=0)
    if xerr is not None:
        if isinstance(xerr, float):
            xerr = np.full(len(x), xerr)
        ax1.errorbar(x, y, xerr=xerr, fmt='x', color=color, linewidth=1.,\
                        markersize=3, alpha=0.8)
        err_xm = np.sqrt(sum(i*i for i in xerr)/len(x))
    else:
        ax1.scatter(x, y, c=color, alpha=0.8, s=s)
        err_xm = np.std(x)/np.sqrt(len(x))
    ax1.set_ylim(y.min()*0.8, y.max()*1.2)
    #Media scritta
    ax1.text(x.max()*0.1, y.max(),
                r'$<${}$>=(${:.2f}$\,\pm\,${:.2f}$)$'.format(xlabel, mean_x, err_xm),
                color='black', fontsize=13.5)
    vh = ax0.hist(x, bins=nb2, color=color, histtype='step', orientation='vertical')
    #linea
    ax1.plot([mean_x, mean_x], [min(y)*(0.75), max(y)*(1.25)], linestyle='--',\
                c='black', alpha=0.85)
    ax0.plot([mean_x, mean_x], [0, vh[0].max()],  linestyle='--', c='black',\
                alpha=0.85)
    #labels
    ax1.set_ylabel(ylabel, fontdict=label_font)
    ax1.set_xlabel(xlabel, fontdict=label_font)
    ax0.set_ylabel('Counts')
    #minor ticks
    ax0.minorticks_on()
    ax1.minorticks_on()
    ax0.tick_params(axis="both",direction="in", size=6)
    ax0.tick_params(axis='y', which="minor",direction="in", size=0)
    ax0.tick_params(axis='x', which="minor",direction="in", size=3)
    ax1.tick_params(axis="both",direction="in", size=6)
    ax1.tick_params(which="minor",direction="in", size=3)
    title = xlabel+' distribution '
    fig = plt.suptitle(title, fontsize=23.5)
    plt.show()
    return [mean_x, err_xm]

def errorbar(data, dataerr, x=None, xerr=None, **kwargs):
    """
    Scatter plot with error bars.

    Parameters
    ----------
    data : ndarray
        Data to be plotted.
    dataerr : ndarray
        Errors associated with the data.
    x : ndarray, optional
        X-axis data. The default is None.
    xerr : ndarray, optional
        Errors associated with the x-axis data. The default is None.

    Other Parameters
    ----------------
    **kwargs : Additional callbacks for matplotlib (see matplotlib.pyplot.errorbar documentation).
        fmt : str
            Scatter point shape.
        color : str
            Scatter point color.
            Aliases - 'sc' ; 'scolor' ; 'scatcol'.
        ecolor : str
            Error bar color.
            Aliases - 'ec' ; 'errc' ; 'errcolor' ; 'errorcolor'.
        markersize : float
            Scatter point size.
            Aliases - 's' ; 'ms'.
        capsize : float
            Error bar cap length.
        elinewidth : float
            Error bar thickness.
            Aliases - 'elw' ; 'errlinew'.
        barsabove : bool
            If True, the error bars will be plotted over the scatter points.
        title : str
            Title of the plot.
        xlabel : str
            Label of the x-axis.
        ylabel : str
            Label of the y-axis.
        figsize : tuple
            Size of the figure.

    """
    ec     = osu.get_kwargs(('ecolor', 'ec', 'errc', 'errcolor','error_color'), 'red', kwargs)
    sc     = osu.get_kwargs(('color', 'scolor', 'sc', 'scatter_col'), 'black', kwargs)
    elw    = osu.get_kwargs(('elinewidth', 'elw', 'errlinew'), 1, kwargs)
    ms     = osu.get_kwargs(('markersize', 'ms', 's'), 2, kwargs)
    fsize  = kwargs.get('figsize', default_figure_size)
    ba     = kwargs.get('barsabove', False)
    cs     = kwargs.get('capsize', 1.5)
    xlabel = kwargs.get('xlabel', '')
    ylabel = kwargs.get('ylabel', '')
    title  = kwargs.get('title', '')
    fmt    = kwargs.get('fmt', 'x')
    x = np.linspace(0, 1, len(data)) if x is None else x
    plt.figure(figsize=fsize)
    plt.errorbar(x, data, yerr=dataerr, xerr=xerr, fmt=fmt, capsize=cs, ecolor=ec, \
                elinewidth=elw, markersize=ms, color=sc, barsabove=ba)
    plt.xlabel(xlabel, fontdict=label_font)
    plt.ylabel(ylabel, fontdict=label_font)
    plt.title(title, fontdict=title_font)
    plt.show()


def _kde_labels(kind: str, *coeffs):
    """
    Return the labels for the KDE plot.
    """
    if kind == 'gaussian':
        A, mu, sigma2 = coeffs
        label = f"""Gaussian KDE
$A$   = {_format_number(A)}
$\mu$   = {_format_number(mu)}
$\sigma^2$  = {_format_number(sigma2)}"""
    elif kind == 'exponential':
        A, lmbda = coeffs
        label = f"""Exponential KDE
$A$   = {_format_number(A)}
$\lambda$ = {_format_number(lmbda)}"""
    elif kind == 'boltzmann':
        A1, A2, x0, dx = coeffs
        label = f"""Boltzmann KDE
$A1$   = {_format_number(A1)}
$A2$   = {_format_number(A2)}
$x_0$   = {_format_number(x0)}
$dx$   = {_format_number(dx)}"""
    elif kind == 'king':
        A, ve, sigma = coeffs
        label = f"""King KDE
$A$   = {_format_number(A)}
$v_e$   = {_format_number(ve)}
$\sigma$  = {_format_number(sigma)}"""
    elif kind == 'maxwell':
        A, sigma = coeffs
        label = f"""Maxwell KDE
$A$   = {_format_number(A)}
$\sigma$  = {_format_number(sigma)}"""
    return label

def _format_number(num):
    """
    Format the number using scientific notation if it is too large or too small.
    """
    if abs(num) < 1e-3 or abs(num) > 1e3:
        return f"{num:.2e}"
    else:
        return f"{num:.2f}"