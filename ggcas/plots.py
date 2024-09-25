"""
Author(s)
---------
    - Pietro Ferraiuolo : Written in 2024

Description
-----------
Module which contains useful plots for visualizing globular cluster's kinetic d
ata.
Refer to the individual functions documentation for more information about thei
r use.

How to Use
----------
Just import the module

    >>> from ggcas import plots as gplt
    >>> gplt.scatter_2hist(...) # your data

"""
from typing import Optional, Union
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from ggcas._utility import osutils as osu

label_font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }
title_font = {'family': 'coursive',
        'style': 'italic',
        'color':  'black',
        'weight': 'semibold',
        'size': 20,
        }

def scatter_2hist(x, y, kde=False, **kwargs):
    """
    Make a 2D scatter of two data arrays, with the respective histogram distrib
    utions projected on each axis. The kde option allows for gaussian regresssi
    on on the plotted data.

    Parameters
    ----------
    x : ArrayLike
        Data to be scattered on x-axis
    y : ArrayLike
        Data to be scattered on y-axis
    kde : boolean, optional
        Option to show the Gaussian regression on the data. Default is false.

    Other Parameters
    ----------------
    **kwargs : Additional parameters for customizing the plot
        xlabel : str
            xlabel of the plot
        ylabel : str
            ylabel of the plot
        title : str
            title of the figure
        alpha : float
            transparency of the data points of the scatter
        colorx : str
            color of the histogram on x-axis
        colory : str
            color of the histogram on y-axis
        scatter_color : str
            color of the scattered dots
    """
    xlabel=kwargs.get('xlabel','')
    ylabel=kwargs.get('ylabel','')
    title=kwargs.get('title','')
    alpha=kwargs.get('alpha', 0.7)
    colorx=kwargs.get('colorx', 'green')
    colory=kwargs.get('colory', 'blue')
    sc=kwargs.get('scatter_color', 'black')
    s=kwargs.get('size', 5)
    fig = plt.figure(figsize=(8,8))
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
    ax_histx.set_ylabel('Counts')
    ax_histy.set_xlabel('Counts')
    ax.set_xlabel(xlabel, fontdict=label_font)
    ax.set_ylabel(ylabel, fontdict=label_font)
    binwidth = 0.25
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    lim = (int(xymax/binwidth) + 1) * binwidth
    bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_histx.hist(x, bins=bins, color=colorx, alpha=0.6)
    ax_histy.hist(y, bins=bins, orientation='horizontal', color=colory, alpha=0.6)
    plt.suptitle(title, size=21, weight='semibold')
    if kde:
        kdex = gaussian_kde(x)
        xk = np.linspace(min(x), max(x), 10000)
        kdex_values = kdex(xk)*len(x)*binwidth
        kdey = gaussian_kde(y)
        yk = np.linspace(min(y), max(y), 10000)
        kdey_values = kdey(yk)*len(y)*binwidth
        ax_histx.plot(xk, kdex_values, color='r')
        ax_histy.plot(kdey_values, yk, color='r')
    plt.show()

def colorMagnitude(g, b_r, teff_gspphot, **kwargs):
    """
    Make a scatter plot to create a color-magnitude diagram of the sample, usin
    g BP and RP photometry and temperature information.

    Parameters
    ----------
    g : float | ArrayLike
        Gaia mean magnitude in the G band.
    b_r : fliat | ArrayLike
        Gaia color, defined and the BP mean magnitude minus the RP mean magnitude.
    t : float | ArrayLike
        Gaia computed effective surface temperature.

    Other Parameters
    ----------------
    **kwargs : dict
        bgc : tuple
            A tuple with three float values, indicating the RGB gradient which
            define a color (placeholder for the ax.set_facecolor function).
            Aliases: 'bgc', 'bgcolor'
        alpha : float
            Transparency of the scatter points.
        cmap : str
            All accepted matplotlib colormap strings.
    """
    bgc = osu.get_kwargs(('bgc', 'bgcolor'), (0.9,0.9,0.9), kwargs)
    a = kwargs.get('alpha', 0.8)
    cmap = kwargs.get('cmap', 'rainbow_r')
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7,7))
    ax.set_facecolor(bgc)
    plt.scatter(b_r, g, c=teff_gspphot, alpha=a, cmap=cmap)
    plt.colorbar(label=r'$T_{eff}$')
    plt.ylim(max(g)+0.51, min(g)-0.51)
    plt.xlabel(r"$G_{BP} - G_{RP}$", fontdict=label_font)
    plt.ylabel(r"$G$", fontdict=label_font)
    plt.title("Color-Magnitude Diagram", fontsize=17)
    plt.show()

def properMotion(pmra, pmdec, **kwargs):
    """
    Make a scatter plot in the proper motion space of the sample.

    Parameters
    ----------
    pmra : ArrayLike
        Right Ascension direction of the proper mition.
    pmdec : ArrayLike
        Declination direction for the proper motion.

    Other Parameters
    ----------------
    **kwargs : additional parameters for customizing the plot.
        color : str
            Color of the scattered data points.
        s : int or float
            Size of the scattered data points.
        alpha : float
            Transparency of the scattered data points.
    """
    col  = osu.get_kwargs(('color','c'), 'black', kwargs)
    size = osu.get_kwargs(('s','size'), 3, kwargs)
    alpha= kwargs.get('alpha', 0.5)
    fig, ax = plt.subplots(figsize=(8,8))
    plt.xlabel(r'$\mu_{\alpha*}$ [deg]', fontdict=label_font)
    plt.ylabel(r'$\mu_\delta$ [deg]', fontdict=label_font)
    plt.title('Proper Motion Distribution', fontdict=title_font)
    ax.axis('equal')
    plt.scatter(pmra, pmdec, c=col, alpha=alpha, s=size)

def spatial(ra, dec, **kwargs):
    """
    Make a scatter plot in the spatial plot, that is in the Ra-Dec plane.

    Parameters
    ----------
    ra : TYPE
        DESCRIPTION.
    dec : TYPE
        DESCRIPTION.

    Other Parameters
    ----------------
    **kwargs : additional parameters for customizing the plot
        color : str
            Color of the scattered data points.
        s : int or float
            Size of the scattered data points.
        alpha : float
            Transparency of the scattered data points.
    """
    col  = osu.get_kwargs(('color','c'), 'black', kwargs)
    size = osu.get_kwargs(('s','size'), 5, kwargs)
    alpha= kwargs.get('alpha', 0.5)
    fig, ax = plt.subplots(figsize=(8,8))
    plt.xlabel(r'$DEC$ [deg]', fontdict=label_font)
    plt.ylabel(r'$RA$ [deg]', fontdict=label_font)
    plt.title('Spatial Distribution', fontdict=title_font)
    ax.axis('equal')
    plt.scatter(ra, dec, c=col, alpha=alpha, s=size)

def histogram(data, kde=False, **kwargs):
    """
    Plots the data distribution with a histogram. The number of bins is defined
    as 1.5*sqrt(N). If kde is True, the kernel density estimation will be compu
    ted and plotted over the histogram.

    Parameters
    ----------
    data : ArrayLike
        Input dataset for the histogram.
    kde : boolean
        Option for the computation of the Gaussian Kernel density estimation of
        the histogram. The default is False.

    Other Parameters
    ----------------
    **kwargs : Additional parameters for customizing the plot.
        xlabel : str
            Label of the plot's x-axis.
        alpha : float
            Transparency of the bins.
        hcolor : str
            Color of the histogram's bins.
        kcolor : str
            Color of the kde curve.

    Returns
    -------
    hist_result : dict
        Dictionary containing the hinstogram results and the distribution mean
        and standard deviation. The keys are 'h' and 'kde'.

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
    hcolor = kwargs.get('hcolor','gray')
    kcolor = kwargs.get('kcolor', 'red')
    if 'xlim' in kwargs :
        if isinstance(kwargs['xlim'], tuple):
            xlim = kwargs['xlim']
        else:
            raise TypeError("'xlim' arg must be a tuple")
    else:
        xlim = None
    n_bin = int(1.5*np.sqrt(len(data)))
    plt.figure(figsize=(9,8))
    h = plt.hist(data, bins=n_bin, color=hcolor, alpha=alpha)
    plt.ylabel('counts')
    plt.xlabel(xlabel, fontdict=label_font)
    title = xlabel+' Distribution'
    plt.title(title, fontdict=label_font)
    bins = h[1][:len(h[0])]
    binwidth = bins[1] - bins[0]
    counts = h[0]
    mean=np.mean(data)
    std=np.std(data)
    res={'h': [bins, counts]}
    res['kde'] = [mean, std]
    if kde:
        kde = gaussian_kde(data)
        xk = np.linspace(min(data), max(data), 10000)
        kde_values = kde(xk)*len(data)*binwidth
        label=r"""Gaussian KDE
$\mu$   = {:.2e}
$\sigma^2$  = {:.2e}"""
        plt.plot(xk, kde_values, c=kcolor, label=label.format(mean, std))
        plt.legend(loc='best', fontsize='large')
    if xlim is not None:
        plt.xlim(xlim)
    plt.show()
    return res

def scat_xhist(x, y, xerr: Optional[Union[float, np.ndarray]] = None, **kwargs):
    """
    Make a scatter plot of a quantity 'x', with its projected histogram, relative
    to a quantity 'y'.

    Parameters
    ----------
    x : ndarray
        Input data.
    y : ndarray
        Related data.
    xerr : int of ndarray, optional
        Error on the input data points. Can be either a single value or an arra
        y with the same size as the input data. The default is None.

    Other Parameters
    ----------------
    **kwarg : additional arguments for customizing the plot.
        xlabel : str
            Label on x axis. The default is 'x'.
        ylabel : str
            Label on y axis. The default is 'y'.
        color : str
            Color of the scattered data points.
    """
    xlabel = kwargs.get('xlabel', 'x')
    ylabel = kwargs.get('ylabel', 'y')
    color = osu.get_kwargs(('c','color'), 'gray', kwargs)
    s = osu.get_kwargs(('s','size'), 7.5, kwargs)
    nb2 = int(1.5*np.sqrt(len(x)))
    mean_x = np.mean(x)
    fig, (ax0, ax1) = plt.subplots(nrows=2, ncols=1, height_ratios=[1,3.5], \
                                   figsize=(8.5,8), sharex=True)
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
    Scatter plot with errorbars.

    Parameters
    ----------
    data : ndarray
        DESCRIPTION.
    dataerr : ndarray
        DESCRIPTION.
    x : ndarray, optional
        DESCRIPTION. The default is None.
    xerr : ndarray, optional
        DESCRIPTION. The default is None.

    Other Parameters
    ----------------
    **kwargs : Additional callbacks for matplotlib (see matplotlib.pyplot.errorbar documentation).
        fmt : str
            Scatter point shape.
        color : str
            Scatter point color.
            Aliases - 'sc' ; 'scolor' ; 'scatcol'
        ecolor : str
            Error bar color.
            Aliases - 'ec' ; 'errc' ; 'errcolor' ; 'errorcolor'
        markersize : float
            Scatter point size.
            Aliases - 's' ; 'ms'
        capsize : float
            Errorbar cap length.
        elinewidth : float
            Errorbar thickness.
            Aliases - 'elw' ; 'errlinew'
        barsabove : bool
            If Ture, the errorbars will be plotted over the scatter points
        title : str
        xlabel : str
        ylabel : str
    """
    ec     = osu.get_kwargs(('ecolor', 'ec', 'errc', 'errcolor','errorcolor'), 'red', kwargs)
    sc     = osu.get_kwargs(('color', 'scolor', 'sc', 'scatcol'), 'black', kwargs)
    ms     = osu.get_kwargs(('markersize', 'ms', 's'), 2, kwargs)
    cs     = kwargs.get('capsize', 1.5)
    ba     = kwargs.get('barsabove', False)
    elw    = osu.get_kwargs(('elinewidth', 'elw', 'errlinew'), 1, kwargs)
    fmt    = kwargs.get('fmt', 'x')
    title  = kwargs.get('title', '')
    xlabel = kwargs.get('xlabel', '')
    ylabel = kwargs.get('ylabel', '')
    x = np.linspace(0, 1, len(data)) if x is None else x
    plt.figure()
    plt.errorbar(x, data, yerr=dataerr, xerr=xerr, fmt=fmt, capsize=cs, ecolor=ec, \
                 elinewidth=elw, markersize=ms, color=sc, barsabove=ba)
    plt.xlabel(xlabel, fontdict=label_font)
    plt.ylabel(ylabel, fontdict=label_font)
    plt.title(title, fontdict=title_font)
    plt.show()
