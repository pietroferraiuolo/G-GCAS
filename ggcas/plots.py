"""
Created on May 2024
    -Author: P.Ferraiuolo
"""
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, Union
from scipy.stats import gaussian_kde

label_font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 18,
        }

title_font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'semibold',
        'size': 21,
        }


def scatter_2hist(x, y, kde=False, **kwargs):
    """
    Make a 2D scatter of two quantities, with the respective histogram distributions projected on each axis.

    Parameters
    ----------
    x : ArrayLike
        Data to be scattered on x-axis
    y : ArrayLike
        Data to be scattered on y-axis

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
    if 'xlabel' in kwargs:
        xlabel=kwargs['xlabel']
    else: xlabel=''
    if 'ylabel' in kwargs:
        ylabel=kwargs['ylabel']
    else: ylabel=''
    if 'title' in kwargs:
        title=kwargs['title']
    else: title=''
    if 'alpha' in kwargs:
        alpha=kwargs['alpha']
    else: alpha=0.7
    if 'colorx' in kwargs:
        colorx=kwargs['colorx']
    else: colorx='green'
    if 'colory' in kwargs:
        colory=kwargs['colory']
    else: colory='blue'
    if 'scatter_color' in kwargs:
        sc=kwargs['scatter_color']
    else: sc='black'
    if 's' in kwargs:
        s=kwargs['s']
    elif 'size' in kwargs:
        s=kwargs['size']
    else:
        s=5
    
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
        
def colorMagnitude(g, b_r, teff_gspphot):
    """
    Perform a scatter plot to create a color-magnitude diagram of the sample, using photometry and temperature information.

    Parameters
    ----------
    g : float | ArrayLike
        Gaia mean magnitude in the G band.
    b_r : fliat | ArrayLike
        Gaia color, defined and the BP mean magnitude minus the RP mean magnitude.
    t : float | ArrayLike
        Gaia computed effective surface temperature.
    """
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7,7))
    
    ax.set_facecolor((0.9,0.9,0.9))
    plt.scatter(b_r, g, c=teff_gspphot, alpha=0.8, cmap='rainbow_r')
    plt.colorbar(label=r'$T_{eff}$')
    
    plt.ylim(max(g)+0.51, min(g)-0.51)
    
    plt.xlabel(r"$G_{BP} - G_{RP}$", fontdict=label_font)
    plt.ylabel(r"$G$", fontdict=label_font)
    plt.title("Color-Magnitude Diagram", fontsize=17)
    
    plt.show()
    
    
def properMotion(pmra, pmdec):
    '''
    

    Parameters
    ----------
    pmra : TYPE
        DESCRIPTION.
    pmdec : TYPE
        DESCRIPTION.

    '''
    fig, ax = plt.subplots(figsize=(8,8))
    plt.xlabel(r'$\mu_{\alpha*}$ [deg]', fontdict=label_font)
    plt.ylabel(r'$\mu_\delta$ [deg]', fontdict=label_font)
    plt.title('Proper Motion Distribution', fontdict=title_font)
    
    ax.axis('equal')
    plt.scatter(pmra, pmdec, c='black', alpha=0.5, s=3)
    
def raDec(ra, dec):
    '''
    

    Parameters
    ----------
    ra : TYPE
        DESCRIPTION.
    dec : TYPE
        DESCRIPTION.

    '''
    fig, ax = plt.subplots(figsize=(8,8))
    plt.xlabel(r'$DEC$ [deg]', fontdict=label_font)
    plt.ylabel(r'$RA$ [deg]', fontdict=label_font)
    plt.title('Spatial Distribution', fontdict=title_font)
    
    ax.axis('equal')
    plt.scatter(ra, dec, c='black', alpha=0.4, s=3)

def histogram(data, xlabel='x', kde=False, **kwargs):
    '''
    Plots the data distribution with a histogram. The number of bins is defined as 1.5*sqrt(N). If kde is set on Tue, the kerned desdity estimation will be computed and plotted over the histogram.

    Parameters
    ----------
    data : ArrayLike
        Imput dataset for the histogram.
    kde : Boolean
        Gaussian Kernel density estimation of the histogram. The default value is False.

    Other Parameters
    ----------------
    **kwargs : Additional parameters for customizing the plot.
    
        xlabel : str
            Label of the plot's x-axis.
        alpha : float
            Transparency of the bins.
        color : str
            Color of the histogram's bins.

    Returns
    -------
    bins : ArrayLike
        DESCRIPTION.
    counts : ArrayLike
        DESCRIPTION.

    '''
    if 'xlabel' in kwargs:
        xlabel=kwargs['xlabel']
    else: xlabel=''
    if 'alpha' in kwargs:
        alpha=kwargs['alpha']
    else: alpha=1
    if 'color' in kwargs:
        color=kwargs['color']
    else: color='gray'

    n_bin = int(1.5*np.sqrt(len(data)))
    
    plt.figure(figsize=(9,8))
    
    h = plt.hist(data, bins=n_bin, color=color, alpha=alpha)
    plt.ylabel('counts')
    plt.xlabel(xlabel, fontdict=label_font)
    title = xlabel+' Distribution'
    plt.title(title, fontdict=label_font)
    
    bins = h[1][:len(h[0])]
    binwidth = bins[1] - bins[0]
    counts = h[0]

    if kde:
        kde = gaussian_kde(data)
        xk = np.linspace(min(data), max(data), 10000)
        kde_values = kde(xk)*len(data)*binwidth
        mean=np.mean(data)
        std=np.std(data)
        
        label=r"""Gaussian KDE
mean={:.2e}
$\sigma^2$={:.2e}"""
        
        plt.plot(xk, kde_values, c='r', label=label.format(mean, std))
        plt.legend(loc='best', fontsize='large')
    
    plt.show()
    
    return bins, counts
    
def scat_xhist(x, y, xerr: Optional[Union[float, np.ndarray]] = None, xlabel: str='x', ylabel: str='y'):
    """
    
    
    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    xerr : TYPE, optional
        DESCRIPTION. The default is None.
    xlabel : TYPE, optional
        DESCRIPTION. The default is 'x'.
    ylabel : TYPE, optional
        DESCRIPTION. The default is 'y'.

    Returns
    -------
    None.

    """
    nb2 = int(1.5*np.sqrt(len(x)))
    mean_x = np.mean(x)
    
    fig, (ax0, ax1) = plt.subplots(nrows=2, ncols=1, height_ratios=[1,3.5], figsize=(8.5,8), sharex=True)
    fig.subplots_adjust(hspace=0)

    if xerr is not None:
        ax1.errorbar(x, y, xerr=xerr, fmt='x', color='red', linewidth=1., markersize=3, alpha=0.8)
        err_xm = np.sqrt(sum(i*i for i in xerr)/len(x))
    else:
        ax1.scatter(x, y, c='red', alpha=0.8, s=10)
        err_xm = np.std(x)/np.sqrt(len(x))

    ax1.set_ylim(y.min()*0.8, y.max()*1.2)

    #Media scritta
    ax1.text(x.max()*0.1, y.max(), r'$<${}$>=(${:.2f}$\,\pm\,${:.2f}$)$'.format(xlabel, mean_x, err_xm), color='black', fontsize=15)

    vh = ax0.hist(x, bins=nb2, color='red', histtype='step', orientation='vertical')

    #linea
    ax1.plot( [mean_x, mean_x], [min(y)*(0.75), max(y)*(1.25)], linestyle='--', c='black', alpha=0.85)
    ax0.plot( [mean_x, mean_x], [0, vh[0].max()],  linestyle='--', c='black', alpha=0.85)

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
    
    fig = plt.suptitle(title, fontsize=25)
    plt.show()
    
    return mean_x, err_xm
