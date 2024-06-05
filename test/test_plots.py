"""
Created on May 2024
    -Author: P.Ferraiuolo
"""
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, Union

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


def scatter_2hist(x, y):
    '''
    

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    show : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    '''
    fig = plt.figure(figsize=(8,8))

    gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)

    ax = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y, color='black', alpha=0.8)
    ax_histx.set_ylabel('Counts')
    ax_histy.set_xlabel('Counts')
   
    binwidth = 0.25
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    lim = (int(xymax/binwidth) + 1) * binwidth

    bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_histx.hist(x, bins=bins, color='blue', alpha=0.6)
    ax_histy.hist(y, bins=bins, orientation='horizontal', color='green', alpha=0.6)

    plt.show()
        
def colorMagnitude(g, b_r, teff_gspphot):
    '''
    

    Parameters
    ----------
    g : TYPE
        DESCRIPTION.
    b_r : TYPE
        DESCRIPTION.
    t : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
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

    Returns
    -------
    None.

    '''
    fig, ax = plt.subplots(figsize=(8,8))
    plt.xlabel(r'$\mu_{\alpha*}$ [deg]', fontdict=label_font)
    plt.ylabel(r'$\mu_\delta$ [deg]', fontdict=label_font)
    plt.title('Proper Motion Distribution', fontdict=title_font)
    
    ax.axis('equal')
    plt.scatter(pmra, pmdec, c='black', alpha=0.4, s=3)
    
def raDec(ra, dec):
    '''
    

    Parameters
    ----------
    ra : TYPE
        DESCRIPTION.
    dec : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    fig, ax = plt.subplots(figsize=(8,8))
    plt.xlabel(r'$DEC$ [deg]', fontdict=label_font)
    plt.ylabel(r'$RA$ [deg]', fontdict=label_font)
    plt.title('Spatial Distribution', fontdict=title_font)
    
    ax.axis('equal')
    plt.scatter(ra, dec, c='black', alpha=0.4, s=3)

def histogram(data, xlabel='x'):
    '''
    

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    xlabel : str
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.

    '''
    n_bin = int(1.5*np.sqrt(len(data)))
    
    plt.figure(figsize=(9,8))
    
    h = plt.hist(data, bins=n_bin, color='black', alpha=0.85)
    plt.ylabel('counts')
    plt.xlabel(xlabel, fontdict=label_font)
    title = xlabel+' Distribution'
    plt.title(title, fontdict=label_font)
    plt.show()
    
    return h[1][:len(h[0])], h[0]
    
def scat_xhist(x, y, xerr: Optional[Union[float, np.ndarray]] = None, xlabel: str='x', ylabel: str='y'):
    '''
    

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

    '''
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
