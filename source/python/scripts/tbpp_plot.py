# ============================================================================
# Copyright (c) 2016-2017 Giacomo Resta
#
# This file is part of TightBinding++.
#
# TightBinding++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# TightBinding++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ============================================================================

import matplotlib.pyplot as plt
import numpy as np
import matplotlib

__all__ = ['plot_kpath','plot_model','narray_read','cimshow']

def plot_kpath(kp, Ef=None, save=None, close=False, over=False):
    """Plot the results of the tbpp::KPath solver

    Parameters
    ----------
    kp : tbpp.KPath
        The tbpp::KPath object to plot
    E_f : float
        Fermi Energy
    save : str
        If not None will save figure with the above name
    close : bool
        Whether to close the figure after saving
    over : bool
        Whether to create a new figure or plot on active figure

    Returns
    ----------
    None
    """
    kp.make_ready()

    if not over:
        plt.figure()

    for i in range(kp.eigval.shape[1]):
        plt.plot(kp.eigval[:,i], 'k-')
    plt.xlabel('k')
    plt.ylabel('E [Ha]')

    label_k = kp.steps*np.arange(0, kp.kpoints.shape[0])
    plt.xticks(label_k, [x.decode('utf-8') for x in kp.kpoints['label']])

    for k in label_k:
        plt.axvline(k, -10, 10, c='0.6', ls='--')

    if Ef != None:
        plt.axhline(Ef, -10, 10, c='k', ls='--')

    if save:
        plt.savefig(save, bbox_inches='tight')
        if close:
            plt.close()


def plot_model(m, axis='xy', trans=True, hops=False,
        label='index', colors=None, legend=True,
        save=None, close=False, over=False):
    """Plot a tbpp::Model object

    Parameters
    ----------
    m : tbpp.Model
        The Model object to plot
    axis : str
        Must be a string with two character which are either 'x','y','z'
    trans : bool
        Whether to plot translational symmetries
    hops : bool
        Whether to plot hoping terms
    label : str
        Can either be 'index' to label site index, 'kind' to label site kind
        or 'none' to disable labeling
    colors : dict or None
        Dictionary with site kind as keys and colors as values.
    legend : bool
        Whether to show the legend
    save : str
        If not None will save figure with the above name
    close : bool
        Whether to close the figure after saving
    over : bool
        Whether to create a new figure or plot on active figure


    Returns
    ----------
    None
    """
    m.make_ready()

    # determine sx,sy,ix,iy
    sx = axis[0]
    sy = axis[1]
    s2i = {'x':0,'y':1,'z':2}
    if sx == sy:
        raise ValueError("Incorrect value for axis")
    if sx not in s2i or sy not in s2i:
        raise ValueError("Incorrect value for axis")
    ix = s2i[sx]
    iy = s2i[sy]

    if not over:
        plt.figure()

    # Translational coordinate axis
    if trans:
        for ti in range(m.kdim):
            plt.annotate("", xytext=(0,0), xy=(m.a[ti,ix], m.a[ti,iy]),
                    arrowprops=dict(facecolor='k', edgecolor='k',arrowstyle="->"))

    # Range for translational cells
    n_range = [0,0,0]
    for i in range(m.kdim):
        if hops and m.R.shape[0] > 0:
            n_range[i] = max(1, max(m.R[:,i]))
        else:
            n_range[i] = 1

    # Colors
    if colors == None:
        default_colors=['r','b','g','c','m','y','k'] + list(matplotlib.colors.cnames.keys())
        site_types = set([x.decode('utf-8') for x in m.site_info["kind"]])
        colors = dict(zip(site_types, default_colors[:len(site_types)]))

    # function to draw a site
    def draw_site(x,y, kind, alpha=1.0, size=100):
        c = 'red'
        if kind in colors:
            c = colors[kind]
        plt.scatter(x,y, s=size, c=c, alpha=alpha)
        if label == 'index':
            plt.annotate(str(si), xy=(x,y), xytext=(5,5),
                    textcoords='offset points', fontsize=12, alpha=alpha)
        elif label == 'kind':
            plt.annotate(m.site_info[si]['kind'].decode('utf-8'), xy=(x,y), xytext=(5,5),
                    textcoords='offset points', fontsize=12, alpha=alpha)
        elif label == 'none':
            pass
        else:
            raise ValueError("Invalid value for label")

    # sites
    for si in range(m.sites):
        x = m.site_info[si][sx]
        y = m.site_info[si][sy]
        kind = m.site_info[si]['kind'].decode('utf-8')
        draw_site(x,y,kind)
        if trans:
            # draw translational copies of home cell
            for i in range(-n_range[0], n_range[0]+1):
                for j in range(-n_range[1], n_range[1]+1):
                    for k in range(-n_range[2], n_range[2]+1):
                        if i == j == k == 0:
                            # skip home cell (dR = {0,0,0})
                            continue
                        x = m.site_info[si][sx] + i*m.a[0,ix] + j*m.a[1,ix] + k*m.a[2,ix]
                        y = m.site_info[si][sy] + i*m.a[0,iy] + j*m.a[1,iy] + k*m.a[2,iy]
                        draw_site(x, y, kind, alpha=0.2)

    # function to draw a hop
    def draw_hop(x1,y1,x2,y2,color):
        plt.annotate("", xytext=(x1, y1), xy=(x2, y2),
                arrowprops=dict(facecolor=color, alpha=1.0,
                    edgecolor=color, arrowstyle="->",
                    connectionstyle="arc3,rad=.2"))

    # hops
    if hops:
        for s1 in range(m.sites):
            for s2 in range(m.sites):
                s1s = int(m.site_info[s1]["si"])
                s2s = int(m.site_info[s2]["si"])
                s1e = int(m.site_info[s1]["ei"]+1)
                s2e = int(m.site_info[s2]["ei"]+1)
                kind = m.site_info[s1]["kind"].decode('utf-8')
                color = 'r'
                if kind in colors:
                    color = colors[kind]

                if s1 != s2:
                    # hops from s1 to s2 in home unit cell
                    if np.sum(np.abs(m.V[s2s:s2e,s1s:s1e])) != 0.0:
                        draw_hop(m.site_info[s1][sx], m.site_info[s1][sy],
                                 m.site_info[s2][sx], m.site_info[s2][sy],
                                 color)

                # hops from s1 to s2 where s2 is outside home unit cell
                for ri in range(m.hops):
                    if np.sum(np.abs(m.T[ri,s2s:s2e,s1s:s1e])) != 0.0:
                        i,j,k = m.R[ri,:]
                        x = m.site_info[s2][sx] + i*m.a[0,ix] + j*m.a[1,ix] + k*m.a[2,ix]
                        y = m.site_info[s2][sy] + i*m.a[0,iy] + j*m.a[1,iy] + k*m.a[2,iy]
                        draw_hop(m.site_info[s1][sx], m.site_info[s1][sy], x, y, color)

    # draw legend
    if legend:
        handles = []
        types = []
        for k,v in colors.items():
            handles.append(plt.Line2D([0],[0], ls='', marker='o', c=v))
            types.append(k)
        plt.legend(handles, types, numpoints=1)

    plt.xlabel(sx + ' [$a_0$]')
    plt.ylabel(sy + ' [$a_0$]')
    plt.axis('equal')

    if save:
        plt.savefig(save, bbox_inches='tight')
        if close:
            plt.close()


def cimshow(m, *args, **kwargs):
    """Plot a complex matrix

    Parameters
    ----------
    m : array_like, shape (n, m)
    *args : all args are passed to matplotlib `imshow` function
    **kwargs : all kwargs are passed to matplotlib `imshow` function

    Returns
    ----------
    figure : matplotlib figure
    """
    fig = plt.figure()
    plt.subplot(121)
    plt.imshow(np.real(m), *args, **kwargs)
    plt.title('$\Re$')
    plt.colorbar()
    plt.subplot(122)
    plt.title('$\Im$')
    plt.imshow(np.imag(m), *args, **kwargs)
    plt.colorbar()
    return fig


def narray_read(filename):
    """Load TBPP NArray dense file format

    Parameters
    ----------
    filename : str
        Path to NArray file to load

    Returns
    ----------
    data : ndarray
        A NumPy ndarray with the dimension, size and data of the saved NArray
    """
    f = open(filename, 'r')
    lines = f.readlines()
    if lines[0] != 'NArray dense\n':
        raise RuntimeError("Not NArray format")
    dim = int(lines[1][3:])
    sizes = [int(x) for x in lines[2][5:].strip().split(" ")]

    elem = 1
    for i in sizes:
        elem *= i

    if '(' in lines[3]:
        dtype = complex
    else:
        dtype = float

    data = np.ndarray(elem, dtype=dtype)
    for i in range(elem):
        if dtype is complex:
            re,im = lines[i+3].strip().split(',')
            data[i] = complex(float(re[1:]), float(im[:-1]))
        else:
            data[i] = float(lines[i+3])
    data.resize(sizes)
    return data

