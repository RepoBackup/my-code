# -*- coding: utf-8 -*-

import numpy as np

from pymatgen.electronic_structure.core import Spin, OrbitalType

from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga

from pymatgen.io.vasp.outputs import BSVasprun, Vasprun

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter
from matplotlib.collections import LineCollection

from scipy.interpolate import splev, splrep

########################################## FUNCTIONS ##############################################

class My_Axes_Y(mpl.axes.Axes):
    name = 'My_Axes_Y'

    def drag_pan(self, button, key, x, y):
        mpl.axes.Axes.drag_pan(self, button, '', x, y)


class My_Axes_X(mpl.axes.Axes):
    name = 'My_Axes_X'

    def drag_pan(self, button, key, x, y):
        mpl.axes.Axes.drag_pan(self, button, '', x, y)


def plot(locs, parse=False, figsize=None, dpi=300, lw=1.0, save=False):

    for p in (My_Axes_X, My_Axes_Y): mpl.projections.register_projection(p)

    def setDos(axis, minor, spin=[Spin.up], lim=None):
        """ Set up plotting area for DOS. If bands are plotted, set y-axis to vertical"""

        def spines(zero,hide):
            for z in zero: axD.spines[z].set_position('zero')
            for h in hide: axD.spines[h].set_visible(False)

        if len(spin) == 1:
            if spin[0].value == 1:
                axD.set_title(u'DOS $\u2191$', **title)
                if axis == 'y':
                    axD.set_xlim(0, lim[1])
                    spines(zero=['left', 'top'], hide=['right', 'bottom'])
                    axD.yaxis.tick_left()
                else:
                    axD.set_ylim(0, lim[1])
                    spines(zero=['right', 'bottom'], hide=['left', 'top'])
                    axD.xaxis.tick_bottom()
            else:
                axD.set_title(u'DOS $\u2193$', **title)
                if axis == 'y':
                    axD.set_xlim(lim[0], 0)
                    spines(zero=['right', 'top'], hide=['left', 'bottom'])
                    axD.yaxis.tick_right()
                else:
                    axD.set_ylim(lim[0], 0)
                    spines(zero=['right', 'top'], hide=['left', 'bottom'])
                    axD.xaxis.tick_top()
                    axD.xaxis.set_label_coords(0.5, -0.1)
        else:
            axD.set_title('DOS', **title)
            if axis == 'y':
                axD.set_xlim(lim)
                spines(zero=['left' ,'top'], hide=['right', 'bottom'])
                axD.set_xlabel(u'$\u2193$     $\u2191$', **arrows)
            else:
                axD.set_ylim(lim)
                spines(zero=['right', 'bottom'], hide=['left','top'])
                axD.set_ylabel(u'$\u2190$     $\u2192$', **arrows)
            axD.xaxis.set_label_coords(0.5, -0.04)

        # define tick labels for density (in arbitrary units)
        if 'DOSticks' in dosopt:
            func = lambda x, pos: "" \
                if np.isclose(x,0) else "%d" % x if np.equal(np.mod(x, 1), 0) else x
        else:
            func = lambda x, pos: ""

        if axis == 'y': # if bands are plotted...
            axD.xaxis.tick_top()
            axD.xaxis.set_major_formatter(FuncFormatter(func))
            axD.yaxis.set_minor_locator(MultipleLocator(minor))
            plt.subplots_adjust(
                    top=0.92,
                    bottom=0.08,
                    left=0.12,
                    right=0.96,
                    hspace=0.2,
                    wspace=0.2
                    )
        else:
            axD.yaxis.tick_right()
            axD.set_xlabel("E - E$\mathrm{_f}$ (eV)", **axes)
            axD.yaxis.set_major_formatter(FuncFormatter(func))
            axD.xaxis.set_minor_locator(MultipleLocator(minor))
            plt.tight_layout()
            plt.subplots_adjust(
                    top=0.85,
                    bottom=0.21,
                    left=0.08,
                    right=0.91,
                    hspace=0.2,
                    wspace=0.2
                    )

        global _axes; _axes = [axD]
        if 'zoom' in dosopt and dosopt['zoom']:
            if len(dosopt['zoom']) > 2:
                axz = axD.inset_axes(dosopt['zoom'][2])
            elif axis == 'y':
                axz = axD.inset_axes([0.5, 0., 0.5, 0.25])
            else:
                axz = axD.inset_axes([0., 0.5, 0.25, 0.6])
            axz.set_xlim(dosopt['zoom'][0])
            axz.set_ylim(dosopt['zoom'][1])
            axz.set_xticklabels('')
            axz.set_yticklabels('')
            axz.set_label('zoom')
            _axes.append(axz)
            axD.indicate_inset_zoom(axz)

    ###########################################################################

    # set general linewidth for plot borders
    for i in ['axes.linewidth', 'xtick.major.width', 'xtick.minor.width',
              'ytick.major.width', 'ytick.minor.width']:
        mpl.rcParams[i] = lw

    # define plot axes
    global bandruns, dosruns, bands
    
    if parse: 
        bandruns = []; dosruns = []
    
    for i, loc in enumerate(locs):
        
#        proj = bandopt['projections'] if 'projections' in bandopt else False
        
        if 'bands' in plots or 'BZ' in plots:
            
            if parse: 
                bandruns.append(BSVasprun(loc + "/bands.xml", parse_projected_eigen=True))
                
            bands = bandruns[i].get_band_structure(loc + "/KPOINTS", line_mode=True,
                                                        efermi=bandruns[i].efermi)

        if 'bands' in plots:
            figsize = figsize if figsize else (7, 8)
            fig = plt.figure(figsize=figsize)
            gs = fig.add_gridspec(1,5, wspace=0.3) if 'dos' in plots \
                                                    else fig.add_gridspec(1,2)

#            axB = fig.add_subplot(gs[0:3], projection='My_Axes_Y')
            axB = fig.add_subplot(gs[0:3])
            axB.set_title('Band Structure', **title)
            axB.set_xlabel(r'Symmetry Points', **axes)
            axB.set_ylabel("E - E$\mathrm{_f}$ (eV)", **axes)
            axB.tick_params(axis='x', which='both', length=0, pad=5)
            axB.tick_params(axis='y', which='both',
                            labelsize=ticks['fontsize'])

            plotBands(axB, bandruns[i], **bandopt)

        if 'dos' in plots:
            
            if parse: dosruns.append(Vasprun(loc + "\..\dos.xml"))

            if 'spin' in dosopt:
                spin = dosopt['spin']
            else:
                if dosruns[i].is_spin:
                    spin = [Spin.up, Spin.down]
                else:
                    spin = [Spin.up]

            xlim = dosopt['xlim'] if 'xlim' in dosopt else (-4,4)
            ylim = dosopt['ylim'] if 'ylim' in dosopt else (-15,15)
            shift = dosopt['shift'] if 'shift' in dosopt else 0.
            if 'ylim' in dosopt:
                ylim = dosopt['ylim']
            else:
                e = dosruns[-1].tdos.energies - \
                    dosruns[-1].tdos.get_cbm_vbm()[1] - shift
                ind = [e.tolist().index(i) for i in 
                       e[(e > xlim[0]) & (e < xlim[1])]]
                lim = max([max(a[ind[0]:ind[-1]]) for a in 
                               dosruns[-1].tdos.densities.values()])
                ylim = (-lim,lim)
            
            vertical = dosopt['vertical'] if 'vertical' in dosopt else None

            if 'bands' in plots:
#                axD = fig.add_subplot(gs[3:5], sharey=axB, projection='My_Axes_Y')
                axD = fig.add_subplot(gs[3:5], sharey=axB)
                axD.tick_params(labelleft=False)
                setDos('y', 0.25, spin=spin, lim=ylim)
                plotDOS(fig, dosruns[i], **dosopt)
            else:
                if figsize:
                    figsize = figsize
                elif vertical:
                    figsize = (5,10)
                else:
                    figsize = (14,5)
                fig = plt.figure(figsize=figsize)

                if vertical:
#                    axD = fig.add_subplot(111,projection='My_Axes_Y')
                    axD = fig.add_subplot(111)
                    axD.set_ylim(xlim)
                    setDos('y', 0.25, spin=spin, lim=ylim)
                    plotDOS(fig, dosruns[i], orient='y', **dosopt)
                else:
#                    axD = fig.add_subplot(111, projection='My_Axes_X')
                    axD = fig.add_subplot(111)
                    axD.set_xlim(xlim)
                    setDos('x', 0.25, spin=spin, lim=ylim)
                    plotDOS(fig, dosruns[i], orient='x', **dosopt)

        if 'BZ' in plots: plotBZ()

        if save: 
            plt.savefig(loc + '/%s' % '_'.join(plots) + \
                        '.png', format='png', dpi=dpi)


def plotBands(ax, bandrun, **kargs):
    """ Plot band structure """

    def setTicks(ax, dk):
        """ Collect symmetry k-point labels and positions """

        if k.label:
            if k.label != prevK[0]:
                if k.label == branch['name'].split('-')[0]:
                    ax.axvline(prevK[1], color='k', ls='-', lw=1.5)
                    if k.label != k0: dk += 0.2
                    ax.axvline(d + dk, color='k', ls='-', lw=1.5)
                if k.label.startswith("\\") or k.label.find("_") != -1:
                    kLabels['label'].append('$' + k.label + '$')
                else:
                    kLabels['label'].append(k.label)
                kLabels['distance'].append(d + dk)
            else:
                ax.axvline(d + dk, color='k', ls='-', lw=0.4)

        return dk

    def plotBranch(ax, d, e, tol, i, c, ls):
        """ Plot each band over each k-point branch """

        global ax1
        ax1 = ax

        def contrib(species, data, el, orb):
            """ Calculate orbital contributions """

            # SPEED THIS UP BY UTILIZING NUMPY'S VECTORIZED CODE

            p = [i for i, s in enumerate(species) if s == el]

            contrib = []
            c = np.array(orbitals[el][orb][1])/255
            o = eval('OrbitalType.' + orb + '.value')
            for b in range(*bandRange):
                for k in range(len(d)-1):
                    contrib.append(
                        np.concatenate(
                            [c, [np.average([data[b,k+i,o,p[0]:p[-1]+1].sum()
                                             for i in range(2)])]]))
                    if smooth:
                        for i in range(int((tol - len(d)) / len(d))):
                            contrib.append(contrib[-1])
                if smooth:
                    for i in range(int((tol - len(d)) / len(d))):
                        contrib.append(contrib[-1])
            return contrib

        segs = [] # store segmented bands if plotting projections

        if smooth:

            if not tol: tol = len(d) - 1

            for b in range(*bandRange):

                tck = splrep(d, e[:,b])
                step = (d[-1] - d[0]) / (tol-1)
                xs = np.array([x * step + d[0] for x in range(tol)])
                ys = np.array([splev(x * step + d[0], tck, der=0)
                                            for x in range(tol)])
                for y in ys:
                    if np.isnan(y):
                        print(" -> Smoothing failure :( ")
                        break

                if projections:
                    pts = np.array([xs, ys]).T.reshape(-1, 1, 2)
                    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)
                    segs.append(seg)
                else:
                    ax.plot(xs, ys, c, lw=lw)

        else:
            if projections:
                for b in range(*bandRange):
                    pts = np.array([d, e[:,b]]).T.reshape(-1, 1, 2)
                    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)
                    segs.append(seg)
            else:
                ax.plot(d, e[:,bandRange[0]:bandRange[1]], c, lw=lw, picker=5)

        if projections:
            species = [str(i) for i in bands.structure.species]
            data = bands.projections[s][:,i[0]:i[1]+1]
            for el in elements:
                for orb in orbitals[el].keys():
                    colors = contrib(species, data, el, orb)
                    ax.add_collection(LineCollection(np.array(segs).reshape(-1,2,2),
                                                     colors=colors, lw=lw, linestyles=ls))

        if showVBM: ax.plot([d[0],d[-1]], [0,0], color='k', ls=':', lw=0.5)
        if showCBM: ax.plot([d[0],d[-1]], [Eg,Eg], color='k', ls=':', lw=0.5)

    ###########################################################################

    def onpick(event):
        print('Band Number -> ', event.artist.get_label())
    
    fig = ax.get_figure()
    fig.canvas.mpl_connect('pick_event', onpick)

    ###########################################################################

    bandRange = kargs['bandRange'] if 'bandRange' in kargs else [0,bands.nb_bands-1]
    xlim = kargs['xlim'] if 'xlim' in kargs else (0,len(bands.branches))
    ylim = kargs['ylim'] if 'ylim' in kargs else [-4,4]
    spin = kargs['spin'] if 'spin' in kargs else [Spin.up, Spin.down]
    zero_to_fermi = kargs['zero_to_fermi'] if 'zero_to_fermi' in kargs else None
    showVBM = kargs['showVBM'] if 'showVBM' in kargs else None
    showCBM = kargs['showCBM'] if 'showCBM' in kargs else None
    bandMarks = kargs['bandMarks'] if 'bandMarks' in kargs else None
    projections = kargs['projections'] if 'projections' in kargs else None
    lw = kargs['lw'] if 'lw' in kargs else 1
    if 'smooth' in kargs:
        smooth = kargs['smooth'][0]
        tol = kargs['smooth'][1]
    else:
        smooth = None
        tol = None

    ###########################################################################
    
    Eg = bands.get_band_gap()['energy']
    try:
        VBM = bands.get_vbm()['energy']
        CBM = bands.get_cbm()['energy']
        print('\nBand Gap = %0.3f eV\nVBM = %0.3f eV\nCBM = %0.3f eV' % (Eg, VBM, CBM))
    except:
        pass

    ###########################################################################

    if no_shift:
        zero_energy = 0.
    elif bands.is_metal():
        zero_energy = bands.efermi
    elif zero_to_fermi:
        zero_energy = bands.efermi
    else:
        zero_energy = bands.get_vbm()['energy']

    elements = bands.structure.symbol_set

    ###########################################################################

    k0 = bands.kpoints[0].label
    prevK = ('',0); dk = 0
    kLabels = {'label': [], 'distance': []}
    totD = list()
    for branch in bands.branches[xlim[0]:(xlim[1]+1)]:

        i = (branch['start_index'], branch['end_index'])

        # calculate x-axis data
        distances = []
        for k, d in zip(bands.kpoints[i[0]:i[1]+1], bands.distance[i[0]:i[1]+1]):
            dk = setTicks(ax, dk)
            distances.append(d + dk)

        # calculate y-axis data
        e = bandrun.eigenvalues
        energies = {'1': None, '-1': None}
        for s in spin:
            energies[str(s)] = e[s][i[0]:i[1]+1,:,0]
            plotBranch(ax, np.array(distances), energies[str(s)] - zero_energy,
                       tol, i, 'k' if s.value == 1 else 'r', '-' if s.value == 1 else '--')

        prevK = (k.label, d + dk)
        totD.append(distances)

    if bandMarks:
        shift = bands.branches[xlim[0]]['start_index']
        totD = np.array(totD).flatten()
        cMarks = 'k' if projections else ['r']
        for i in [bands.get_vbm(), bands.get_cbm()]:
            ax.scatter(totD[i['kpoint_index'][0] - shift], i['energy'] - zero_energy, color = cMarks)

    ax.set_xticks(kLabels['distance'])
    ax.set_xticklabels(kLabels['label'], **ticks)
    ax.set_xlim(kLabels['distance'][0], kLabels['distance'][-1])

    if bandRange[0] != 0:
        yMin = 0; yMax = 0
        for s in spin:
            eMin = (energies[str(s)][:,bandRange[0]] - zero_energy).min()
            eMax = (energies[str(s)][:,bandRange[1]] - zero_energy).max()
            if yMin > eMin: yMin = eMin
            if yMax < eMax: yMax = eMax
        ax.set_ylim(np.floor(yMin),np.ceil(yMax))
    else:
        ax.set_ylim(ylim)
    

def plotDOS(fig, dosrun, orient='y', **kargs):

    from scipy.interpolate import interp1d

    def formatDos(x):
        
        l = np.array([0 if i == 0 else 1 for i in interp1d(e, x)(y)])
        start = []; stop = []; prevL = None
        for i, j in enumerate(l):
            if j == 1 and prevL == 0: start.append(i)
            if len(start) > 0 and (j == 0 or i == len(l)-1) and prevL == 1: stop.append(i)
            prevL = j
        
        if smooth:
            func = splev(y, splrep(e,x)) * l
        else:
            func = x * l

        return func, [[i, j + 1] for i, j in zip(start, stop)]

    def plotIntervals(ax, x, y, o, intervals, t, c, a, lw):
        
        if orient != 'y': x, y = y, x
        
#        global tempX, tempY
#        tempX = x
#        tempY = y

        if shade: # and ax.get_label() != 'zoom':
            ax.fill_between(x, y, edgecolor=c/255, facecolor=c/255, alpha=a, lw=0.)
            handles.append(ax.collections[-1])
        else:
            if smooth:
                for i, j in intervals:
                    ax.plot(x[i:j][eval(o)], y[i:j][eval(o)], color=c/255, lw=lw)
            else:
                ax.plot(x, y, color=c/255, lw=lw, label=t)
            
            handles.append(ax.lines[-1])
#            ax.fill_between(x[i:j][eval(o)], y[i:j][eval(o)],
#                            edgecolor=c*0, facecolor=c/255, alpha=a, lw=0)

        ax.tick_params(axis='both', which='both',
                            labelsize=ticks['fontsize'], 
#                            pad=ticks['pad']
                            )

        labels.append(t)

    ###########################################################################

    if 'spin' in kargs:
        spin = kargs['spin']
    else:
        if dosrun.is_spin:
            spin = [Spin.up, Spin.down]
        else:
            spin = [Spin.up]
    shade = None
    smooth = None
    if 'shade' in kargs: shade = kargs['shade']
    if 'smooth' in kargs: smooth = kargs['smooth']
    shift = kargs['shift'] if 'shift' in kargs else 0.
    plotOrbs = kargs['plotOrbs'] if 'plotOrbs' in kargs else False

    elements = list(dict.fromkeys(dosrun.atomic_symbols))
    labels = list(); handles = list()
    lw = kargs['lw'] if 'lw' in kargs else 1

    e = dosrun.tdos.energies - dosrun.tdos.get_cbm_vbm()[1] - shift; y = e
    if no_shift: e = dosrun.tdos.energies; y = e
    if smooth: y = np.arange(e[0],e[-1],0.0005)

    for sp in spin:

        i = int(sp)
        if orient == 'y':
            o = 'x[i:j]>=0' if i == 1 else 'x[i:j]<=0'
        else:
            o = 'y[i:j]>=0' if i == 1 else 'y[i:j]<=0'

        # total DOS
        c = np.array([0,0,0]); a = 0.1
        x, inter = formatDos(dosrun.tdos.densities[sp] * i)
#        for ax in _axes: plotIntervals(ax, x, y, o, inter, 'Total', c, a, lw)

        # partial DOS
        if plotOrbs:
            for name in elements:
    
                # spd projected DOS
                spd_dos = dosrun.complete_dos.get_element_spd_dos(name)
    
                # plot selected orbitals
                for orb, props in orbitals[name].items():
                    c = np.array(props[1]); a = 0.25
                    x, inter = formatDos(
                        spd_dos[eval('OrbitalType.' + orb)].densities[sp] * i)
                    for ax in _axes: plotIntervals(
                        ax, x, y, o, inter, f"{name} {props[0]}", c, a, lw)

    from collections import OrderedDict
    legend = OrderedDict(zip(labels, handles))

    fig.legend(legend.values(), legend.keys(), bbox_to_anchor=(0.98,0.93), fontsize=16)


def plotBZ():

    from pymatgen.electronic_structure.plotter import \
    plot_brillouin_zone_from_kpath as pltBZ

    st = bands.structure
    st_sym = sga(st)
    prim = sga(st).get_primitive_standard_structure(
            international_monoclinic=False)

    t = '{0} - {1} ({2})'.format(st_sym.get_lattice_type().capitalize(),
                                 st_sym.get_space_group_symbol(),
                                 st_sym.get_space_group_number())

    pltBZ(HighSymmKpath(prim),title=t)


#%%############################### MAIN #######################################

title = {'fontsize': 20, 'pad': 20}
axes = {'fontsize': 20, 'labelpad': 8}
ticks = {'fontsize': 16} # , 'pad': 10
arrows = {'fontsize': 12, 'labelpad': 20}
font={'family': 'times new roman'}
opt = {'font': font}
for i in opt: plt.rc(i,**opt[i])

plt.close()

orbitals = {
    'Al': {
            's': ['3s', [129,178,214]],
            'p': ['3p', [64,128,128]],
            },
    'Ti': {
            'd': ['3d', [182,175,169]]
            },
    'Cu': {
            's': ['4s', [231,102,80]], 
            'd': ['3d', [180,85,15]]
            },
    'Zn': {
            'd': ['3d', [186,196,200]]},
    'Bi': {
            's': ['6s', [248,153,231]], 
            'p': ['6p', [128,0,128]]
            },
    'In': {
            's': ['5s', [5,67,235]],
            'd': ['4d', [0,162,232]]
            },
    'W': {
#            's': ['6s', [91,155,213]], 
#            'p': ['6p', [112,173,71]], 
            'd': ['5d', [120,120,120]]
            },
    'V': {
            'd': ['3d', [215,165,0]]
            },
    'Sb': {
            's': ['5s', [204,204,102]], 
            'p': ['5p' , [102,153,77]]
            },
    'Ta': {
            'd': ['5d', [183,154,86]]
            },
    'O': {
#            's': ['2s', [100,0,0]],
            'p': ['2p', [154,0,0]]
            },
    'S': {
#            's': ['3s', [200,0,0]],
            'p': ['3p', [251,205,111]]
            },
    'N': {
            'p': ['3p', 'white']
            }
    }

dosopt = {
    'xlim': (-4,4),
#     'ylim': (-15,15),
#     'DOSticks': True,
#     'shade': True,
     'spin': [Spin.up], # comment out for both spins
#     'zoom': [(-1,0),(0,5)], #,[0.55, 0.35, 0.2, 0.6]], 
     'lw': 0.5,
    'smooth': True,
#     'shift': 0.04, # positive shifts down
#     'vertical': True,
     'plotOrbs': True
    }

bandopt = {
    'ylim': (-8,4),
#     'xlim': (0,3), # numbers represent branches (0 is first)
    'spin': [Spin.up], # comment out for both spins
    'showVBM': True,
    'showCBM': True,
#     'zero_to_fermi': True,
#     'bandRange': [70,80],
#     'bandMarks': True,
#    'projections': True,
    'lw': .5,
#    'smooth': [True, 100]
    }

no_shift = False

plots = [
        'bands',
        'dos',
#        'BZ'
        ]

if 'zoom' in dosopt:
    if 'bands' in plots or 'vertical' in dosopt: 
        dosopt['zoom'].reverse()

locs = [
#        r'E:\Research\VASP Data\Bi2WO6\bands\full',
#        r'E:\Research\VASP Data\Sb-Bi2WO6\SYM\1\bands\full',
#        r'E:\Research\VASP Data\Sb-Bi2WO6\SYM\2\bands\full',
#        r'E:\Research\VASP Data\Sb-Bi2WO6\SYM\3\bands\full',
#        r'E:\Research\VASP Data\Sb-Bi2WO6\SYM\4\bands\full',
#        r'E:\Research\VASP Data\Sb-Bi2WO6\SYM\5\bands\full',
#        r'E:\Research\VASP Data\Sb-Bi2WO6\SYM\6\bands\full',
#        r'E:\Research\VASP Data\Sb-Bi2WO6\SYM\full\bands\full'
#        r'E:\Research\VASP Data\Bi2WO6\prim\bands\128\H-Y-Gm-Z-D-Y',
#        r'E:\Research\VASP Data\S-Bi2WO6\plotData\D_1-Y_1-H_2-X-A'
#        r'E:\Research\VASP Data\WO3\Pristine\full',
#        r'E:\Research\VASP Data\WO3\157-tilt\full',
#        r'E:\Research\VASP Data\WO3\145-tilt\full',
#        r'E:\Research\VASP Data\BiVO4\mono\plots\Gm-A-L-M-Y-V-Gm',
        r'E:\Research\VASP Data\BiVO4\ortho\plots\Gm-A-L-M-Y-V-Gm',
#        r'E:\Research\VASP Data\BIVO\1\plots\Gm-A-L-M-Y-V-Gm',
#        r'E:\Research\VASP Data\BIVO\2\plots\Gm-A-L-M-Y-V-Gm',
#        r'E:\Research\VASP Data\BIVO\3\plots\Gm-A-L-M-Y-V-Gm'
#        r'E:\Research\VASP Data\SbTaO4\mono\plots'
        r'E:\Research\VASP Data\BAVO\1\plots\full',
        r'E:\Research\VASP Data\BAVO\2\plots\full',
        r'E:\Research\VASP Data\BAVO\3\plots\full'
       ]

plot(locs,
     parse=True,
#     figsize=(8,4),
     dpi=1200,
     lw=0.4,
#     save=True
     )