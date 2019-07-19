"""
read in all the histograms and save dmdts
optionally make plots. plot raw with alpha to compare. go through these plots to
decide what window lens, tmin, tmax to choose
comment out some ds if just ran a few new simulations and don't want to do everything

run like:
python load_histograms_save_dmdts.py --cluster=fend
"""
## TODO need to shift dmdes


import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.io import ascii
from scipy.ndimage.filters import gaussian_filter
from scipy.integrate import simps
from scipy.integrate import cumtrapz
from scipy import interpolate
import os
import argparse, sys

parser = argparse.ArgumentParser()
parser.add_argument('--cluster', help='cluster, e.g., hyades/fend/pfe', type=str, default='fend')
args = parser.parse_args()

if args.cluster == 'hyades': execfile('/pfs/lawsmith/jYT/my_settings.py')
elif args.cluster == 'fend': execfile('/groups/dark/lawsmith/jYT/my_settings.py')
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['font.size'] = 18

M_bh = 1e6*M_sun

# if choosing what smoothing to use
LOOP_THRU_SIGMAS = True
sigmas = [1,5,10,15,20,25,30,35,40,45,50]

PLOT_DMDES = True

"""
ds = [
    # run,                  chk,        sigma,  tmin,	tmax,     slope, split, sigma2
    ['m1.0_p1_b1.0',        '0100',     20,     -2,	    2,        0,      0,  0],
    ['m1.0_p1_b1.0_pfe',    '0100',     20,     -2,	2, 0, 0, 0], #too messy
    ['m1.0_p1_b1.5',		'0090',		25,		-2,	    2,        0,      0,  0],
    ['m1.0_p1_b1.75',		'0090',		50,		-2,     2,        0,      0,  0],
    ['m1.0_p1_b2.0',		'0080',		30,		-2,	    2,     0, 0,  0],
    ['m1.0_p1_b3.0',		'0060',		30,		-2,	    2,        0,  0,  0],
    ['m1.0_p10_b1.0',		'0100',     20,     -2,	    2, 0, 0, 0], # too messy
    ['m1.0_p10_b1.0_48k',	'0100',     35,     -2,	  2, 0, 0, 0], # too messy
    ['m1.0_p10_b1.5',		'0090',     20,     -2,	    2,        0,  0,  0],
    ['m1.0_p10_b2.0',		'0075',     25,     -2,	    2,        0, 0,  0],
    ['m1.0_p10_b2.0_48k',	'0080',     50,     -2,	    2,        0,  0,  0],
    ['m1.0_p10_b2.5',		'0060',     50,     -2,	    2,        0,  0,  0],
    ['m1.0_p10_b3.0',		'0060',     30,     -2,	    2,        0,  0,  0],
    ['m1.0_p10_b3.0_48k',   '0057',     20,     -2,	    2,        0,  0,  0], # very different shape
    ['m1.0_p10_b4.0',		'0050',     40,     -2,	    2,        0,  0,  0],
    ['m1.0_p16_b1.0',		'0100',     50,     -2,	2, 0, 0, 0], # too messy
    ['m1.0_p16_b1.0_pfe',	'0052',     50,     -2,	    2,        0,  0,  0], # bit different shape
    ['m1.0_p16_b1.5',		'0090',     50,     -2,	    2,        0,  0,  0],
    ['m1.0_p16_b2.0',		'0075',     35,     -2,   2,        0,  0,  0],
    ['m1.0_p16_b3.0',		'0060',     50,     -2,   2,        0,  0,  0],
    ['m1.0_p16_b4.0',		'0050',     50,     -2,	    2,        0,  0,  0],
]

"""
ds = [
    #['m1.0_p1_b1.0',		'0100',     140,     -2,	2,        0,  0,  0],
#    ['m1.0_p1_b1.0',        '0060',     80,     -2,    2,        0,  0,  0],
    #['m1.0_p1_b1.5',		'0090',     20,     -2,	2,        0,  0,  0],
#    ['m1.0_p1_b2.0',		'0080',     20,     -2,	2,        0,  0,  0],
#    ['m1.0_p1_b3.0',		'0060',     20,     -2,	2,        0,  0,  0],

    #['m1.0_p10_b1.0',		'0100',     140,     -2,	2,        0,  0,  0],
#    ['m1.0_p10_b1.0',       '0060',     160,     -2,    2,        0,  0,  0],
#    ['m1.0_p10_b2.0',		'0075',     80,     -2,	2,        0,  0,  0],
#    ['m1.0_p10_b3.0',		'0060',     30,     -2,	2,        0,  0,  0],
    #['m1.0_p10_b3.0_24k',		'0060',     20,     -2,	2,        0,  0,  0],

    #['m1.0_p16_b1.0',		'0100',     140,     -2,	2,        0,  0,  0],
#    ['m1.0_p16_b1.0',       '0060',     160,     -2,    2,        0,  0,  0],
#    ['m1.0_p16_b2.0',		'0075',     80,     -2,	2,        0,  0,  0],
#    ['m1.0_p16_b3.0',		'0060',     60,     -2,	2,        0,  0,  0],
    ['m1.0_p16_b3.0_96k_1000', '0060',  60,     -2, 2,        0,  0,  0],
    #['m1.0_p16_b3.0_24k',		'0060',     20,     -2,	2,        0,  0,  0],
#     ['m1.0_p16_b4.0',		'0050',     20,     -2,	2,        0,  0,  0],

#    ['m3.0_p1_b2.0_48k',       '0070',     20,     -2, 2,        0,  0,  0],
#    ['m3.0_p16_b4.0_48k',       '0050',     20,     -2, 2,        0,  0,  0],

    #['43_b1.0_48k',		'0100',     40,     -2,	2,        0,  0,  0],
    #['43_b1.5_48k',		'0090',     20,     -2,	2,        0,  0,  0],
    #['43_b1.5_48k',		'0157',     20,     -2,	2,        0,  0,  0],
    #['43_b1.5_48k',		'0200',     20,     -2,	2,        0,  0,  0],
    #['43_b1.5_300k',     '0090',     20,     -2, 2,        0,  0,  0],
    #['43_b1.5_300k',     '0200',     20,     -2, 2,        0,  0,  0],
    #['43_b1.5_300k',     '0200',     20,     -2, 2,        0,  0,  0],
    #['43_b2.0_48k',		'0080',     20,     -2,	2,        0,  0,  0],
    #['43_b2.0_24k',		'0080',     20,     -2,	2,        0,  0,  0],
]
"""
ds = [
    #['m1.0_p10_b3.0',		'0043',     20,     -2,	2,        0,  0,  0],
    #['m1.0_p10_b3.0',		'0030',     1,     -2,	2,        0,  0,  0],
    #['43_b2.0_6k',		'0040',     5,     -2,	2,        0,  0,  0],
    #['43_b2.0_6k',		'0060',     10,     -2,	2,        0,  0,  0],
    #['43_b2.0_6k',		'0080',     20,     -2,	2,        0,  0,  0],
    #['43_b2.0_12k',		'0041',     1,     -2,	2,        0,  0,  0],
    #['43_b2.0_12k',		'0080',     20,     -2,	2,        0,  0,  0],
    ['43_b2.0_48k',		'0080',     20,     -2,	2,        0,  0,  0],
    #['m1.0_p10_b3.0',		'0060',     20,     -2,	2,        0,  0,  0],
    #['m1.0_p10_b3.0_24k',		'0060',     20,     -2,	2,        0,  0,  0],
    #['m1.0_p16_b3.0',		'0060',     20,     -2,	2,        0,  0,  0],
    #['m1.0_p16_b3.0_24k',		'0060',     20,     -2,	2,        0,  0,  0],
    #['m1.0_p16_b3.0_12k',		'0060',     20,     -2,	2,        0,  0,  0],
    #['m1.0_p16_b3.0_12k',		'0040',     5,     -2,	2,        0,  0,  0],
    #['43_b2.0_12k_flash',		'0040',     5,     -2,	2,        0,  0,  0],
    #['43_b2.0_12k_flash',		'0080',     20,     -2,	2,        0,  0,  0]
    ]

ds = [
    ['43_b2.0_24k_intel',		'0080',     20,     -2,	2,        0,  0,  0],
    ['43_b2.0_48k_intel',		'0080',     20,     -2,	2,        0,  0,  0],
    ]

ds = [
    # run,                  chk,        sigma,  tmin,	tmax,     slope, split, sigma2
    # p1
    ['m1.0_p1_b1.0',        '0040',     15,     -1.4,	-0.1,     0.2,    0,  0],
    ['m1.0_p1_b1.0',        '0100',     20,     -2,	    2,        0,      0,  0],
    #['m1.0_p1_b1.0_pfe',    '0100',     15,     -1.4,	-0.1], too messy
	['m1.0_p1_b1.5',		'0040',		10,		-1.65,	0.1,      0,      0,  0],
    ['m1.0_p1_b1.5',		'0090',		25,		-2,	    2,        0,      0,  0],
	['m1.0_p1_b1.75',		'0040',		15,		-1.75,	-0.35,    0,      0,  0],
    ['m1.0_p1_b1.75',		'0090',		50,		-2,     2,        0,      0,  0],
	['m1.0_p1_b2.0',		'0040',		25,		-1.8,	-0.65,    0,   0,  0],
    ['m1.0_p1_b2.0',		'0080',		30,		-2,	    -0.2,     0, 0,  0],
	['m1.0_p1_b3.0',		'0040',		15,		-1.9,	-0.3,     0,  0,  0],
    ['m1.0_p1_b3.0',		'0060',		30,		-2,	    2,        0,  0,  0],
    # p10
	['m1.0_p10_b1.0',		'0040',     20,     -1.3,	-0.3,     0.3,    -0.5, 5],
    #['m1.0_p10_b1.0',		'0100',     20,     -1.3,	-0.3], # too messy
	['m1.0_p10_b1.0_256',	'0040',     35,     -1.3,	-0.15,    0.5,    -0.5, 15],
    #['m1.0_p10_b1.0_256',	'0100',     35,     -1.3,	-0.15], # too messy
	['m1.0_p10_b1.5',		'0040',     20,     -1.65,	-0.3,     0.5, 0,  0],
    ['m1.0_p10_b1.5',		'0090',     20,     -2,	    2,        0,  0,  0],
	['m1.0_p10_b2.0',		'0040',     25,     -1.8,	-0.1,     0.3,    0,  0],
    ['m1.0_p10_b2.0',		'0075',     25,     -2,	    2,        0, 0,  0],
	['m1.0_p10_b2.0_256',	'0040',     50,     -1.8,	2,        0,  0,  0],
    ['m1.0_p10_b2.0_256',	'0080',     50,     -2,	    2,        0,  0,  0],
	['m1.0_p10_b2.5',		'0040',     20,     -1.9,	-0.3,     0,  0,  0],
    ['m1.0_p10_b2.5',		'0060',     50,     -2,	    2,        0,  0,  0],
	['m1.0_p10_b3.0',		'0040',     20,     -2,	    -0.4,     0,  0,  0],
    ['m1.0_p10_b3.0',		'0060',     30,     -2,	    2,        0,  0,  0],
    ['m1.0_p10_b3.0_256',   '0040',     20,     -2,	    2,        0,  0,  0],
    ['m1.0_p10_b3.0_256',   '0057',     20,     -2,	    2,        0,  0,  0], # very different shape
	['m1.0_p10_b4.0',		'0040',     40,     -2,	    2,        0,  0,  0],
    ['m1.0_p10_b4.0',		'0050',     40,     -2,	    2,        0,  0,  0],
    # p16
	['m1.0_p16_b1.0',		'0040',     50,     -1.1,	2,        0,  0,  0],
    #['m1.0_p16_b1.0',		'0100',     50,     -1.1,	2], # too messy
    ['m1.0_p16_b1.0_pfe',	'0052',     50,     -1.1,	2,        0,  0,  0], # bit different shape
	['m1.0_p16_b1.5',		'0040',     50,     -1.5,	2,        0,  0,  0],
    ['m1.0_p16_b1.5',		'0090',     50,     -2,	    2,        0,  0,  0],
	['m1.0_p16_b2.0',		'0040',     30,     -1.65,	2,        0,  0,  0],
    ['m1.0_p16_b2.0',		'0075',     35,     -1.6,   2,        0,  0,  0],
	['m1.0_p16_b3.0',		'0040',     50,     -1.85,	2,        0,  0,  0],
    ['m1.0_p16_b3.0',		'0060',     50,     -1.8,   2,        0,  0,  0],
	['m1.0_p16_b4.0',		'0040',     30,     -2,	    2,        0,  0,  0],
    ['m1.0_p16_b4.0',		'0050',     50,     -2,	    2,        0,  0,  0],
]
"""


#els = ['ev', 'h1', 'he4', 'o16', 'c12', 'ne20', 'n14']
els = ['ev', 'h1', 'he3', 'he4', 'c12', 'c13', 'n14', 'o16', 'ne20', 'na23', 'mg24', 'al27', 'si28', 's34']
#els = ['ev']
"""
if LOOP_THRU_SIGMAS == True:
    for d in ds:
        for el in els:
            for sigma in sigmas:
                fig2, ax2 = plt.subplots()
                e, dm = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/'+d[0]+
                    '/b10000_'+el+'_bhbound_histogram_multitidal_hdf5_chk_'+d[1]+'.dat',
                    skiprows=4)
                de = e[1]-e[0]
                dm_de = dm/de
                log_dm_de = np.log10(dm_de)
                # smooth linear value of dmde first
                s_dm_de = gaussian_filter(dm_de, sigma, mode='wrap')
                slog_dm_de = np.log10(s_dm_de)
                e_bound = e[np.where(e<0.)]
                t = 2.*np.pi*G*M_bh/((2*np.abs(e_bound))**(3./2))
                de_dt = (1./3)*((2*np.pi*G*M_bh)**(2./3))*t**(-5./3)
                dm_de_bound = 10**slog_dm_de[np.where(e<0)]
                dm_de_bound_orig = dm_de[np.where(e<0)]
                mdot = dm_de_bound*de_dt
                mdot_orig = dm_de_bound_orig*de_dt

                log_t_yr = np.log10(t/yr)
                log_mdot_moyr = np.log10(mdot*yr/M_sun)
                log_mdot_moyr_orig = np.log10(mdot_orig*yr/M_sun)

                # plot unsmoothed
                ax2.scatter(log_t_yr, log_mdot_moyr_orig, c='C0', rasterized=True,
                	alpha=0.5, s=1, edgecolors='none')

                sel = np.where((log_t_yr>d[3]) & (log_t_yr<d[4]))[0]
                x = log_t_yr[sel]
                y = log_mdot_moyr[sel]
                ax2.plot(x, y, c='C1')

                # extend dmdt from slope near. could maybe improve slope function
                if x[-1] < 2.0:
                    ninf = (y[-20] - y[-1]) / (x[-20] - x[-1])
                    ax2.plot([x[-1], 2.0], [y[-1], y[-1] + (2.0 - x[-1])*ninf], c='C2')

                if PLOT_DMDES:
                    fig, ax = plt.subplots()
                    ax.scatter(e/1e17, log_dm_de, c='C0', rasterized=True,
        				alpha=0.5, s=1, edgecolors='none')
                    ax.plot(e/1e17, slog_dm_de, c='C1')
                    ax.set_xlim(-10, 10)
                    ax.set_ylim(9, 15)
                    ax.set_xlabel(r'$E\ \mathrm{[10^{17}\ erg\ g^{-1}]}$')
                    ax.set_ylabel(r'$\log\ dM/dE\ \mathrm{[g^2\ erg^{-1}]}$')
                    fig.tight_layout(pad=0.3)
                    directory = '/pfs/lawsmith/FLASH4.3/runs/results/dmdes/gaussian_wrap/' + d[0].replace(".","_")
                    if not os.path.exists(directory): os.makedirs(directory)
                    fig.savefig(directory + '/dmde_'
                        + d[0].replace(".","_") + '_' + d[1] + '_' + el + '_'
                        + str(sigma) + '.pdf')

                ax2.set_xlim(-2, 2)
                ax2.set_ylabel(r'$\log\ \dot M\ {\rm [M_\odot/yr]}$')
                ax2.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
                fig2.tight_layout(pad=0.3)
                directory = '/pfs/lawsmith/FLASH4.3/runs/results/dmdts/dmde_first/' + d[0].replace(".","_")
                if not os.path.exists(directory): os.makedirs(directory)
                fig2.savefig(directory + '/dmdt_'
                	+ d[0].replace(".","_") + '_' + d[1] + '_' + el + '_'
                	+ str(sigma) + '.pdf')

                plt.close('all')
"""

# commenting above out because of some weird behavior where it wouldn't loop
# over els; would only do 6

# plot final versions
#else:
for d in ds:
    for i, el in enumerate(els):
        fig2, ax2 = plt.subplots()
        fig3, ax3 = plt.subplots()
        if args.cluster == 'hyades': directory = '/pfs/lawsmith/FLASH4.3/runs/'
        elif args.cluster == 'fend': directory = '/storage/dark/lawsmith/FLASH4.3/runs/'
        e, dm = np.loadtxt(directory+d[0]+'/b10000_'+el+'_bhbound_histogram_multitidal_hdf5_chk_'+d[1]+'.dat', skiprows=4)
        de = e[1]-e[0]
        dm_de = dm/de
        log_dm_de = np.log10(dm_de)
        # smooth linear value of dmde first
        s_dm_de = gaussian_filter(dm_de, d[2], mode='wrap')
        slog_dm_de = np.log10(s_dm_de)

        # shift dmde histogram to be centered at 0
        if d[0] == 'm1.0_p1_b1.0':
            e = e - e[(-1<e/1e17) & (e/1e17<1)][np.argmax(slog_dm_de[(-1<e/1e17) & (e/1e17<1)])]
        if d[0] == 'm1.0_p1_b2.0':
            e = e + 0.1e17
        if d[0] == 'm1.0_p1_b3.0':
            e = e - e[(-1<e/1e17) & (e/1e17<1)][np.argmin(slog_dm_de[(-1<e/1e17) & (e/1e17<1)])]
        if d[0] == 'm1.0_p10_b1.0':
            e = e - e[(-1<e/1e17) & (e/1e17<1)][np.argmin(slog_dm_de[(-1<e/1e17) & (e/1e17<1)])]
        if d[0] == 'm1.0_p10_b2.0':
            e = e - e[np.argmax(slog_dm_de)]
        if d[0] == 'm1.0_p10_b3.0':
            pass
        if d[0] == 'm1.0_p16_b1.0':
            e = e - e[(-1<e/1e17) & (e/1e17<1)][np.argmin(slog_dm_de[(-1<e/1e17) & (e/1e17<1)])]
        if d[0] == 'm1.0_p16_b2.0':
            e = e - e[(-1<e/1e17) & (e/1e17<1)][np.argmax(slog_dm_de[(-1<e/1e17) & (e/1e17<1)])]
        if d[0] == 'm1.0_p16_b3.0':
            e = e - e[(-1<e/1e17) & (e/1e17<1)][np.argmax(slog_dm_de[(-1<e/1e17) & (e/1e17<1)])]
        if d[0] == 'm1.0_p16_b3.0_96k_1000':
            e = e - e[(0<e/1e17) & (e/1e17<1)][np.argmax(slog_dm_de[(0<e/1e17) & (e/1e17<1)])]

        if d[0] == 'm1.0_p16_b4.0':
            e = e - e[(0<e/1e17) & (e/1e17<0.3)][np.argmin(slog_dm_de[(0<e/1e17) & (e/1e17<0.3)])]
        if d[0] == 'm3.0_p1_b2.0_48k':
            e = e - e[(-1<e/1e17) & (e/1e17<1)][np.argmin(slog_dm_de[(-1<e/1e17) & (e/1e17<1)])]
        if d[0] == 'm3.0_p16_b4.0_48k':
            e = e - e[np.argmax(slog_dm_de)]

        e_bound = e[np.where(e<0.)]
        t = 2.*np.pi*G*M_bh/((2*np.abs(e_bound))**(3./2))
        de_dt = (1./3)*((2*np.pi*G*M_bh)**(2./3))*t**(-5./3)
        dm_de_bound = 10**slog_dm_de[np.where(e<0)]
        dm_de_bound_orig = dm_de[np.where(e<0)]
        mdot = dm_de_bound*de_dt
        mdot_orig = dm_de_bound_orig*de_dt

        log_t_yr = np.log10(t/yr)
        log_mdot_moyr = np.log10(mdot*yr/M_sun)
        log_mdot_moyr_orig = np.log10(mdot_orig*yr/M_sun)

        # plot unsmoothed
        ax2.scatter(log_t_yr, log_mdot_moyr_orig, c='C0', rasterized=True,
           alpha=0.5, s=1, edgecolors='none')

        # hacky thing to fix extended slope for a few dmdts
        # TODO change this
        if d[6] != 0:
            slog_dm_de1 = gaussian_filter(log_dm_de, d[2], mode='wrap')
            dm_de_bound1 = 10**slog_dm_de1[np.where(e<0)]
            mdot1 = dm_de_bound1*de_dt
            log_mdot_moyr1 = np.log10(mdot1*yr/M_sun)

            slog_dm_de2 = gaussian_filter(log_dm_de, d[7], mode='wrap')
            dm_de_bound2 = 10**slog_dm_de2[np.where(e<0)]
            mdot2 = dm_de_bound2*de_dt
            log_mdot_moyr2 = np.log10(mdot2*yr/M_sun)

            # join them
            sel1 = np.where(log_t_yr < d[6])[0]
            sel2 = np.where(log_t_yr > d[6])[0]
            log_mdot_moyr = np.append(log_mdot_moyr1[sel1], log_mdot_moyr2[sel2])


        sel = np.where((log_t_yr>d[3]) & (log_t_yr<d[4]))[0]
        x = log_t_yr[sel]
        y = log_mdot_moyr[sel]
        ax2.plot(x, y, c='C1')
        ax3.plot(x, (y-np.min(y[np.isfinite(y)]))/(np.max(y)-np.min(y)), c='C1')

		# extend dmdt from slope near end. could maybe improve slope function
        # this should just apply to the else statement above. could combine?
        if d[4] < 2.0:
            ninf = ((y[-10] - y[-1]) / (x[-10] - x[-1])) + d[5]
            yext = np.linspace(y[-1], y[-1] + (2.0 - x[-1])*ninf, num=100)[1:]
            xext = np.linspace(x[-1], 2.0, num=100)[1:]
            ax2.plot(xext, yext, c='C2')
            ax3.plot(xext, (yext-np.min(yext))/(np.max(yext)-np.min(yext)), c='C2')
            x = np.append(x, xext)
            y = np.append(y, yext)
            ax3.plot(x, (y-np.min(y[np.isfinite(y)]))/(np.max(y)-np.min(y)), c='C3')

        # interpolate for smoother integrating/interpolating
        # doesn't change integration answer though.
        f = interpolate.interp1d(x, y)
        xnew = np.linspace(min(x), max(x), num=10000)
        ynew = f(xnew)
        x = xnew
        y = ynew

        # write for later plotting
        if args.cluster == 'hyades': directory = '/pfs/lawsmith/results/dmdts/data/'
        elif args.cluster == 'fend': directory = '/groups/dark/lawsmith/results/'+d[0].replace(".","_")+'/'
        ascii.write([x, y], directory+d[0].replace(".","_")+'_'+d[1]+'_'+el+'.dat', overwrite=True, names=['log_t_yr','log_mdot_moyr'])

        """
        # integrate mdot curves, to compare with delta m, and save for later
        if el == 'ev':
            int_trapz = np.trapz(10**y, 10**x)
            int_simps = simps(10**y, 10**x)
            if args.cluster == 'hyades':
                ascii.write([[int_trapz], [int_simps]],
                    '/pfs/lawsmith/results/dmdts/integrals/'+d[0].replace(".","_")+
                    '_'+d[1]+'_'+el+'_'+str(d[2])+'.dat', overwrite=True, names=['int_trapz','int_simps'])
            elif args.cluster == 'fend':
                ascii.write([[int_trapz], [int_simps]],
                    '/groups/dark/lawsmith/results/'+d[0].replace(".","_")+'/'+
                    'ints_'+d[0].replace(".","_")+
                    '_'+d[1]+'_'+el+'_'+str(d[2])+'.dat', overwrite=True, names=['int_trapz','int_simps'])

            # cumulative plots
            int_cumtrapz = cumtrapz(10**y, 10**x, initial=0)
            ax3.plot(x, int_cumtrapz)
            ax3.set_xlim(-2, 2)
            ax3.set_ylim(0, 1)
            #ax3.set_ylabel(r'$\log\ \dot M\ {\rm [M_\odot/yr]}$')
            ax3.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
            fig3.tight_layout(pad=0.3)
            if args.cluster == 'hyades':
                directory = '/pfs/lawsmith/results/dmdts/integrals/'
            elif args.cluster == 'fend':
                directory = '/groups/dark/lawsmith/results/'+d[0].replace(".","_")+'/'+'ints_'
            if not os.path.exists(directory): os.makedirs(directory)
            fig3.savefig(directory \
            	+ d[0].replace(".","_") + '_' + d[1] + '_' + el + '_' + str(d[2]) + '.pdf')
        """

        if PLOT_DMDES:
            fig, ax = plt.subplots()
            ax.scatter(e/1e17, log_dm_de, c='C0', rasterized=True,
            	alpha=0.5, s=1, edgecolors='none')
            ax.plot(e/1e17, slog_dm_de, c='C1')
            ax.set_xlim(-10, 10)
            ax.set_ylim(9, 18)
            ax.set_xlabel(r'$E\ \mathrm{[10^{17}\ erg\ g^{-1}]}$')
            ax.set_ylabel(r'$\log\ dM/dE\ \mathrm{[g^2\ erg^{-1}]}$')
            ax.grid()
            fig.tight_layout(pad=0.3)
            if args.cluster == 'hyades': directory = '/pfs/lawsmith/results/dmdes/gaussian_wrap/' + d[0].replace(".","_")
            elif args.cluster == 'fend': directory = '/groups/dark/lawsmith/results/'+d[0].replace(".","_")
            if not os.path.exists(directory): os.makedirs(directory)
            fig.savefig(directory+'/dmde_'+d[0].replace(".","_")+'_'+d[1]+'_'+el+'_'+str(d[2])+'.pdf')

        #ax2.set_ylim(-6, 1)
        ax2.set_xlim(-2, 2)
        ax2.set_ylabel(r'$\log\ \dot M\ {\rm [M_\odot/yr]}$')
        ax2.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
        fig2.tight_layout(pad=0.3)
        if args.cluster == 'hyades': directory = '/pfs/lawsmith/results/dmdts/final'
        elif args.cluster == 'fend': directory = '/groups/dark/lawsmith/results/'+d[0].replace(".","_")
        fig2.savefig(directory+'/dmdt_'+d[0].replace(".","_")+'_'+d[1]+'_'+el+'_'+str(d[2])+'.pdf')

        plt.close('all')
