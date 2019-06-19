"""
trying to combine all files with plot_ having to read in histograms
maybe a faster/better way would be to load all the files first, into some
dictionary or something. because right now repeats.
"""
## TODO need to shift dmdes


import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.io import ascii
import argparse, sys

parser = argparse.ArgumentParser()
parser.add_argument('--cluster', help='cluster, e.g., hyades/fend/pfe', type=str, default='fend')
args = parser.parse_args()

if args.cluster == 'hyades':
	execfile('/pfs/lawsmith/jYT/my_settings.py')
	path = '/pfs/lawsmith/FLASH4.3/runs/results/'   # + 'dmdts/data/'
elif args.cluster == 'fend':
	execfile('/groups/dark/lawsmith/jYT/my_settings.py')
	path = '/groups/dark/lawsmith/results/'

# option for seaborn style background
#mpl.style.use('seaborn')

LW = 2
LWD = 2.0
S = 30
S2 = 200


### for total mdot, not split into elements
def myplot_ev(ds=None, text=None, savename=None,
                x2=None, y2=None, l2=None,
                x3=None, y3=None, l3=None):
    el = 'ev'
    fig, ax = plt.subplots()
    for d in ds:
        # log_t_yr, log_mdot_moyr
        x, y = np.loadtxt(path + d[0].replace(".","_") + '/'
                + d[0].replace(".","_") + '_' + d[1] + '_' + el + '.dat',
                skiprows=1, unpack=True)
        ax.plot(x[x>d[3]], y[x>d[3]], lw=LW, label=d[2])

    	#if beta == '1.000':
        #	print x[np.argmax(y)], max(y)

    if x2 is not None:
        ax.plot(x2, y2, lw=LWD, label=l2, ls=':', c='k')
    if x3 is not None:
        ax.plot(x3, y3, lw=LWD, label=l3, ls=':', c='gray')

    ax.set_xlim(-2, 0.25)
    ax.set_ylim(-2, 1)
    if savename == 'mdot_1Msun_vs_3Msun.pdf':
        ax.set_ylim(-2, 1.5)
    ax.set_ylabel(r'$\log\ \dot M\ {\rm [M_\odot/yr]}$')
    ax.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
    ax.legend(loc=1)
    ax.text(-1.95,0.75,text)
    ax.grid()
    fig.tight_layout(pad=0.3)
    fig.savefig(path + 'paper/' + savename)
    plt.close('all')


"""
### for same age different betas
# p1
ds = [
    ['m1.0_p1_b1.0',        '0100',		r'$\beta=1.0$'],
	['m1.0_p1_b1.5',		'0090',		r'$\beta=1.5$'],
	#['m1.0_p1_b1.75',		'0040',		r'$\beta=1.75$'],
	['m1.0_p1_b2.0',		'0080',		r'$\beta=2.0$'],
	#['m1.0_p1_b3.0',		'0040',		r'$\beta=3.0$'],
]
text = 'age=0Gyr'
savename = 'mdot_p1_allbetas.pdf'
myplot_ev(ds=ds, text=text, savename=savename)

# p10
ds = [
    ['m1.0_p10_b1.0',       '0100',		r'$\beta=1.0$'],
	#['m1.0_p10_b1.0_256',	'0040',		r'$\beta=1.0 (256)$'],
	#['m1.0_p10_b1.5',		'0040',		r'$\beta=1.5$'],
	['m1.0_p10_b2.0',		'0075',		r'$\beta=2.0$'],
	#['m1.0_p10_b2.0_256',	'0040',		r'$\beta=2.0 (256)$'],
	#['m1.0_p10_b2.5',		'0040',		r'$\beta=2.5$'],
	['m1.0_p10_b3.0',		'0060',		r'$\beta=3.0$'],
]
text = 'age=4.8Gyr'
savename = 'mdot_p10_allbetas.pdf'
myplot_ev(ds=ds, text=text, savename=savename)

# p16
ds = [
    ['m1.0_p16_b1.0',       '0100',		r'$\beta=1.0$'],
	#['m1.0_p16_b1.5',		'0040',		r'$\beta=1.5$'],
	['m1.0_p16_b2.0',		'0075',		r'$\beta=2.0$'],
	['m1.0_p16_b3.0',		'0060',		r'$\beta=3.0$'],
	['m1.0_p16_b4.0',		'0050',		r'$\beta=4.0$'],
]
text = 'age=8.4Gyr'
savename = 'mdot_p16_allbetas.pdf'
myplot_ev(ds=ds, text=text, savename=savename)
"""


### for same beta different ages
# b1.0
# Guillochon 2013 dmdts
beta = '1.000'
g13_43 = np.loadtxt('/groups/dark/lawsmith/Guillochon2013_dmdts/4-3/' + beta + '.dat')
#g13_53 = np.loadtxt('/pfs/lawsmith/dmdts/5-3/' + beta + '.dat')
# Lodato/Monica dmdts
# only have 4/3 for beta = 1.0. can ask Monica or do myself for other betas
m18_43 = np.loadtxt('/groups/dark/lawsmith/Gallegos2018/poly_dmdt.data', skiprows=1)
#m17_53 = np.loadtxt('/pfs/lawsmith/dmdts/poly53-b' + beta2 + 'dmdT.data')
ds = [
    ['m1.0_p1_b1.0',        '0060',     'age=0Gyr',     -1.5],
    ['m1.0_p10_b1.0',       '0060',     'age=4.8Gyr',   -1.5],
	['m1.0_p16_b1.0',       '0060',     'age=8.4Gyr',   -4],
]
text = r'$\beta=1.0$'
savename = 'mdot_b1_0_allages.pdf'
myplot_ev(ds=ds, text=text, savename=savename,
    # scale GRR2013 to ZAMS sun
	x2=g13_43[:,0] + 1.5*np.log10(0.9012), y2=g13_43[:,1] - 1.5*np.log10(0.9012), l2='GRR13 4/3')
	#x3=np.log10(m18_43[:,5]/yr), y3=np.log10(m18_43[:,4]*yr/M_sun), l3='analytic 4/3')
#ax.plot(g13_53[:,0], g13_53[:,1], lw=lw, label='GRR13 5/3', ls='-.', c='k')
#ax2.plot(np.log10(m18_53[:,1]/yr), np.log10(m18_53[:,0]*yr/M_sun), lw=lw, label='analytic 5/3', color='gray')

# b2.0
beta = '2.000'
g13_43 = np.loadtxt('/groups/dark/lawsmith/Guillochon2013_dmdts/4-3/' + beta + '.dat')
ds = [
    ['m1.0_p1_b2.0',        '0080',     'age=0Gyr',     -4],
    ['m1.0_p10_b2.0',       '0075',     'age=4.8Gyr',   -1.82],
	['m1.0_p16_b2.0',       '0075',     'age=8.4Gyr',   -1.72],
]
text = r'$\beta=2.0$'
savename = 'mdot_b2_0_allages.pdf'
myplot_ev(ds=ds, text=text, savename=savename,
    # scale GRR2013 to ZAMS sun
    x2=g13_43[:,0] + 1.5*np.log10(0.9012), y2=g13_43[:,1] - 1.5*np.log10(0.9012), l2='GRR13 4/3')

# b3.0
beta = '3.000'
g13_43 = np.loadtxt('/groups/dark/lawsmith/Guillochon2013_dmdts/4-3/' + beta + '.dat')
ds = [
    ['m1.0_p1_b3.0',        '0060',     'age=0Gyr',     -1.9],
    ['m1.0_p10_b3.0',       '0060',     'age=4.8Gyr',   -1.9],
	#['m1.0_p16_b3.0',       '0060',     'age=8.4Gyr',   -1.85],
    ['m1.0_p16_b3.0_96k_1000',       '0060',     'age=8.4Gyr',   -1.85],
]
text = r'$\beta=3.0$'
savename = 'mdot_b3_0_allages.pdf'
myplot_ev(ds=ds, text=text, savename=savename,
    # scale GRR2013 to ZAMS sun
    x2=g13_43[:,0] + 1.5*np.log10(0.9012), y2=g13_43[:,1] - 1.5*np.log10(0.9012), l2='GRR13 4/3')


# comparing 1Msun and 3Msun
ds = [
    ['m1.0_p1_b2.0',        '0080',     'ZAMS '+ r'$1M_\odot, \beta=2$',     -4],
    #['m1.0_p10_b3.0',       '0060',    'middle-age '+ r'$1M_\odot, \beta=3$',   -1.9],
    ['m1.0_p16_b4.0',       '0050',     'TAMS '+ r'$1M_\odot, \beta=4$',   -4],
    ['m3.0_p1_b2.0_48k',       '0070',     'ZAMS '+ r'$3M_\odot, \beta=2$',   -4],
    ['m3.0_p16_b4.0_48k',       '0050',     'TAMS '+ r'$3M_\odot, \beta=4$',   -4],
]
text = ''
savename = 'mdot_1Msun_vs_3Msun.pdf'
myplot_ev(ds=ds, text=text, savename=savename)


### for delta_m vs beta
# in plot_deltaM_vs_beta for now



### for line ratios (plot_comp_solar)
f_h1_ZAMS = 0.7153
f_he4_ZAMS = 0.2704
f_c12_ZAMS = 2.435e-3
f_n14_ZAMS = 7.509e-4
f_o16_ZAMS = 6.072e-3
f_ne20_ZAMS = 1.232e-3
def myplot_lr(ds):
    els = ['ev', 'h1', 'he4', 'c12', 'n14', 'o16']
    colors = ['', '', 'C0', 'C1', 'C2', 'C3']
    for d in ds:
        fig, ax = plt.subplots()
        for el, col in zip(els, colors):
            # log_t_yr, log_mdot_moyr
            x, y = np.loadtxt(path + d[0].replace(".","_") + '/'
            + d[0].replace(".","_") + '_' + d[1] + '_' + el + '.dat',
            skiprows=1, unpack=True)
            if el == 'ev':
                tpeak = 10**x[np.argmax(y)]
                continue
            elif el == 'h1':
                h1_mdot = 10**y
                continue
            elif el == 'he4':
                el_lr = (10**y/h1_mdot)/(f_he4_ZAMS/f_h1_ZAMS)
            elif el == 'o16':
                #if d[0] == 'm1.0_p16_b4.0': continue
                el_lr = (10**y/h1_mdot)/(f_o16_ZAMS/f_h1_ZAMS)
            elif el == 'c12':
                el_lr = (10**y/h1_mdot)/(f_c12_ZAMS/f_h1_ZAMS)
            elif el == 'n14':
                el_lr = (10**y/h1_mdot)/(f_n14_ZAMS/f_h1_ZAMS)

            ax.plot(10**x/tpeak, el_lr, lw=LW, label=el, color=col)
            #ax.plot(np.log10(10**x/tpeak), el_lr, lw=LW, label=el, color=col)

        if d[0] == 'm1.0_p10_b2.0':
            monica_name = 'M1.0-p10dmdT.data'
            ax.set_ylim(0, 3)
        if d[0] == 'm1.0_p10_b3.0':
            monica_name = 'M1.0-p10dmdT.data'
            ax.set_ylim(0, 3)
            #ax.text(-0.38, 3.72, d[2])

        if d[0] == 'm1.0_p16_b1.0':
            monica_name = 'M1.0-p16dmdT.data'
            ax.set_ylim(0, 4)
        if d[0] == 'm1.0_p16_b2.0':
            monica_name = 'M1.0-p16dmdT.data'
            ax.set_ylim(0, 4)
        if d[0] == 'm1.0_p16_b3.0':
            monica_name = 'M1.0-p16dmdT.data'
            ax.set_ylim(0, 4)

        if d[0] == 'm1.0_p16_b3.0_96k_1000':
            monica_name = 'M1.0-p16dmdT.data'
            ax.set_ylim(0, 4)

        if d[0] == 'm1.0_p16_b4.0':
            monica_name = 'M1.0-p16dmdT.data'
            ax.set_ylim(0, 4)
            #ax.text(-0.38, 6.5, d[2])

        if d[0] == 'm3.0_p1_b2.0_48k':
            monica_name = 'M3.0-p1dmdT.data'
            ax.set_ylim(0, 5)
        if d[0] == 'm3.0_p16_b4.0_48k':
            monica_name = 'M3.0-p16dmdT.data'
            ax.set_ylim(0, 5)
            #ax.text(-0.38, 6.5, d[2])

        ax.set_title(d[2])

        dat = np.genfromtxt('/groups/dark/lawsmith/Gallegos2018/Archive/' + monica_name,
            skip_header=1, usecols=(0,1,2,4,6,10,12),
            dtype=[('dmdt',float), ('t',float), ('dmdt_he4',float), ('dmdt_o16',float),
            ('dmdt_n14',float), ('dmdt_h1',float), ('dmdt_c12',float)])
        t_peak = dat['t'][np.argmax(dat['dmdt'])]
        ax.plot(dat['t']/t_peak, (dat['dmdt_he4']/dat['dmdt_h1'])/(f_he4_ZAMS/f_h1_ZAMS), ls='--', color='C0')
        ax.plot(dat['t']/t_peak, (dat['dmdt_c12']/dat['dmdt_h1'])/(f_c12_ZAMS/f_h1_ZAMS), ls='--', color='C1')
        ax.plot(dat['t']/t_peak, (dat['dmdt_n14']/dat['dmdt_h1'])/(f_n14_ZAMS/f_h1_ZAMS), ls='--', color='C2')
        ax.plot(dat['t']/t_peak, (dat['dmdt_o16']/dat['dmdt_h1'])/(f_o16_ZAMS/f_h1_ZAMS), ls='--', color='C3')
        ax.plot([98,99],[98,99],ls='--',label='analytic',color='C7')
        ax.plot([98,99],[98,99],label='simulation',color='C7')

        ax.set_xlim(0, 5)

        ax.grid()
        #ax.set_xlabel(r'$\log\ t/t_{\rm peak}$')
        ax.set_xlabel(r'$t/t_{\rm peak}$')
        ax.set_ylabel(r'$X/X_\odot$')

        ax.legend(loc=1)
        fig.tight_layout(pad=0.3)
        fig.savefig(path + 'paper/mdot_comp_solar_' + d[0].replace(".","_") + '_' + d[1] + '.pdf')
        plt.close('all')

ds = [
    #['m1.0_p1_b1.0',       '0040',		'age=0Gyr, '+r'$\beta=1.0$'],
	#['m1.0_p1_b1.5',		'0040',		'age=0Gyr, '+r'$\beta=1.5$'],
	#['m1.0_p1_b1.75',		'0040',		'age=0Gyr, '+r'$\beta=1.75$'],
	#['m1.0_p1_b2.0',		'0040',		'age=0Gyr, '+r'$\beta=2.0$'],
	#['m1.0_p1_b3.0',		'0040',		'age=0Gyr, '+r'$\beta=3.0$'],
	#['m1.0_p10_b1.0',		'0040',		'age=4.8Gyr, '+r'$\beta=1.0$'],
	#['m1.0_p10_b1.0_256',	'0040',		'age=4.8Gyr, '+r'$\beta=1.0$'],
	#['m1.0_p10_b1.5',		'0040',		'age=4.8Gyr, '+r'$\beta=1.5$'],
	['m1.0_p10_b2.0',		'0075',		r'$1M_\odot$'+', age=4.8Gyr, '+r'$\beta=2.0$'],
	#['m1.0_p10_b2.5',		'0040',		'age=4.8Gyr, '+r'$\beta=2.5$'],
	#['m1.0_p10_b3.0',		'0040',		'age=4.8Gyr, '+r'$\beta=3.0$'],
    ['m1.0_p10_b3.0',		'0060',		r'$1M_\odot$'+', age=4.8Gyr, '+r'$\beta=3.0$'],
	#['m1.0_p10_b4.0',		'0040',		'age=4.8Gyr, '+r'$\beta=4.0$'],
	#['m1.0_p10_b5.0',		'0040',		'age=4.8Gyr, '+r'$\beta=5.0$'],
	['m1.0_p16_b1.0',		'0060',		r'$1M_\odot$'+', age=8.4Gyr, '+r'$\beta=1.0$'],
	#['m1.0_p16_b1.5',		'0040',		'age=8.4Gyr, '+r'$\beta=1.5$'],
	['m1.0_p16_b2.0',		'0075',		r'$1M_\odot$'+', age=8.4Gyr, '+r'$\beta=2.0$'],
    #['m1.0_p16_b2.0',		'0075',		'age=8.4Gyr, '+r'$\beta=2.0$'],
    ['m1.0_p16_b3.0',		'0060',		r'$1M_\odot$'+', age=8.4Gyr, '+r'$\beta=3.0$'],
    ['m1.0_p16_b3.0_96k_1000',       '0060',     r'$1M_\odot$'+', age=8.4Gyr, '+r'$\beta=3.0$'],
    ['m1.0_p16_b4.0',		'0050',		r'$1M_\odot$'+', age=8.4Gyr, '+r'$\beta=4.0$'],
	#['m1.0_p16_b5.0',		'0040',		'age=8.4Gyr, '+r'$\beta=5.0$'],

	['m3.0_p1_b2.0_48k',		'0070',		r'$3M_\odot$'+', age=0Gyr, '+r'$\beta=2.0$'],
	['m3.0_p16_b4.0_48k',		'0050',		r'$3M_\odot$'+', age=0.3Gyr, '+r'$\beta=4.0$'],
]
myplot_lr(ds)


"""
### for summary figure
# 1st way of doing it. my parameter grid
mpl.rcParams['ytick.right'] = False
fig, ax = plt.subplots()
dss = [
        [
    ['m1.0_p1_b1.0',		'0040',     1.0],
	['m1.0_p1_b1.5',		'0040',		1.5],
	['m1.0_p1_b1.75',		'0040',		1.75],
	['m1.0_p1_b2.0',		'0040',		2.0],
	['m1.0_p1_b3.0',		'0040',		3.0]
        ],
        [
	['m1.0_p10_b1.0',		'0040',		1.0],
	#['m1.0_p10_b1.0_256',	'0040',		1.0],
	['m1.0_p10_b1.5',		'0040',		1.5],
	['m1.0_p10_b2.0',		'0040',		2.0],
	#['m1.0_p10_b2.0_256',	'0040',		2.0],
	['m1.0_p10_b2.5',		'0040',		2.5],
	['m1.0_p10_b3.0',		'0040',		3.0],
	['m1.0_p10_b4.0',		'0040',		4.0],
	#['m1.0_p10_b5.0',		'0040',		5.0]
        ],
        [
	['m1.0_p16_b1.0',		'0040',		1.0],
	['m1.0_p16_b1.5',		'0040',		1.5],
	['m1.0_p16_b2.0',		'0040',		2.0],
	['m1.0_p16_b3.0',		'0040',		3.0],
	['m1.0_p16_b4.0',		'0040',		4.0],
	#['m1.0_p16_b5.0',		'0040',		5.0]
        ]
]
labels = ['age=0Gyr', 'age=4.8Gyr', 'age=8.4Gyr']
carr = ['C0', 'C1', 'C2']
f_h1_ZAMS = 0.7153
f_he4_ZAMS = 0.2704
f_c12_ZAMS = 2.435e-3
f_n14_ZAMS = 7.509e-4
f_o16_ZAMS = 6.072e-3
f_ne20_ZAMS = 1.232e-3
els = ['ev', 'h1', 'n14']
age_array = []
beta_array = []
g_array = []
for i, ds in enumerate(dss):
	for d in ds:
		for el in els:
			# log_t_yr, log_mdot_moyr
			x, y = np.loadtxt(path + d[0].replace(".","_") + '/'
					+ d[0].replace(".","_") '_' + d[1] + '_' + el + '_20.dat',
					skiprows=1, unpack=True)
			if el == 'ev':
				tpeak = 10**x[np.argmax(y)]
				continue
			elif el == 'h1':
			    h1_mdot = 10**y
			    continue
			elif el == 'n14':
				el_lr = (10**y/h1_mdot)/(f_n14_ZAMS/f_h1_ZAMS)
				if i == 0: age_array.append(0.)
				if i == 1: age_array.append(4.8)
				if i == 2: age_array.append(8.4)
				beta_array.append(d[2])
				# TODO can change definition of g
				g_array.append(el_lr[(np.abs(np.log10(10**x/tpeak) - 1.0)).argmin()])

# one way
#cbaxes = fig.add_axes([0.4, 0.5, 0.5, 0.03])
#sc = ax.scatter(beta_array, age_array, c=g_array, cmap='viridis', s=S2)
#cb = fig.colorbar(sc, cax=cbaxes, orientation='horizontal') #ticks=[0.8,0.9,1]
# another way
sc = ax.scatter(beta_array, age_array, c=g_array, cmap='viridis', s=S2)
cb = fig.colorbar(sc)
#cb.set_label(label=r'$(X/X_\odot)_{^{14}{\rm N}}$', fontsize=22)
cb.set_label(label=r'$g$')#, fontsize=22)
ax.set_xlabel(r'$\beta$')
ax.set_ylabel('age [Gyr]')
fig.tight_layout(pad=0.3)
fig.savefig(path + 'paper/mdot_comp_solar_summary1.pdf')
plt.close('all')
mpl.rcParams['ytick.right'] = True



# 2nd way of doing it. just plotting N14, a few betas, comparing ages
fig, ax = plt.subplots()
dss = [
# no p1 because straight lines
        [
	#['m1.0_p10_b1.0',		1.0],
	['m1.0_p10_b2.0_256',	2.0,    '0080'],
	['m1.0_p10_b3.0_256',	3.0,    '0057'],
    ['m1.0_p10_b4.0',		4.0,    '0050'],
        ],
        [
	#['m1.0_p16_b1.0',		1.0],
	['m1.0_p16_b2.0',		2.0,    '0075'],
	['m1.0_p16_b3.0',		3.0,    '0060'],
    ['m1.0_p16_b4.0',		4.0,    '0050'],
        ]
]
labels = ['age=4.8Gyr', 'age=8.4Gyr']
carr = ['C0', 'C1']
lss = ['-', ':', '-.']
lws = [1.5, 2.5, 2]
f_h1_ZAMS = 0.7153
f_he4_ZAMS = 0.2704
f_c12_ZAMS = 2.435e-3
f_n14_ZAMS = 7.509e-4
f_o16_ZAMS = 6.072e-3
f_ne20_ZAMS = 1.232e-3
els = ['ev', 'h1', 'n14']
for i, ds in enumerate(dss):
	for k, d in enumerate(ds):
		for el in els:
			# log_t_yr, log_mdot_moyr
			x, y = np.loadtxt(path + d[0].replace(".","_") + '/'
					+ d[0].replace(".","_") + '_' + d[2] + '_' + el + '_20.dat',
					skiprows=1, unpack=True)
			if el == 'ev':
				tpeak = 10**x[np.argmax(y)]
				continue
			elif el == 'h1':
			    h1_mdot = 10**y
			    continue
			elif el == 'n14':
				el_lr = (10**y/h1_mdot)/(f_n14_ZAMS/f_h1_ZAMS)
                if k == 0:
                    label = labels[i]
                else:
                    label = None
                ax.plot(np.log10(10**x/tpeak), el_lr, ls=lss[k], lw=lws[k],
                        c=carr[i], label=label)

ax.plot(-9,-9, ls=lss[0], c='k', label=r'$\beta=2.0$')
ax.plot(-9,-9, ls=lss[1], c='k', label=r'$\beta=3.0$')
ax.plot(-9,-9, ls=lss[2], c='k', label=r'$\beta=4.0$')
ax.set_xlim(-0.3, 3)
ax.set_ylim(0, 4)
ax.set_xlabel(r'$\log\ t/t_{\rm peak}$')
ax.set_ylabel(r'$\left(X/X_\odot\right)_{^{14}{\rm N}}$')
ax.legend(loc=2)
fig.tight_layout(pad=0.3)
fig.savefig(path + 'paper/mdot_comp_solar_summary2.pdf')
plt.close('all')




# 3rd way of doing it. like Monica's figure
fig, ax = plt.subplots()
dss = [
        [
    ['m1.0_p1_b1.0',		'0040',     1.0],
	['m1.0_p1_b1.5',		'0040',		1.5],
	#['m1.0_p1_b1.75',		'0040',		1.75],
	['m1.0_p1_b2.0',		'0040',		2.0],
	#['m1.0_p1_b3.0',		'0040',		3.0]
        ],
        [
	['m1.0_p10_b1.0',		'0040',		1.0],
	#['m1.0_p10_b1.0_256',	'0040',		1.0],
	#['m1.0_p10_b1.5',		'0040',		1.5],
	['m1.0_p10_b2.0',		'0040',		2.0],
	#['m1.0_p10_b2.0_256',	'0040',		2.0],
	#['m1.0_p10_b2.5',		'0040',		2.5],
	['m1.0_p10_b3.0',		'0040',		3.0],
	#['m1.0_p10_b4.0',		'0040',		4.0],
        ],
        [
	['m1.0_p16_b1.0',		'0040',		1.0],
	#['m1.0_p16_b1.5',		'0040',		1.5],
	['m1.0_p16_b2.0',		'0040',		2.0],
	#['m1.0_p16_b3.0',		'0040',		3.0],
	['m1.0_p16_b4.0',		'0040',		4.0],
        ]
]


labels = ['age=0Gyr', 'age=4.8Gyr', 'age=8.4Gyr']
carr = ['C0', 'C1', 'C2']
f_h1_ZAMS = 0.7153
f_he4_ZAMS = 0.2704
f_c12_ZAMS = 2.435e-3
f_n14_ZAMS = 7.509e-4
f_o16_ZAMS = 6.072e-3
f_ne20_ZAMS = 1.232e-3
els = ['ev', 'h1', 'n14']
tpeak_array = []
mdotpeak_array = []
g_array = []
for i, ds in enumerate(dss):
    tp = []
    mp = []
    for d in ds:
        for el in els:
            # log_t_yr, log_mdot_moyr
            x, y = np.loadtxt(path + d[0].replace(".","_") + '/'
                + d[0].replace(".","_") '_' + d[1] + '_' el + '.dat',
                skiprows=1, unpack=True)
            if el == 'ev':
                tpeak = 10**x[np.argmax(y)]
                continue
            elif el == 'h1':
                h1_mdot = 10**y
                continue
            elif el == 'n14':
                el_lr = (10**y/h1_mdot)/(f_n14_ZAMS/f_h1_ZAMS)
                tpeak_array.append(np.log10(tpeak))
                mdotpeak_array.append(max(y))
                # TODO change 1.0 value
                g_array.append(el_lr[(np.abs(np.log10(10**x/tpeak) - 0.5)).argmin()])
                tp.append(np.log10(tpeak))
                mp.append(max(y))
    ax.plot(tp, mp, label=labels[i])

# one way
#cbaxes = fig.add_axes([0.4, 0.5, 0.5, 0.03])
#sc = ax.scatter(beta_array, age_array, c=g_array, cmap='viridis', s=S2)
#cb = fig.colorbar(sc, cax=cbaxes, orientation='horizontal') #ticks=[0.8,0.9,1]
# another way
sc = ax.scatter(tpeak_array, mdotpeak_array, c=g_array, cmap='viridis', s=S2)
cb = fig.colorbar(sc)
#cb.set_label(label=r'$(X/X_\odot)_{^{14}{\rm N}}$', fontsize=22)
cb.set_label(label=r'$g$')#, fontsize=22)
ax.set_xlabel(r'$\log\ t_{\rm peak}\ {\rm [yr]}$')
ax.set_ylabel(r'$\log\ \dot M_{\rm peak}\ {\rm [M_\odot/yr]}$')
ax.legend(loc=1)
fig.tight_layout(pad=0.3)
fig.savefig(path + 'paper/mdot_comp_solar_summary3_05.pdf')
plt.close('all')



### for power law index, tpeak vs beta, mdotpeak vs beta
dss = [
        [
    ['m1.0_p1_b1.0',        1.0,	'0100'],
	['m1.0_p1_b1.5',		1.5,	'0090'],
	['m1.0_p1_b1.75',		1.75,	'0090'],
	['m1.0_p1_b2.0',		2.0,	'0080'],
	['m1.0_p1_b3.0',		3.0,	'0060']
	#['m1.0_p1_b1.0',        1.0,	'0040'],
	#['m1.0_p1_b1.5',		1.5,	'0040'],
	#['m1.0_p1_b1.75',		1.75,	'0040'],
	#['m1.0_p1_b2.0',		2.0,	'0040'],
	#['m1.0_p1_b3.0',		3.0,	'0040']
        ],
        [
	['m1.0_p10_b1.0_256',	1.0,	'0040'],
	['m1.0_p10_b1.5',		1.5,	'0090'],
	['m1.0_p10_b2.0_256',	2.0,	'0080'],
	['m1.0_p10_b2.5',		2.5,	'0060'],
	['m1.0_p10_b3.0_256',	3.0,	'0057'],
	['m1.0_p10_b4.0',		4.0,	'0050'],
	#['m1.0_p10_b1.0',		1.0,	'0040'],
	#['m1.0_p10_b1.5',		1.5,	'0040'],
	#['m1.0_p10_b2.0',		2.0,	'0040'],
	#['m1.0_p10_b2.5',		2.5,	'0040'],
	#['m1.0_p10_b3.0',		3.0,	'0040'],
	#['m1.0_p10_b4.0',		4.0,	'0040'],
        ],
        [
	['m1.0_p16_b1.0',		1.0,	'0040'],
	['m1.0_p16_b1.5',		1.5,	'0090'],
	['m1.0_p16_b2.0',		2.0,	'0075'],
	['m1.0_p16_b3.0',		3.0,	'0060'],
	['m1.0_p16_b4.0',		4.0,	'0050'],
	#['m1.0_p16_b1.0',		1.0,	'0040'],
	#['m1.0_p16_b1.5',		1.5,	'0040'],
	#['m1.0_p16_b2.0',		2.0,	'0040'],
	#['m1.0_p16_b3.0',		3.0,	'0040'],
	#['m1.0_p16_b4.0',		4.0,	'0040'],
        ]
]
labels = ['age=0Gyr', 'age=4.8Gyr', 'age=8.4Gyr']
carr = ['C0', 'C1', 'C2']
fig, ax = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()
fig4, ax4 = plt.subplots()
fig5, ax5 = plt.subplots()
for i, ds in enumerate(dss):
	b_array = []
	ninf_array = []
	mdotpeak_array = []
	tpeak_array = []
	for k, d in enumerate(ds):
		x, y = np.loadtxt(path + d[0].replace(".","_") + '/'
					+ d[0].replace(".","_") + '_' + d[2] + '_ev_20.dat',
					skiprows=1, unpack=True)
		# mdotpeak, tpeak
		sel_y = np.isfinite(y)
		mdotpeak_array.append(max(y[sel_y]))
		tpeak_array.append(x[sel_y][np.argmax(y[sel_y])])

		# n infinity
		# take this at 1 yr. TODO can change this
		ix = (np.abs(x)).argmin() # log t = 0
		ninf = (y[ix-1] - y[ix]) / (x[ix-1] - x[ix])
		#ninf = (y[-2] - y[-1]) / (x[-2] - x[-1])
		b_array.append(d[1])
		ninf_array.append(ninf)

		# n instantaneous
		n = []
		for j in range(len(x)):
			if j == 0:
				continue
			s = (y[j] - y[j-1])/(x[j] - x[j-1])
			n.append(s)

		if k == 0:
			label = labels[i]
		else:
			label = None
		ax.plot(x[:-1], n, lw=LW, c=carr[i], alpha=0.4, label=label)

	ax2.plot(b_array, ninf_array, lw=LW, label=labels[i], alpha=0.5)
	ax2.scatter(b_array, ninf_array, s=S)

	# TODO add fitting function
	ax3.plot(b_array, 10**np.array(tpeak_array), lw=LW, label=labels[i], alpha=0.5)
	ax3.scatter(b_array, 10**np.array(tpeak_array), s=S)

	ax4.plot(b_array, mdotpeak_array, lw=LW, label=labels[i], alpha=0.5)
	ax4.scatter(b_array, mdotpeak_array, s=S)

	ax5.plot(10**np.array(tpeak_array), mdotpeak_array, label=labels[i], alpha=0.5)
	ax5.scatter(10**np.array(tpeak_array), mdotpeak_array, s=S)


ax.set_xlim(-2, 2)
ax.set_ylim(-4, 4)
ax.axhline(-5/3., ls=':', lw=LWD, c='k')
ax.set_ylabel(r'$n$')
ax.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
ax.legend()
fig.tight_layout(pad=0.3)
fig.savefig(path + 'paper/power_law_index_instantaneous.pdf')

ax2.set_xlim(0.9, 4.1)
ax2.axhline(-5/3., ls=':', lw=LWD, c='k')
ax2.set_ylabel(r'$n_\infty$')
ax2.set_xlabel(r'$\beta$')
ax2.legend()
fig2.tight_layout(pad=0.3)
fig2.savefig(path + 'paper/power_law_index_infinity.pdf')

def A_43(b):
    result = np.exp((27.261 - 27.516*b + 3.8716*b**2)/(1 - 3.2605*b - 1.3865*b**2))
    return result

def B_43(b):
    result = (-0.38670 + 0.57291*np.sqrt(b) - 0.31231*b)/(1 - 1.2744*np.sqrt(b) - 0.90053*b)
    return result

b43 = np.linspace(0.6, 4.0)
ax3.plot(b43, B_43(b43) * (0.9)**(3./2), lw=1.5, ls='--', label='GRR13 4/3, 0.9'+r'$R_\odot$', c='C0', alpha=0.5)
#ax3.plot(b43, B_43(b43), lw=1.5, ls='--', label='GRR13 4/3', c='C1', alpha=0.5)
#ax3.plot(b43, B_43(b43) * (1.29)**(3./2), lw=1.5, ls='--', label='GRR13 4/3', c='C2', alpha=0.5)
ax3.set_xlim(0.9, 4.1)
ax3.set_ylim(0.05, 0.19)
ax3.set_ylabel(r'$t_{\rm peak}\ \mathrm{[yr]}$')
ax3.set_xlabel(r'$\beta$')
ax3.legend()
fig3.tight_layout(pad=0.3)
fig3.savefig(path + 'paper/tpeak_vs_beta.pdf')

ax4.set_xlim(0.9, 4.1)
ax4.set_ylim(-1.2, 0.75)
ax4.plot(b43, np.log10(A_43(b43) * (0.9)**(-3./2)), lw=1.5, ls='--', label='GRR13 4/3, 0.9'+r'$R_\odot$', c='C0', alpha=0.5)
#ax4.plot(b43, np.log10(A_43(b43)), lw=1.5, ls='--', label='GRR13 4/3', c='C1', alpha=0.5)
#ax4.plot(b43, np.log10(A_43(b43) * (1.29)**(-3./2)), lw=1.5, ls='--', label='GRR13 4/3', c='C2', alpha=0.5)
ax4.set_ylabel(r'$\log\ \dot M_{\rm peak}\ {\rm [M_\odot/yr]}$')
ax4.set_xlabel(r'$\beta$')
ax4.legend()
fig4.tight_layout(pad=0.3)
fig4.savefig(path + 'paper/mdotpeak_vs_beta.pdf')


ax5.set_xlabel(r'$t_{\rm peak}\ \mathrm{[yr]}$')
ax5.set_ylabel(r'$\log\ \dot M_{\rm peak}\ {\rm [M_\odot/yr]}$')
ax5.legend()
fig5.tight_layout(pad=0.3)
fig5.savefig(path + 'paper/mdotpeak_vs_tpeak.pdf')
"""



"""
for old orange plot
bs = [
    '0.600',
    '0.650',
    '0.700',
    '0.750',
    '0.800',
    '0.850',
    '0.900',
    '1.000',
    '1.100',
    '1.200',
    '1.300',
    '1.400',
    '1.500',
    '1.600',
    '1.700',
    '1.800',
    '1.850',
    '1.900',
    '2.000',
    '2.500',
    '3.000',
    '3.500',
    '4.000'
    ]
for b in bs:
    #g13_43 = np.loadtxt('/nobackup/jlawsmit/dmdts/4-3/' + b + '.dat')
    g13_43 = np.loadtxt('/pfs/lawsmith/dmdts/4-3/' + b + '.dat')

    if b == beta:
        ax2.plot(g13_43[:,0], g13_43[:,1], ls='-', lw=1.5, color='orange', label='4/3')
    else:
        ax2.plot(g13_43[:,0], g13_43[:,1], ls='-', lw=0.5, color='orange')

    ax2.plot(g13_53[:,0], g13_53[:,1], ls=':', lw=2, color=blue, label='5/3')
"""
