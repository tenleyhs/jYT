"""
trying to combine all files with plot_ having to read in histograms
TODO in progress
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.io import ascii
execfile('/pfs/lawsmith/jYT/my_settings.py')
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['font.size'] = 18

LW = 1.5
LWD = 2.0

ds = [
    # run,                  chk
    ['m1.0_p1_b1.0',        '0040'],
	['m1.0_p1_b1.5',		'0040'],
	['m1.0_p1_b1.75',		'0040'],
	['m1.0_p1_b2.0',		'0040'],
	['m1.0_p1_b3.0',		'0040'],
	['m1.0_p10_b1.0',		'0040'],
	['m1.0_p10_b1.0_256',	'0040'],
	['m1.0_p10_b1.5',		'0040'],
	['m1.0_p10_b2.0',		'0040'],
	['m1.0_p10_b2.0_256',	'0040'],
	['m1.0_p10_b2.5',		'0040'],
	['m1.0_p10_b3.0',		'0040'],
	['m1.0_p10_b4.0',		'0040'],
	['m1.0_p10_b5.0',		'0040'],
	['m1.0_p16_b1.0',		'0040'],
	['m1.0_p16_b1.5',		'0040'],
	['m1.0_p16_b2.0',		'0040'],
	['m1.0_p16_b3.0',		'0040'],
	['m1.0_p16_b4.0',		'0040'],
	['m1.0_p16_b5.0',		'0040'],
]

els = ['ev', 'h1', 'he4', 'o16', 'c12', 'ne20', 'n14']







### for total mdot, not split into elements
def myplot_ev(ds=None, text=None, savename=None):
	el = 'ev'
	fig, ax = plt.subplots()
	for d in ds:
		# log_t_yr, log_mdot_moyr
		x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/'
				+ d[0].replace(".","_") + '_' + d[1] + '_' + el + '.dat', skiprows=1)
		ax.plot(log_t_yr, log_mdot_moyr, lw=LW, label=d[2])

	ax.set_xlim(-2, 0.25)
	ax.set_ylim(-2, 1)
	ax.set_ylabel(r'$\log\ \dot M\ {\rm [M_\odot/yr]}$')
	ax.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
	ax.legend(loc=1)
	ax.text(-1.95,0.8,text)
	fig.tight_layout()
	fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/paper/' + savename)
	plt.close('all')


### for same age different betas
# p1
ds = [
    ['m1.0_p1_b1.0',        '0040',		r'$\beta=1.0$'],
	['m1.0_p1_b1.5',		'0040',		r'$\beta=1.5$'],
	['m1.0_p1_b1.75',		'0040',		r'$\beta=1.75$'],
	['m1.0_p1_b2.0',		'0040',		r'$\beta=2.0$'],
	['m1.0_p1_b3.0',		'0040',		r'$\beta=3.0$'],
]
text = 'age=0Gyr'
savename = 'mdot_p1_allbetas.pdf'
myplot_ev(ds=ds, text=text, savename=savename)

# p10
ds = [
    ['m1.0_p10_b1.0',       '0040',		r'$\beta=1.0$'],
	['m1.0_p10_b1.0_256',	'0040',		r'$\beta=1.0 (256)$'],
	['m1.0_p10_b1.5',		'0040',		r'$\beta=1.5$'],
	['m1.0_p10_b2.0',		'0040',		r'$\beta=2.0$'],
	['m1.0_p10_b2.0_256',	'0040',		r'$\beta=2.0 (256)$'],
	['m1.0_p10_b2.5',		'0040',		r'$\beta=2.5$'],
	['m1.0_p10_b3.0',		'0040',		r'$\beta=3.0$'],
	['m1.0_p10_b4.0',		'0040',		r'$\beta=4.0$'],
	['m1.0_p10_b5.0',		'0040',		r'$\beta=5.0$'],
]
text = 'age=4.8Gyr'
savename = 'mdot_p10_allbetas.pdf'
myplot_ev(ds=ds, text=text, savename=savename)

# p16
ds = [
    ['m1.0_p16_b1.0',       '0040',		r'$\beta=1.0$'],
	['m1.0_p16_b1.5',		'0040',		r'$\beta=1.5$'],
	['m1.0_p16_b2.0',		'0040',		r'$\beta=2.0$'],
	['m1.0_p16_b3.0',		'0040',		r'$\beta=3.0$'],
	['m1.0_p16_b4.0',		'0040',		r'$\beta=4.0$'],
	['m1.0_p16_b5.0',		'0040',		r'$\beta=5.0$'],
]
text = 'age=8.4Gyr'
savename = 'mdot_p16_allbetas.pdf'
myplot_ev(ds=ds, text=text, savename=savename)



### for same beta different ages
# TODO this is where I left off


# Guillochon 2013 dmdts
g13_43 = np.loadtxt('/pfs/lawsmith/dmdts/4-3/' + beta + '.dat')
g13_53 = np.loadtxt('/pfs/lawsmith/dmdts/5-3/' + beta + '.dat')

ax2.plot(g13_43[:,0], g13_43[:,1], lw=lwd, label='GRR13 4/3', ls=':', c='k')
#ax.plot(g13_53[:,0], g13_53[:,1], lw=lw, label='GRR13 5/3', ls='-.', c='k')

# Lodato/Monica dmdt
"""
m17_53 = np.loadtxt('/pfs/lawsmith/dmdts/poly53-b' + beta2 + 'dmdT.data')
ax2.plot(np.log10(m17_53[:,1]/yr), np.log10(m17_53[:,0]*yr/M_sun), lw=lw, label='analytic 5/3', color='gray')
"""

if beta2 == '1.0':
    m17_43 = np.loadtxt('/pfs/lawsmith/dmdts/poly_dmdt.data', skiprows=1)
    ax2.plot(np.log10(m17_43[:,5]/yr), np.log10(m17_43[:,4]*yr/M_sun), lw=lwd, label='analytic 4/3', color='gray', ls=':')


# for combined plot
beta = '2.000'
beta2 = '2.0'
#text = r'$\beta=2.0$'
ds = [
    # run,                  chk,        label,          wlen,   log_t_cutoff
    # these lines are for same beta, different age
    # b1.0
    ['m1.0_p1_b1.0',        '0100',     'age=0Gyr',     300,    0.],
    ['m1.0_p10_b1.0',       '0100',     'age=4.8Gyr',   400,    -0.4],
    ['m1.0_p10_b1.0_256',   '0100',     'age=4.8Gyr',   300,    -0.4],
    ['m1.0_p16_b1.0_pfe',   '0052',     'age=8.4Gyr',   500,    1],
    # b2.0
    'm1.0_p1_b2.0',         '0080',     'age=0Gyr',     250,    -0.1],
    'm1.0_p10_b2.0',        '0075',     'age=4.8Gyr',   200,    0.25],
    'm1.0_p16_b2.0',        '0075',     'age=8.4Gyr',   300,    3],
    # b3.0
    ['m1.0_p1_b3.0',        '0060',     'age=0Gyr',     250,    0.3],
    ['m1.0_p10_b3.0',       '0060',     'age=4.8Gyr',   200,    0.4],
    ['m1.0_p16_b3.0',       '0060',     'age=8.4Gyr',   500,    0],
]

# b1.0
ds = [
    ['m1.0_p1_b1.0',        '0040',     'age=0Gyr'],
    ['m1.0_p10_b1.0',       '0040',     'age=4.8Gyr'],
	['m1.0_p16_b1.0',       '0040',     'age=8.4Gyr'],
]
text = r'$\beta=1.0$'
savename = 'mdot_b1.0_allages.pdf'
myplot_ev(ds=ds, text=text, savename=savename)













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
