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
from scipy.interpolate import UnivariateSpline
execfile('/pfs/lawsmith/jYT/my_settings.py')
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['font.size'] = 18

# for dmde's
def smooth(x,window_len=11,window='hanning'):
	if x.ndim != 1:
		raise ValueError, "smooth only accepts 1 dimension arrays."
	if x.size < window_len:
		raise ValueError, "Input vector needs to be bigger than window size."
	if window_len<3:
		return x
	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError, "Window is one of the above"
	s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
	if window == 'flat': #moving average
		w=np.ones(window_len,'d')
	else:
		w=eval('np.'+window+'(window_len)')
	y=np.convolve(w/w.sum(),s,mode='same')
	return y[window_len:-window_len+1]


# first time. after running new simulations


# for combined plot
beta = '2.000'
beta2 = '2.0'
text = 'age=0Gyr'
#text = r'$\beta=2.0$'
SAVENAME = 'test.png'

SMOOTH = True
EXPORT_FILE = False

LW = 1.5
LWD = 2.0
M_bh = 1e6*M_sun


# first, save all the dmdts



ds = [
    # run,                  chk,        label,          wlen,   log_t_cutoff
    #
    # these lines for comparison of same run, different chks
    #['m1.0_p1_b1.0',       '0040',     'chk40',        300,    0.],
    ['m1.0_p1_b1.0',        '0100',     '0100',         5,      2.],
    ['m1.0_p1_b1.0',        '0040',     '0040',         5,      2.],
    #
    # these lines for same age, different beta
    # p1
    ['m1.0_p1_b1.0',        '0100',     r'$\beta=1.0$', 300,    0.],
    ['m1.0_p1_b2.0',        '0080',     r'$\beta=2.0$', 250,    -0.1],
    ['m1.0_p1_b3.0',        '0060',     r'$\beta=3.0$', 250,    0.3],
    # p10
    ['m1.0_p10_b1.0',       '0100',     r'$\beta=1.0$', 400,    -0.4],
    ['m1.0_p10_b1.5',       '0090',     r'$\beta=1.5$', 200,    -0.15],
    ['m1.0_p10_b2.0',       '0075',     r'$\beta=2.0$', 200,    0.25],
    ['m1.0_p10_b3.0',       '0060',     r'$\beta=3.0$', 200,    0.4],
    ['m1.0_p10_b4.0',       '0050',     r'$\beta=4.0$', 200,    0.25],
    # p16
    ['m1.0_p16_b2.0',       '0075',     r'$\beta=2.0$', 300,    3],
    ['m1.0_p16_b3.0',       '0060',     r'$\beta=3.0$', 500,    0],
    ['m1.0_p16_b4.0',       '0050',     r'$\beta=4.0$', 300,    0.1],
    ['m1.0_p16_b5.0',       '0050',     r'$\beta=5.0$', 300,    0.1],
    # these lines are for same beta, different age
    #
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


fig, ax = plt.subplots()
fig2, ax2 = plt.subplots()

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


for d in ds:
    e, dm = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/'+d[0]+
        '/b10000_ev_bhbound_histogram_multitidal_hdf5_chk_'+d[1]+'.dat',
        skiprows=4)
    de = e[1]-e[0]
    dm_de = dm/de

    # calculate mdot
    e_bound = e[np.where(e<0.)]
    t = 2.*np.pi*G*M_bh/((2*np.abs(e_bound))**(3./2))
    de_dt = (1./3)*((2*np.pi*G*M_bh)**(2./3))*t**(-5./3)
    dm_de_bound = dm_de[np.where(e<0)]
    mdot = dm_de_bound*de_dt

    # take logs. better to smooth what I'm actually plotting (log)
    log_t = np.log10(t)
    log_t_yr = np.log10(t/yr)
    log_mdot = np.log10(mdot)
    log_mdot_moyr = np.log10(mdot*yr/M_sun)
    log_dm_de = np.log10(dm_de)

    # mdot eddington
    eta = 0.1
    mdot_edd = (2.2e-8)*(eta/0.1)*(M_bh/M_sun)*M_sun/yr
    log_mdot_edd = np.log10(mdot_edd)

    if SMOOTH:
        # smooth hanning dmde
        slog_dm_de = smooth(log_dm_de, window_len=d[3])
        slog_mdot_moyr = smooth(log_mdot_moyr, window_len=d[3])
        ax.plot(e/1e17, slog_dm_de, lw=LW, label=d[2])
        sel = np.where(log_t_yr<d[4])[0]
        ax2.plot(log_t_yr[sel], slog_mdot_moyr[sel], lw=LW, label=d[2])

        if EXPORT_FILE:
            ascii.write([log_t_yr[sel], slog_mdot_moyr[sel]],
                '/pfs/lawsmith/FLASH4.3/runs/results/'+d[0].replace(".","_")+
                '_'+d[1]+'.dat',
                names=['log_t_yr','log_mdot_moyr'])

    else:
        ax.plot(e/1e17, log_dm_de, lw=LW, rasterized=True, label=d[2])
        ax2.plot(log_t_yr, log_mdot_moyr, lw=LW, rasterized=True, label=d[2])

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

ax.set_xlim(-10, 10)
#ax.set_ylim(9, 15)
ax.set_xlabel(r'$E\ \mathrm{[10^{17}\ erg\ g^{-1}]}$')
ax.set_ylabel(r'$\log\ dM/dE\ \mathrm{[g^2\ erg^{-1}]}$')
ax.legend()
fig.tight_layout()
#fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/dmde_'+d[0].replace(".","_")+
#    '_'+d[1]+'.png')

#ax2.set_xlim(-2, 0)
#ax2.set_xlim(-2, 0.25)
#ax2.set_ylim(-2, 1)
ax2.set_ylabel(r'$\log\ \dot M\ {\rm [M_\odot/yr]}$')
ax2.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
ax2.legend(loc=1)
#ax2.text(-1.95,0.8,text)
#ax2.text(-1.45,-0.25,text)
fig2.tight_layout()
#fig2.savefig('/pfs/lawsmith/FLASH4.3/runs/results/mdot_'+d[0].replace(".","_")+
#    '_'+d[1]+'.png')
fig2.savefig('/pfs/lawsmith/FLASH4.3/runs/results/'+SAVENAME)
