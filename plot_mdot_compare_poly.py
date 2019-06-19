import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import UnivariateSpline
#execfile('/nobackup/jlawsmit/jYT/my_settings.py')
execfile('/pfs/lawsmith/jYT/my_settings.py')
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['font.size'] = 18

# PARAMS
M_bh = 1e6*M_sun
lw = 1.5
lwd = 2.0
beta = '2.000'
beta2 = '2.0'
fname = 'p11016_b20_compare_poly'
text = r'$\beta=2.0$'
do_smoothing = True

ds = [
	#['m1.0_p1_b1.0/b10000_ev_bhbound_histogram_multitidal_hdf5_chk_0100.dat','age=0Gyr',300,0.],
	#['m1.0_p10_b1.0/b10000_ev_bhbound_histogram_multitidal_hdf5_chk_0100.dat','age=4.8Gyr',400,-0.4],
	['m1.0_p1_b2.0/b10000_ev_bhbound_histogram_multitidal_hdf5_chk_0080.dat','age=0Gyr',250,-0.1],
	['m1.0_p10_b2.0_128/b10000_ev_bhbound_histogram_multitidal_hdf5_chk_0075.dat','age=4.8Gyr',200,0.25],
	['m1.0_p16_b2.0/b10000_ev_bhbound_histogram_multitidal_hdf5_chk_0075.dat','age=8.4Gyr',300,3],
	#['m1.0_p1_b3.0/b10000_ev_bhbound_histogram_multitidal_hdf5_chk_0060.dat','age=0Gyr',250,0.3],
	#['m1.0_p10_b3.0/b10000_ev_bhbound_histogram_multitidal_hdf5_chk_0060.dat','age=4.8Gyr',200,0.4],
	#['m1.0_p16_b3.0/b10000_ev_bhbound_histogram_multitidal_hdf5_chk_0060.dat','age=8.4Gyr',500,0],
	]

# for dmde's
def smooth(x,window_len=11,window='hanning'):
	if x.ndim != 1:
		raise ValueError, "smooth only accepts 1 dimension arrays."
	if x.size < window_len:
		raise ValueError, "Input vector needs to be bigger than window size."
	if window_len<3:
		return x
	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
	s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
	if window == 'flat': #moving average
		w=np.ones(window_len,'d')
	else:
		w=eval('np.'+window+'(window_len)')
	y=np.convolve(w/w.sum(),s,mode='same')
	return y[window_len:-window_len+1]


fig, ax = plt.subplots()

# Guillochon 2013 dmdts
g13_43 = np.loadtxt('/pfs/lawsmith/dmdts/4-3/' + beta + '.dat')
g13_53 = np.loadtxt('/pfs/lawsmith/dmdts/5-3/' + beta + '.dat')

ax.plot(g13_43[:,0], g13_43[:,1], lw=lwd, label='GRR13 4/3', ls=':', c='k')
#ax.plot(g13_53[:,0], g13_53[:,1], lw=lw, label='GRR13 5/3', ls='-.', c='k')

# Lodato/Monica dmdt
"""
m17_53 = np.loadtxt('/pfs/lawsmith/dmdts/poly53-b' + beta2 + 'dmdT.data')
ax.plot(np.log10(m17_53[:,1]/yr), np.log10(m17_53[:,0]*yr/M_sun), lw=lw, label='analytic 5/3', color='gray')
"""

if beta2 == '1.0':
    m17_43 = np.loadtxt('/pfs/lawsmith/dmdts/poly_dmdt.data', skiprows=1)
    ax.plot(np.log10(m17_43[:,5]/yr), np.log10(m17_43[:,4]*yr/M_sun), lw=lwd, label='analytic 4/3', color='gray', ls=':')

for d in ds:
    e, dm = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/'+d[0], skiprows=4)
    de = e[1]-e[0]
    dm_de = dm/de

    # calculate mdot
    e_bound = e[np.where(e<0.)]
    t = 2.*np.pi*G*M_bh/((2*np.abs(e_bound))**(3./2))
    de_dt = (1./3)*((2*np.pi*G*M_bh)**(2./3))*t**(-5./3)
    dm_de_bound = dm_de[np.where(e<0)]
    mdot = dm_de_bound*de_dt

    # take logs. better to smooth what I'm actually plotting (log-linear or log-log)
    log_t = np.log10(t)
    log_t_yr = np.log10(t/yr)
    log_mdot = np.log10(mdot)
    log_mdot_moyr = np.log10(mdot*yr/M_sun)

    if do_smoothing:
        # smooth hanning dmde
        slog_mdot_moyr = smooth(log_mdot_moyr, window_len=d[2])
        sel = np.where(log_t_yr<d[3])[0]
        ax.plot(log_t_yr[sel], slog_mdot_moyr[sel], lw=2.0, label=d[1])

    else:
        ax.plot(log_t_yr, log_mdot_moyr, lw=lw, rasterized=True, label=d[1])


"""
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
        ax.plot(g13_43[:,0], g13_43[:,1], ls='-', lw=1.5, color='orange', label='4/3')
    else:
        ax.plot(g13_43[:,0], g13_43[:,1], ls='-', lw=0.5, color='orange')

    ax.plot(g13_53[:,0], g13_53[:,1], ls=':', lw=2, color=blue, label='5/3')
"""

ax.set_xlim(-2, 0)
ax.set_ylim(-2, 1)
ax.set_ylabel(r'$\log\ \dot M\ {\rm [M_\odot/yr]}$')
ax.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
ax.legend(loc=1)
ax.text(-1.9,0.8,text)
fig.tight_layout(pad=0.3)
fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/compare_poly/mdot_'+fname+'.png')
