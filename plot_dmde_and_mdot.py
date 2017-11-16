import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import UnivariateSpline
#execfile('/Users/lawsmith/Dropbox/helpers/my_settings.py')
execfile('/nobackup/jlawsmit/jYT/my_settings.py')

plt.rcParams['legend.fontsize'] = 16
plt.rcParams['font.size'] = 18

# PARAMS
M_bh = 1e6*M_sun
lw = 1.5

# dir
d = '/nobackup/jlawsmit/m1.0_p10_b2.0/'

beta = '2.000'

title = 'm=1.0, p=10, b=2.0'

bins = '50000'

hists = [
	'b50000_ev_bhbound_histogram_multitidal_hdf5_chk_0040.dat',
	#'h1_bhbound_histogram_multitidal_hdf5_chk_0030.dat'
	]

labels = [
	'total',
	#'h1'
]

colors = [
	'black',
	#red
]

# smoothing
do_smoothing = False
wl = 10

# limits by eye for e in dmde plot
elim = 2e17

# limits by eye for t in dmdt plot. from jupyter
min_log_t = 6.5
max_log_t = 9.0

# load Guillochon 2013 dmdts
#g13_43 = np.loadtxt('/Users/lawsmith/Dropbox/temp/dmdts/4-3/' + beta + '.dat')
#g13_53 = np.loadtxt('/Users/lawsmith/Dropbox/temp/dmdts/5-3/' + beta + '.dat')
g13_43 = np.loadtxt('/nobackup/jlawsmit/dmdts/4-3/' + beta + '.dat')
g13_53 = np.loadtxt('/nobackup/jlawsmit/dmdts/5-3/' + beta + '.dat')



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
fig2, ax2 = plt.subplots()

for h, l, c in zip(hists, labels, colors):

	e, dm = np.loadtxt(d+h, skiprows=4)
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
	log_dm_de = np.log10(dm_de)

	# mdot eddington
	eta = 0.1
	mdot_edd = (2.2e-8)*(eta/0.1)*(M_bh/M_sun)*M_sun/yr
	log_mdot_edd = np.log10(mdot_edd)

	if do_smoothing:
		# smooth hanning dmde
		slog_dm_de = smooth(log_dm_de, window_len=wl)
		slog_mdot_moyr = smooth(log_mdot_moyr, window_len=wl)
		ax.plot(e/1e17, slog_dm_de, color=c, lw=lw)
		ax2.plot(log_t_yr, slog_mdot_moyr, color=c, lw=lw, label=l)

	else:
		ax.plot(e/1e17, log_dm_de, color=c, lw=lw, rasterized=True)
		ax2.plot(log_t_yr, log_mdot_moyr, color=c, lw=lw, rasterized=True, label=l)


ax.set_xlim(-10, 10)
#ax.set_ylim(9, 15)
ax.set_xlabel(r'$E\ \mathrm{[10^{17}\ erg\ g^{-1}]}$')
ax.set_ylabel(r'$\log\ dM/dE\ \mathrm{[g^2\ erg^{-1}]}$')
fig.tight_layout()
if do_smoothing:
	fig.savefig(d+'dmde_' + bins + '_sm.png')
else:
	fig.savefig(d+'dmde_' + bins + '.png')


# plot Guillochon 2013
ax2.plot(g13_43[:,0], g13_43[:,1], ls=':', lw=2, color=orange, label='4/3')
ax2.plot(g13_53[:,0], g13_53[:,1], ls=':', lw=2, color=blue, label='5/3')

#ax2.axhline(0., c='k', ls='--')
ax2.set_xlim(-3, 3)
ax2.set_ylim(-6, 2)
ax2.set_ylabel(r'$\log\ \dot M\ {\rm [M_\odot/yr]}$')
ax2.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
ax2.legend()
ax2.set_title(title)
fig2.tight_layout()
if do_smoothing:
	fig2.savefig(d+'mdot_' + bins + '_sm.png')
else:
	fig2.savefig(d+'mdot_' + bins + '.png')
