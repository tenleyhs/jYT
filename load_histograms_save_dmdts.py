"""
read in all the histograms and save dmdts
optionally make plots. plot raw with alpha to compare. go through these plots to
decide what window lens, tmin, tmax to choose
comment out some ds if just ran a few new simulations and don't want to do everything
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.io import ascii
from scipy.interpolate import UnivariateSpline
#from scipy.ndimage.filters import gaussian_filter
execfile('/pfs/lawsmith/jYT/my_settings.py')
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['font.size'] = 18
import os

"""
TODO might want to use the gaussian filter smoothing than brenna emailed me
instead.
https://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.ndimage.filters.gaussian_filter.html
"""

def smooth(x,window_len=11,window='hanning'):
	if x.ndim != 1:
		raise ValueError, "smooth only accepts 1 dimension arrays."
	if x.size < window_len:
		raise ValueError, "Input vector needs to be bigger than window size."
	if window_len<3:
		return x
	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError, "Window is not one of the above"
	s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
	if window == 'flat': #moving average
		w=np.ones(window_len,'d')
	else:
		w=eval('np.'+window+'(window_len)')
	y=np.convolve(w/w.sum(),s,mode='same')
	return y[window_len:-window_len+1]

PLOT = True
LW = 1.5
M_bh = 1e6*M_sun


"""
TODO a better way to decide on smoothing kernel to use might be to plot a bunch
at once, and then I just have to decide once, rather than iterate
"""

#wlens = [10,20,30,40,50,60,70,80,90,100] # for hanning
wlens = [100,200,300,400,500,600,700] # for UnivariateSpline

ds = [
    # run,                  chk,        wlen,   tmin,	tmax
    # p1
    ['m1.0_p1_b1.0',        '0040',     40,      -1.4,	0.2],
	['m1.0_p1_b1.5',		'0040',		40,		-1.7,	0.4],
	['m1.0_p1_b1.75',		'0040',		50,		-1.8,	-0.2],
	['m1.0_p1_b2.0',		'0040',		50,		-1.85,	-0.3],
	['m1.0_p1_b3.0',		'0040',		50,		-1.9,	-0.1],
    # p10
	['m1.0_p10_b1.0',		'0040',     50,      -1.3,	-0.2],
	['m1.0_p10_b1.0_256',	'0040',     70,      -1.3,	0.1],
	['m1.0_p10_b1.5',		'0040',     60,      -1.65,	-0.2],
	['m1.0_p10_b2.0',		'0040',     60,      -2.0,	0.1],
	['m1.0_p10_b2.0_256',	'0040',     50,      -2,	3],
	['m1.0_p10_b2.5',		'0040',     60,      -1.9,	-0.1],
	['m1.0_p10_b3.0',		'0040',     60,      -2,	-0.1],
	['m1.0_p10_b4.0',		'0040',     40,      -2,	2],
	['m1.0_p10_b5.0',		'0040',     40,      -2,	2],
    # p16
	['m1.0_p16_b1.0',		'0040',     50,      -1.25,	2],
	['m1.0_p16_b1.5',		'0040',     40,      -1.5,	2],
	['m1.0_p16_b2.0',		'0040',     40,      -1.9,	2],
	['m1.0_p16_b3.0',		'0040',     50,      -2,	0.25],
	['m1.0_p16_b4.0',		'0040',     50,      -2,	1],
	['m1.0_p16_b5.0',		'0040',     50,      -2,	1.5],
]

    #['m1.0_p1_b1.0',        '0100',     5,      2.],
    #['m1.0_p1_b2.0',        '0080',     250,    -0.1],
    #['m1.0_p1_b3.0',        '0060',     250,    0.3],
    #['m1.0_p10_b1.0',       '0100',     400,    -0.4],
    #['m1.0_p10_b1.0_256',   '0100',     300,    -0.4],
    #['m1.0_p10_b1.5',       '0090',     200,    -0.15],
    #['m1.0_p10_b2.0',       '0075',     200,    0.25],
    #['m1.0_p10_b3.0',       '0060',     200,    0.4],
    #['m1.0_p10_b4.0',       '0050',     200,    0.25],
	#['m1.0_p16_b1.0_pfe',   '0052',     500,    1],
    #['m1.0_p16_b2.0',       '0075',     300,    3],
    #['m1.0_p16_b3.0',       '0060',     500,    0],
    #['m1.0_p16_b4.0',       '0050',     300,    0.1],
    #['m1.0_p16_b5.0',       '0050',     300,    0.1],


#els = ['ev', 'h1', 'he4', 'o16', 'c12', 'ne20', 'n14']
els = ['ev']

for d in ds:
    for el in els:
		for wlen in wlens:
			#fig, ax = plt.subplots()
			fig2, ax2 = plt.subplots()
			e, dm = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/'+d[0]+
			    '/b10000_'+el+'_bhbound_histogram_multitidal_hdf5_chk_'+d[1]+'.dat',
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
			#log_dm_de = np.log10(dm_de)

			# mdot eddington
			#eta = 0.1
			#mdot_edd = (2.2e-8)*(eta/0.1)*(M_bh/M_sun)*M_sun/yr
			#log_mdot_edd = np.log10(mdot_edd)

			# plot unsmoothed
			#ax.plot(e/1e17, log_dm_de, lw=LW, rasterized=True, alpha=0.4)
			#ax2.plot(log_t_yr, log_mdot_moyr, lw=LW, rasterized=True, alpha=0.4)
			ax2.scatter(log_t_yr, log_mdot_moyr, c='C0', lw=LW, rasterized=True,
				alpha=0.5, s=1, edgecolors='none')

			# smooth. TODO choose hanning, UnivariateSpline, gaussian_filter
			# TODO choose if dmde first
			#slog_dm_de = smooth(log_dm_de, window_len=d[2])
			#slog_mdot_moyr = smooth(log_mdot_moyr, window_len=wlen)
			#sel = np.where((log_t_yr>d[3]) & (log_t_yr<d[4]))[0]
			#x = log_t_yr[sel]
			#y = slog_mdot_moyr[sel]
			#ax.plot(e/1e17, slog_dm_de, lw=LW)
			#ax2.plot(x, y, c='C1', lw=LW)


			# spline version
			# TODO maybe choose k
			w = np.isinf(log_mdot_moyr)
			log_mdot_moyr[w] = 0.
			spl = UnivariateSpline(log_t_yr, log_mdot_moyr, w=~w, k=3, s=wlen)
			xs = np.linspace(d[3], d[4], 100)
			ax2.plot(xs, spl(xs), c='C1', lw=LW)


			# gaussian_filter version


			# extend dmdt from slope near end
			# TODO make sure I'm doing this robustly. better slope function
			#ninf = (y[-2] - y[-1]) / (x[-2] - x[-1])

			# TODO temp commented out
			#ascii.write([log_t_yr[sel], slog_mdot_moyr[sel]],
			#    '/pfs/lawsmith/FLASH4.3/runs/results/dmdts/'+d[0].replace(".","_")+
			#    '_'+d[1]+'_'+el+'.dat', overwrite=True,
			#    names=['log_t_yr','log_mdot_moyr'])

			#ax.set_xlim(-10, 10)
			#ax.set_ylim(9, 15)
			#ax.set_xlabel(r'$E\ \mathrm{[10^{17}\ erg\ g^{-1}]}$')
			#ax.set_ylabel(r'$\log\ dM/dE\ \mathrm{[g^2\ erg^{-1}]}$')
			#fig.tight_layout()
			#fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/dmde_'+d[0].replace(".","_")+
			#    '_'+d[1]+'_'+el+'.png')

			#ax2.set_xlim(-2, 0)
			#ax2.set_xlim(-2, 0.25)
			#ax2.set_ylim(-2, 1)
			ax2.set_ylabel(r'$\log\ \dot M\ {\rm [M_\odot/yr]}$')
			ax2.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
			fig2.tight_layout()
			#directory = '/pfs/lawsmith/FLASH4.3/runs/results/dmdts/mdot_hanning/' + d[0].replace(".","_")
			directory = '/pfs/lawsmith/FLASH4.3/runs/results/dmdts/mdot_univariatespline/' + d[0].replace(".","_")
			if not os.path.exists(directory):
				os.makedirs(directory)
			fig2.savefig(directory + '/dmdt_'
				+ d[0].replace(".","_") + '_' + d[1] + '_' + el + '_'
				+ str(wlen) + '.pdf')

			plt.close('all')
