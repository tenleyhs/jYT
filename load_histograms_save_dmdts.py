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
from scipy.ndimage.filters import gaussian_filter
import os
execfile('/pfs/lawsmith/jYT/my_settings.py')
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['font.size'] = 18

PLOT = True
LW = 1.5
M_bh = 1e6*M_sun

LOOP_THRU_SIGMAS = False
sigmas = [1,5,10,15,20,25,30,35,40,45,50]

ds = [
    # run,                  chk,        sigma,   tmin,	tmax
    # p1
    ['m1.0_p1_b1.0',        '0040',     15,      -1.4,	-0.1],
	['m1.0_p1_b1.5',		'0040',		10,		-1.65,	0.1],
	['m1.0_p1_b1.75',		'0040',		15,		-1.75,	-0.35],
	['m1.0_p1_b2.0',		'0040',		25,		-1.8,	-0.65],
	['m1.0_p1_b3.0',		'0040',		15,		-1.9,	-0.3],
    # p10
	['m1.0_p10_b1.0',		'0040',     20,      -1.3,	-0.3],
	['m1.0_p10_b1.0_256',	'0040',     50,      -1.3,	0.1],
	['m1.0_p10_b1.5',		'0040',     20,      -1.65,	-0.25],
	['m1.0_p10_b2.0',		'0040',     25,      -1.8,	-0.1],
	['m1.0_p10_b2.0_256',	'0040',     50,      -1.8,	3],
	['m1.0_p10_b2.5',		'0040',     20,      -1.9,	-0.3],
	['m1.0_p10_b3.0',		'0040',     20,      -2,	-0.3],
	['m1.0_p10_b4.0',		'0040',     40,      -2,	2],
	['m1.0_p10_b5.0',		'0040',     40,      -2,	2],
    # p16
	['m1.0_p16_b1.0',		'0040',     50,      -1.1,	2],
	['m1.0_p16_b1.5',		'0040',     50,      -1.5,	2],
	['m1.0_p16_b2.0',		'0040',     30,      -1.65,	2],
	['m1.0_p16_b3.0',		'0040',     50,      -1.85,	2],
	['m1.0_p16_b4.0',		'0040',     30,      -2,	2],
	['m1.0_p16_b5.0',		'0040',     50,      -2,	2],
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

if LOOP_THRU_SIGMAS == True:
	for d in ds:
	    for el in els:
			for sigma in sigmas:
				#fig, ax = plt.subplots()
				fig2, ax2 = plt.subplots()
				e, dm = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/'+d[0]+
				    '/b10000_'+el+'_bhbound_histogram_multitidal_hdf5_chk_'+d[1]+'.dat',
				    skiprows=4)
				de = e[1]-e[0]
				dm_de = dm/de
				log_dm_de = np.log10(dm_de)
				# smooth dmde first
				slog_dm_de = gaussian_filter(log_dm_de, sigma, mode='wrap')
				e_bound = e[np.where(e<0.)]
				t = 2.*np.pi*G*M_bh/((2*np.abs(e_bound))**(3./2))
				de_dt = (1./3)*((2*np.pi*G*M_bh)**(2./3))*t**(-5./3)
				dm_de_bound = 10**slog_dm_de[np.where(e<0)]
				dm_de_bound_orig = dm_de[np.where(e<0)]
				mdot = dm_de_bound*de_dt
				mdot_orig = dm_de_bound_orig*de_dt

				# take logs. better to smooth what I'm actually plotting (log)
				log_t = np.log10(t)
				log_t_yr = np.log10(t/yr)
				log_mdot = np.log10(mdot)
				log_mdot_moyr = np.log10(mdot*yr/M_sun)
				log_mdot_moyr_orig = np.log10(mdot_orig*yr/M_sun)

				# plot unsmoothed
				#ax.scatter(e/1e17, log_dm_de, c='C0', rasterized=True,
				#	alpha=0.5, s=1, edgecolors='none')
				ax2.scatter(log_t_yr, log_mdot_moyr_orig, c='C0', rasterized=True,
					alpha=0.5, s=1, edgecolors='none')

				#ax.plot(e/1e17, slog_dm_de, c='C1', lw=LW)
				sel = np.where((log_t_yr>d[3]) & (log_t_yr<d[4]))[0]
				x = log_t_yr[sel]
				y = log_mdot_moyr[sel]
				ax2.plot(x, y, c='C1', lw=LW)

				# extend dmdt from slope near end
				# could maybe have a better slope function
				if x[-1] < 2.0:
					ninf = (y[-20] - y[-1]) / (x[-20] - x[-1])
					ax2.plot([x[-1], 2.0], [y[-1], y[-1] + (2.0 - x[-1])*ninf], c='C2', lw=LW)

				#ax.set_xlim(-10, 10)
				#ax.set_ylim(9, 15)
				#ax.set_xlabel(r'$E\ \mathrm{[10^{17}\ erg\ g^{-1}]}$')
				#ax.set_ylabel(r'$\log\ dM/dE\ \mathrm{[g^2\ erg^{-1}]}$')
				#fig.tight_layout()
				#directory = '/pfs/lawsmith/FLASH4.3/runs/results/dmdes/gaussian_wrap/' + d[0].replace(".","_")
				#if not os.path.exists(directory):
				#	os.makedirs(directory)
				#fig.savefig(directory + '/dmde_'
				#	+ d[0].replace(".","_") + '_' + d[1] + '_' + el + '_'
				#	+ str(sigma) + '.pdf')

				ax2.set_xlim(-2, 2)
				ax2.set_ylabel(r'$\log\ \dot M\ {\rm [M_\odot/yr]}$')
				ax2.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
				fig2.tight_layout()
				directory = '/pfs/lawsmith/FLASH4.3/runs/results/dmdts/dmde_first/' + d[0].replace(".","_")
				if not os.path.exists(directory):
					os.makedirs(directory)
				fig2.savefig(directory + '/dmdt_'
					+ d[0].replace(".","_") + '_' + d[1] + '_' + el + '_'
					+ str(sigma) + '.pdf')

				plt.close('all')
else:
	for d in ds:
	    for el in els:
			#fig, ax = plt.subplots()
			fig2, ax2 = plt.subplots()
			e, dm = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/'+d[0]+
			    '/b10000_'+el+'_bhbound_histogram_multitidal_hdf5_chk_'+d[1]+'.dat',
			    skiprows=4)
			de = e[1]-e[0]
			dm_de = dm/de
			log_dm_de = np.log10(dm_de)
			# smooth dmde first
			slog_dm_de = gaussian_filter(log_dm_de, d[2], mode='wrap')
			e_bound = e[np.where(e<0.)]
			t = 2.*np.pi*G*M_bh/((2*np.abs(e_bound))**(3./2))
			de_dt = (1./3)*((2*np.pi*G*M_bh)**(2./3))*t**(-5./3)
			dm_de_bound = 10**slog_dm_de[np.where(e<0)]
			dm_de_bound_orig = dm_de[np.where(e<0)]
			mdot = dm_de_bound*de_dt
			mdot_orig = dm_de_bound_orig*de_dt

			# take logs. better to smooth what I'm actually plotting (log)
			log_t = np.log10(t)
			log_t_yr = np.log10(t/yr)
			log_mdot = np.log10(mdot)
			log_mdot_moyr = np.log10(mdot*yr/M_sun)
			log_mdot_moyr_orig = np.log10(mdot_orig*yr/M_sun)

			# plot unsmoothed
			#ax.scatter(e/1e17, log_dm_de, c='C0', rasterized=True,
			#	alpha=0.5, s=1, edgecolors='none')
			ax2.scatter(log_t_yr, log_mdot_moyr_orig, c='C0', rasterized=True,
				alpha=0.5, s=1, edgecolors='none')

			#ax.plot(e/1e17, slog_dm_de, c='C1', lw=LW)
			sel = np.where((log_t_yr>d[3]) & (log_t_yr<d[4]))[0]
			x = log_t_yr[sel]
			y = log_mdot_moyr[sel]
			ax2.plot(x, y, c='C1', lw=LW)

			# extend dmdt from slope near end
			# could maybe have a better slope function
			if x[-1] < 2.0:
				ninf = (y[-20] - y[-1]) / (x[-20] - x[-1])
				ax2.plot([x[-1], 2.0], [y[-1], y[-1] + (2.0 - x[-1])*ninf], c='C2', lw=LW)
				x = np.append(x, 2.0)
				y = np.append(y, (2.0 - x[-1])*ninf)

			ascii.write([x, y],
			    '/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/'+d[0].replace(".","_")+
			    '_'+d[1]+'_'+el+'.dat', overwrite=True,
			    names=['log_t_yr','log_mdot_moyr'])

			#ax.set_xlim(-10, 10)
			#ax.set_ylim(9, 15)
			#ax.set_xlabel(r'$E\ \mathrm{[10^{17}\ erg\ g^{-1}]}$')
			#ax.set_ylabel(r'$\log\ dM/dE\ \mathrm{[g^2\ erg^{-1}]}$')
			#fig.tight_layout()
			#directory = '/pfs/lawsmith/FLASH4.3/runs/results/dmdes/gaussian_wrap/' + d[0].replace(".","_")
			#if not os.path.exists(directory):
			#	os.makedirs(directory)
			#fig.savefig(directory + '/dmde_'
			#	+ d[0].replace(".","_") + '_' + d[1] + '_' + el + '_'
			#	+ str(sigma) + '.pdf')

			ax2.set_xlim(-2, 2)
			ax2.set_ylabel(r'$\log\ \dot M\ {\rm [M_\odot/yr]}$')
			ax2.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
			fig2.tight_layout()
			directory = '/pfs/lawsmith/FLASH4.3/runs/results/dmdts/final'
			fig2.savefig(directory + '/dmdt_'
				+ d[0].replace(".","_") + '_' + d[1] + '_' + el + '.pdf')

			plt.close('all')
