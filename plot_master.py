"""
trying to combine all files with plot_ having to read in histograms
maybe a faster/better way would be to load all the files first, into some
dictionary or something. because right now repeats.
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

# option for seaborn style background
#mpl.style.use('seaborn')

LW = 1.5
LWD = 2.0
S = 30
S2 = 200

# the complete list, for reference
"""
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
"""

"""
### for total mdot, not split into elements
def myplot_ev(ds=None, text=None, savename=None,
				x2=None, y2=None, l2=None,
				x3=None, y3=None, l3=None):
	el = 'ev'
	fig, ax = plt.subplots()
	for d in ds:
		# log_t_yr, log_mdot_moyr
		x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/'
				+ d[0].replace(".","_") + '_' + d[1] + '_' + el + '.dat',
				skiprows=1, unpack=True)
		ax.plot(x, y, lw=LW, label=d[2])

	if x2 is not None:
		ax.plot(x2, y2, lw=LWD, label=l2, ls=':', c='k')
	if x3 is not None:
		ax.plot(x3, y3, lw=LWD, label=l3, ls=':', c='gray')

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
	#['m1.0_p10_b1.0_256',	'0040',		r'$\beta=1.0 (256)$'],
	['m1.0_p10_b1.5',		'0040',		r'$\beta=1.5$'],
	['m1.0_p10_b2.0',		'0040',		r'$\beta=2.0$'],
	#['m1.0_p10_b2.0_256',	'0040',		r'$\beta=2.0 (256)$'],
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
# b1.0
# Guillochon 2013 dmdts
beta = '1.000'
g13_43 = np.loadtxt('/pfs/lawsmith/dmdts/4-3/' + beta + '.dat')
#g13_53 = np.loadtxt('/pfs/lawsmith/dmdts/5-3/' + beta + '.dat')
# Lodato/Monica dmdts
# only have 4/3 for beta = 1.0. can ask Monica or do myself for other betas
m18_43 = np.loadtxt('/pfs/lawsmith/dmdts/poly_dmdt.data', skiprows=1)
#m17_53 = np.loadtxt('/pfs/lawsmith/dmdts/poly53-b' + beta2 + 'dmdT.data')
ds = [
    ['m1.0_p1_b1.0',        '0040',     'age=0Gyr'],
    ['m1.0_p10_b1.0',       '0040',     'age=4.8Gyr'],
	['m1.0_p16_b1.0',       '0040',     'age=8.4Gyr'],
]
text = r'$\beta=1.0$'
savename = 'mdot_b1_0_allages.pdf'
myplot_ev(ds=ds, text=text, savename=savename,
	x2=g13_43[:,0], y2=g13_43[:,1], l2='GRR13 4/3',
	x3=np.log10(m18_43[:,5]/yr), y3=np.log10(m18_43[:,4]*yr/M_sun), l3='analytic 4/3')
#ax.plot(g13_53[:,0], g13_53[:,1], lw=lw, label='GRR13 5/3', ls='-.', c='k')
#ax2.plot(np.log10(m18_53[:,1]/yr), np.log10(m18_53[:,0]*yr/M_sun), lw=lw, label='analytic 5/3', color='gray')

# b2.0
beta = '2.000'
g13_43 = np.loadtxt('/pfs/lawsmith/dmdts/4-3/' + beta + '.dat')
ds = [
    ['m1.0_p1_b2.0',        '0040',     'age=0Gyr'],
    ['m1.0_p10_b2.0',       '0040',     'age=4.8Gyr'],
	['m1.0_p16_b2.0',       '0040',     'age=8.4Gyr'],
]
text = r'$\beta=2.0$'
savename = 'mdot_b2_0_allages.pdf'
myplot_ev(ds=ds, text=text, savename=savename,
	x2=g13_43[:,0], y2=g13_43[:,1], l2='GRR13 4/3')

# b3.0
ds = [
    ['m1.0_p1_b3.0',        '0040',     'age=0Gyr'],
    ['m1.0_p10_b3.0',       '0040',     'age=4.8Gyr'],
	['m1.0_p16_b3.0',       '0040',     'age=8.4Gyr'],
]
text = r'$\beta=3.0$'
savename = 'mdot_b3_0_allages.pdf'
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
	els = ['ev', 'h1', 'he4', 'o16', 'c12', 'ne20', 'n14']
	for d in ds:
		fig, ax = plt.subplots()
		for el in els:
			# log_t_yr, log_mdot_moyr
			x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/'
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
			    el_lr = (10**y/h1_mdot)/(f_o16_ZAMS/f_h1_ZAMS)
			elif el == 'c12':
			    el_lr = (10**y/h1_mdot)/(f_c12_ZAMS/f_h1_ZAMS)
			elif el == 'ne20':
			    el_lr = (10**y/h1_mdot)/(f_ne20_ZAMS/f_h1_ZAMS)
			elif el == 'n14':
			    el_lr = (10**y/h1_mdot)/(f_n14_ZAMS/f_h1_ZAMS)

			ax.plot(np.log10(10**x/tpeak), el_lr, lw=LW, label=el)

		# for t/tpeak
		#ax.set_xlim(-0.5, 1.5)
		#ax.set_ylim(0, 3)
		#ax.text(-0.48, 2.75, d[2])
		ax.set_xlim(-0.5, 3)
		ax.set_ylim(0, 4)
		ax.text(-0.48, 3.75, d[2])

		ax.set_xlabel(r'$\log\ t/t_{\rm peak}$')
		ax.set_ylabel(r'$X/X_\odot$')
		# TODO add t along the top in years
		# for t
	    #ax.set_xlim(-1.5, 0.0)
		#ax.text(-1.45, 2.8, d[2])
		#ax.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
		ax.axhline(1, ls=':', lw=LWD, c='k', alpha=0.5)
		ax.legend(loc=1)
		fig.tight_layout()
		fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/paper/mdot_comp_solar_' + d[3] + '.pdf')
		plt.close('all')

ds = [
    ['m1.0_p1_b1.0',        '0040',		'age=0Gyr, '+r'$\beta=1.0$',		'p1_b1_0'],
	['m1.0_p1_b1.5',		'0040',		'age=0Gyr, '+r'$\beta=1.5$',		'p1_b1_5'],
	['m1.0_p1_b1.75',		'0040',		'age=0Gyr, '+r'$\beta=1.75$',		'p1_b1_75'],
	['m1.0_p1_b2.0',		'0040',		'age=0Gyr, '+r'$\beta=2.0$',		'p1_b2_0'],
	['m1.0_p1_b3.0',		'0040',		'age=0Gyr, '+r'$\beta=3.0$',		'p1_b3_0'],
	['m1.0_p10_b1.0',		'0040',		'age=4.8Gyr, '+r'$\beta=1.0$',		'p10_b1_0'],
	['m1.0_p10_b1.0_256',	'0040',		'age=4.8Gyr, '+r'$\beta=1.0$ (256)',		'p10_b1_0_256'],
	['m1.0_p10_b1.5',		'0040',		'age=4.8Gyr, '+r'$\beta=1.5$',		'p10_b1_5'],
	['m1.0_p10_b2.0',		'0040',		'age=4.8Gyr, '+r'$\beta=2.0$',		'p10_b2_0'],
	['m1.0_p10_b2.0_256',	'0040',		'age=4.8Gyr, '+r'$\beta=2.0$ (256)',		'p10_b2_0_256'],
	['m1.0_p10_b2.5',		'0040',		'age=4.8Gyr, '+r'$\beta=2.5$',		'p10_b2_5'],
	['m1.0_p10_b3.0',		'0040',		'age=4.8Gyr, '+r'$\beta=3.0$',		'p10_b3_0'],
	['m1.0_p10_b4.0',		'0040',		'age=4.8Gyr, '+r'$\beta=4.0$',		'p10_b4_0'],
	['m1.0_p10_b5.0',		'0040',		'age=4.8Gyr, '+r'$\beta=5.0$',		'p10_b5_0'],
	['m1.0_p16_b1.0',		'0040',		'age=8.4Gyr, '+r'$\beta=1.0$',		'p16_b1_0'],
	['m1.0_p16_b1.5',		'0040',		'age=8.4Gyr, '+r'$\beta=1.5$',		'p16_b1_5'],
	['m1.0_p16_b2.0',		'0040',		'age=8.4Gyr, '+r'$\beta=2.0$',		'p16_b2_0'],
	['m1.0_p16_b3.0',		'0040',		'age=8.4Gyr, '+r'$\beta=3.0$',		'p16_b3_0'],
	['m1.0_p16_b4.0',		'0040',		'age=8.4Gyr, '+r'$\beta=4.0$',		'p16_b4_0'],
	['m1.0_p16_b5.0',		'0040',		'age=8.4Gyr, '+r'$\beta=5.0$',		'p16_b5_0'],
]
myplot_lr(ds)
"""


### for summary figure
# 1st way of doing it. my parameter grid
mpl.rcParams['figure.figsize'] = (7,5)
mpl.rcParams['ytick.right'] = False
fig, ax = plt.subplots()
dss = [
        [
    ['m1.0_p1_b1.0',        1.0],
	['m1.0_p1_b1.5',		1.5],
	['m1.0_p1_b1.75',		1.75],
	['m1.0_p1_b2.0',		2.0],
	['m1.0_p1_b3.0',		3.0]
        ],
        [
	['m1.0_p10_b1.0',		1.0],
	#['m1.0_p10_b1.0_256',	1.0],
	['m1.0_p10_b1.5',		1.5],
	['m1.0_p10_b2.0',		2.0],
	#['m1.0_p10_b2.0_256',	2.0],
	['m1.0_p10_b2.5',		2.5],
	['m1.0_p10_b3.0',		3.0],
	['m1.0_p10_b4.0',		4.0],
	#['m1.0_p10_b5.0',		5.0]
        ],
        [
	['m1.0_p16_b1.0',		1.0],
	['m1.0_p16_b1.5',		1.5],
	['m1.0_p16_b2.0',		2.0],
	['m1.0_p16_b3.0',		3.0],
	['m1.0_p16_b4.0',		4.0],
	#['m1.0_p16_b5.0',		5.0]
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
			x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/'
					+ d[0].replace(".","_") + '_0040_' + el + '.dat',
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
				beta_array.append(d[1])
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
fig.tight_layout()
fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/paper/mdot_comp_solar_summary1.pdf')
plt.close('all')
# TODO working on the above. could do a few different versions.
mpl.rcParams['figure.figsize'] = (6,5)
mpl.rcParams['ytick.right'] = True



# 2nd way of doing it. just plotting N14, a few betas, comparing ages
fig, ax = plt.subplots()
dss = [
# no p1 because straight lines
        [
	['m1.0_p10_b1.0',		1.0],
	['m1.0_p10_b2.0',		2.0],
	['m1.0_p10_b3.0',		3.0],
    ['m1.0_p10_b4.0',		4.0],
        ],
        [
	['m1.0_p16_b1.0',		1.0],
	['m1.0_p16_b2.0',		2.0],
	['m1.0_p16_b3.0',		3.0],
    ['m1.0_p16_b4.0',		4.0],
        ]
]
labels = ['age=4.8Gyr', 'age=8.4Gyr']
carr = ['C0', 'C1']
lss = ['-', '--', ':', '-.']
lws = [1.5, 1.5, 2.5, 2]
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
			x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/'
					+ d[0].replace(".","_") + '_0040_' + el + '.dat',
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

ax.plot(-9,-9, ls=lss[0], c='k', label=r'$\beta=1.0$')
ax.plot(-9,-9, ls=lss[1], c='k', label=r'$\beta=2.0$')
ax.plot(-9,-9, ls=lss[2], c='k', label=r'$\beta=3.0$')
ax.plot(-9,-9, ls=lss[3], c='k', label=r'$\beta=4.0$')
ax.set_xlim(-0.5, 3)
ax.set_ylim(0, 4)
ax.set_xlabel(r'$\log\ t/t_{\rm peak}$')
ax.set_ylabel(r'$\left(X/X_\odot\right)_{^{14}{\rm N}}$')
ax.legend(loc=2)
fig.tight_layout()
fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/paper/mdot_comp_solar_summary2.pdf')
plt.close('all')




# 3rd way of doing it. like Monica's figure
fig, ax = plt.subplots()
dss = [
        [
    ['m1.0_p1_b1.0',        1.0],
	['m1.0_p1_b1.5',		1.5],
	['m1.0_p1_b1.75',		1.75],
	['m1.0_p1_b2.0',		2.0],
	['m1.0_p1_b3.0',		3.0]
        ],
        [
	['m1.0_p10_b1.0',		1.0],
	#['m1.0_p10_b1.0_256',	1.0],
	['m1.0_p10_b1.5',		1.5],
	['m1.0_p10_b2.0',		2.0],
	#['m1.0_p10_b2.0_256',	2.0],
	['m1.0_p10_b2.5',		2.5],
	['m1.0_p10_b3.0',		3.0],
	['m1.0_p10_b4.0',		4.0],
	#['m1.0_p10_b5.0',		5.0]
        ],
        [
	['m1.0_p16_b1.0',		1.0],
	['m1.0_p16_b1.5',		1.5],
	['m1.0_p16_b2.0',		2.0],
	['m1.0_p16_b3.0',		3.0],
	['m1.0_p16_b4.0',		4.0],
	#['m1.0_p16_b5.0',		5.0]
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
	for d in ds:
		for el in els:
			# log_t_yr, log_mdot_moyr
			x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/'
					+ d[0].replace(".","_") + '_0040_' + el + '.dat',
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
				g_array.append(el_lr[(np.abs(np.log10(10**x/tpeak) - 1.0)).argmin()])

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
ax.set_ylabel(r'$\log\ \dot M_{\rm peak}\ {\rm [M_\sun/yr]}$')
fig.tight_layout()
fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/paper/mdot_comp_solar_summary3.pdf')
plt.close('all')





"""
### for power law index
dss = [
        [
    ['m1.0_p1_b1.0',        1.0],
	['m1.0_p1_b1.5',		1.5],
	['m1.0_p1_b1.75',		1.75],
	['m1.0_p1_b2.0',		2.0],
	['m1.0_p1_b3.0',		3.0]
        ],
        [
	['m1.0_p10_b1.0',		1.0],
	#['m1.0_p10_b1.0_256',	1.0],
	['m1.0_p10_b1.5',		1.5],
	['m1.0_p10_b2.0',		2.0],
	#['m1.0_p10_b2.0_256',	2.0],
	['m1.0_p10_b2.5',		2.5],
	['m1.0_p10_b3.0',		3.0],
	['m1.0_p10_b4.0',		4.0],
	#['m1.0_p10_b5.0',		5.0]
        ],
        [
	['m1.0_p16_b1.0',		1.0],
	['m1.0_p16_b1.5',		1.5],
	['m1.0_p16_b2.0',		2.0],
	['m1.0_p16_b3.0',		3.0],
	['m1.0_p16_b4.0',		4.0],
	#['m1.0_p16_b5.0',		5.0]
        ]
]
labels = ['age=0Gyr', 'age=4.8Gyr', 'age=8.4Gyr']
carr = ['C0', 'C1', 'C2']
fig, ax = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()
fig4, ax4 = plt.subplots()
for i, ds in enumerate(dss):
	b_array = []
	ninf_array = []
	mdotpeak_array = []
	tpeak_array = []
	for k, d in enumerate(ds):
		x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/'
				+ d[0].replace(".","_") + '_0040_ev.dat',
				skiprows=1, unpack=True)
		# mdotpeak, tpeak
		mdotpeak_array.append(max(y))
		tpeak_array.append(x[np.argmax(y)])

		# n infinity
		ninf = (y[-2] - y[-1]) / (x[-2] - x[-1])
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


ax.set_xlim(-2, 2)
ax.set_ylim(-4, 4)
ax.axhline(-5/3., ls=':', lw=LWD, c='k')
ax.set_ylabel(r'$n$')
ax.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
ax.legend()
fig.tight_layout()
fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/paper/power_law_index_instantaneous.pdf')

ax2.set_xlim(0.9, 4.1)
ax2.axhline(-5/3., ls=':', lw=LWD, c='k')
ax2.set_ylabel(r'$n_\infty$')
ax2.set_xlabel(r'$\beta$')
ax2.legend()
fig2.tight_layout()
fig2.savefig('/pfs/lawsmith/FLASH4.3/runs/results/paper/power_law_index_infinity.pdf')

ax3.set_xlim(0.9, 4.1)
ax3.set_ylabel(r'$t_{\rm peak}\ \mathrm{[yr]}$')
ax3.set_xlabel(r'$\beta$')
ax3.legend()
fig3.tight_layout()
fig3.savefig('/pfs/lawsmith/FLASH4.3/runs/results/paper/tpeak_vs_beta.pdf')

ax4.set_xlim(0.9, 4.1)
ax4.set_ylabel(r'$\log\ \dot M_{\rm peak}\ {\rm [M_\odot/yr]}$')
ax4.set_xlabel(r'$\beta$')
ax4.legend()
fig4.tight_layout()
fig4.savefig('/pfs/lawsmith/FLASH4.3/runs/results/paper/mdotpeak_vs_beta.pdf')

### for mdotpeak vs tpeak
# ^ above
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
