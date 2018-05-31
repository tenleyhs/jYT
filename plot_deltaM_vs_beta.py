"""
first writes, then can load and use .dat files for later plotting
TODO add fitting function
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import ascii
import yt
execfile('/pfs/lawsmith/jYT/my_settings.py')
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['font.size'] = 18
from scipy.optimize import curve_fit

USE_DAT = True

s = 30
lw = 2

# note: have to change a bit if comparing 256
ages = [
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

# old. using different chks. maybe compare later.
"""
ages = [
        [
    ['m1.0_p1_b1.0',1.0,'multitidal_hdf5_chk_0100'],
    ['m1.0_p1_b2.0',2.0,'multitidal_hdf5_chk_0080'],
    ['m1.0_p1_b3.0',3.0,'multitidal_hdf5_chk_0060']
        ],
        [
    ['m1.0_p10_b1.0',1.0,'multitidal_hdf5_chk_0100'],
    ['m1.0_p10_b1.5',1.5,'multitidal_hdf5_chk_0090'],
    ['m1.0_p10_b2.0_128',2.0,'multitidal_hdf5_chk_0075'],
    #['m1.0_p10_b2.5',2.5,'multitidal_hdf5_chk_0080'],
    ['m1.0_p10_b3.0',3.0,'multitidal_hdf5_chk_0060'],
    ['m1.0_p10_b4.0',4.0,'multitidal_hdf5_chk_0050'],
    #['m1.0_p10_b5.0',5.0,'multitidal_hdf5_chk_0060']
        ],
        [
    ['m1.0_p16_b2.0',2.0,'multitidal_hdf5_chk_0075'],
    ['m1.0_p16_b3.0',3.0,'multitidal_hdf5_chk_0060'],
    ['m1.0_p16_b4.0',4.0,'multitidal_hdf5_chk_0050'],
    ['m1.0_p16_b5.0',5.0,'multitidal_hdf5_chk_0050']
        ]
       ]
"""

labels = ['age=0Gyr', 'age=4.8Gyr', 'age=8.4Gyr']
names = ['p0', 'p10', 'p16']

# fitting function. same as GRR2013.
# TODO need at least as many data points as parameters
#def f(x, a, b, c, d, e, f):
#    return np.exp((a + b*x + c*x**2)/(d + e*x + f*x**2))

# my temporary function with 5 params
def f5(x, a, b, c, d, e):
    return a + b*x**0.25 + c*x**0.5 + d*x**0.75 + e*x**1.0

# my temporary function with 4 params
def f4(x, a, b, c, d):
    return a + b*x**0.25 + c*x**0.5 + d*x**0.75


def C_43(b):
    return np.exp((12.996 - 31.149*b + 12.865*b**2)/(1 - 5.3232*b + 6.4262*b**2))

class var(object):
    def __init__(self, name):
        self.name = name
    def mesh(self, field, data):
        if (self.name == 'sb_mass') :
            vel2 = np.zeros(data['x'].shape, dtype='float64')
            for i, ax in enumerate(['velx', 'vely', 'velz']) :
                vel2 += (data[ax].v - obvec[tindex,i+3])**2.
            arr = 0.5*data['gpot'].v + 0.5*vel2
            arr[arr >= 0.] = 1.
            arr[arr < 0.] = 0.
            return arr * data['cell_mass']


fig, ax = plt.subplots()
fig2, ax2 = plt.subplots() # for shifted delta m
fig3, ax3 = plt.subplots() # for shifted delta m ratio

if USE_DAT:
    # plot Guillochon2013
    # Delta M = C_\gamma M_\star
    b43 = np.linspace(0.6, 1.85)
    ax.plot(b43, np.log10(C_43(b43)), lw=lw, ls='--', label='GRR13 4/3', c='k', alpha=0.5)
    ax.plot([1.86, 5], [0, 0], lw=lw, ls='--', c='k', alpha=0.5)

    ax2.plot(b43, np.log10(C_43(b43)), lw=lw, ls='--', label='GRR13 4/3', c='k', alpha=0.5)
    ax2.plot([1.86, 5], [0, 0], lw=lw, ls='--', c='k', alpha=0.5)

    ax3.plot(b43/1.85, np.log10(C_43(b43)), lw=lw, ls='--', label='GRR13 4/3', c='k', alpha=0.5)
    ax3.plot([1.86, 5], [0, 0], lw=lw, ls='--', c='k', alpha=0.5)

    # for comparing mdots at fixed delta ms
    fixed_delta_ms = [-0.5, -0.25, -0.1]

    for i, name in enumerate(names):
        b, log_deltam_m = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/deltaM_vs_beta_' + names[i] + '.dat',
        unpack=True, skiprows=1)

        #ax.plot(b, log_deltam_m, lw=lw, label=labels[i], alpha=0.75)
        ax.scatter(b, log_deltam_m, s=s)

        # fitting
        #popt, pcov = curve_fit(f, b, 10**log_deltam_m)
        xdata = np.linspace(min(b), max(b), num=500)
        # to avoid overfitting, wild oscillation
        if len(b) == 5: f = f4
        if len(b) >= 6: f = f5
        popt, pcov = curve_fit(f, b, log_deltam_m)
        #print names[i], popt
        #ax.plot(xdata, np.log10(f(xdata, *popt)))
        ax.plot(xdata, f(xdata, *popt), label=labels[i])

        # for comparing mdots at fixed delta m
        # choose one or a few delta m at which to compare
        # maybe log delta m/m -0.5, -0.25, -0.1
        for dm in fixed_delta_ms:
            idx = (np.abs(f(xdata, *popt) - dm)).argmin()
            this_fixed_beta = xdata[idx]
            print names[i], dm, this_fixed_beta
            # TODO will want to save this to .dat file later



        # TODO not real critical beta right now.
        # for shifted delta M
        # second number is critical beta for this profile. 1.85 is for GRR13 4/3
        if names[i] == 'p0':
            bc = 1.7
        if names[i] == 'p10':
            bc = 3
        if names[i] == 'p16':
            bc = 4
        shift = 1.85 - bc
        ax2.scatter(b + shift, log_deltam_m, s=s)
        ax2.plot(xdata + shift, f(xdata, *popt), label=labels[i])

        ax3.scatter(b/bc, log_deltam_m, s=s)
        ax3.plot(xdata/bc, f(xdata, *popt), label=labels[i])



else:
    for i, age in enumerate(ages):
        b_array = []
        dm_array = []
        for p in age:
            odata = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/' + p[0] + '/pruned_sinks_evol.dat', dtype='float64')

            if odata[0,14] > odata[1,14]:
                part_tag_pt = odata[0,0]
            else:
                part_tag_pt = odata[1,0]
            odata_pt = odata[np.where(odata[:,0]==part_tag_pt)[0]]
            odata_ob = odata[np.where(odata[:,0]!=part_tag_pt)[0]]

            ptvec = odata_pt[:,2:8]
            obvec = odata_ob[:,2:8]
            boundvec = obvec
            time = odata_pt[:,1]

            edata = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/' + p[0] + '/extras.dat', dtype='float64')
            gmpt = G*edata[6]

            #f = '/pfs/lawsmith/FLASH4.3/runs/' + p[0] + '/' + p[2]
            f = '/pfs/lawsmith/FLASH4.3/runs/' + p[0] + '/multitidal_hdf5_chk_0040'

            pf = yt.load(f)

            tindex = abs(time - pf.current_time.v).argmin()

            myvars = var('sb_mass')

            yt.add_field(("gas","sb_mass"), function=myvars.mesh, take_log=False, force_override=True, units="g")

            # todo might want to add density cut later & see if affects results.
            #ad = pf.all_data()
            #dense_ad = ad.cut_region(['obj["dens"] > 1e-11'])

            sp = pf.sphere("center", (1.0, "Mpc"))
            dm = sp.quantities.total_quantity("sb_mass").v

            b_array.append(p[1])
            dm_array.append(dm)

        ax.plot(b_array, np.log10(np.array(dm_array)/M_sun), lw=lw, label=labels[i], alpha=0.75)
        ax.scatter(b_array, np.log10(np.array(dm_array)/M_sun), s=s)
        ascii.write([b_array, np.log10(np.array(dm_array)/M_sun)],
            '/pfs/lawsmith/FLASH4.3/runs/results/deltaM_vs_beta_' + names[i] + '.dat',
            names=['beta','log_deltaM_Msun'], overwrite=True)

ax.set_xlim(0.9, 4.1)
ax.set_ylim(-1.5, 0.05)
ax.set_xlabel(r'$\beta$')
ax.set_ylabel(r'$\log\ \Delta M/M_\ast$')
ax.legend()
fig.tight_layout()
fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/paper/deltaM_vs_beta.pdf')

ax2.set_xlim(-1.4, 1.9)
ax2.set_ylim(-1.5, 0.05)
ax2.set_xlabel(r'$\beta$' + ' (shifted)')
ax2.set_ylabel(r'$\log\ \Delta M/M_\ast$')
ax2.legend(loc='upper left')
fig2.tight_layout()
fig2.savefig('/pfs/lawsmith/FLASH4.3/runs/results/paper/deltaM_vs_beta_shifted.pdf')

ax3.set_xlim(0, 1.1)
ax3.set_ylim(-1.5, 0.05)
ax3.set_xlabel(r'$\beta$' + ' (normalized to ' + r'$\beta_c$' + ')')
ax3.set_ylabel(r'$\log\ \Delta M/M_\ast$')
ax3.legend(loc='upper left')
fig3.tight_layout()
fig3.savefig('/pfs/lawsmith/FLASH4.3/runs/results/paper/deltaM_vs_beta_shifted2.pdf')
