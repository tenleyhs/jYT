"""
first writes, then can load and use .dat files for later plotting
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

USE_DAT = False

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
	['m1.0_p10_b5.0',		5.0]
        ],
        [
	['m1.0_p16_b1.0',		1.0],
	['m1.0_p16_b1.5',		1.5],
	['m1.0_p16_b2.0',		2.0],
	['m1.0_p16_b3.0',		3.0],
	['m1.0_p16_b4.0',		4.0],
	['m1.0_p16_b5.0',		5.0]
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

if USE_DAT:
    for i, name in enumerate(names):
        b, log_deltam_m = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/deltaM_vs_beta_' + names[i] + '.dat',
        unpack=True, skiprows=1)

        ax.plot(b, log_deltam_m, lw=lw, label=labels[i], alpha=0.75)
        ax.scatter(b, log_deltam_m, s=s)

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
            names=['beta','log_deltaM_Msun'])

ax.set_xlim(0.9, 4.1)
ax.set_ylim(-1.5, 0.05)
ax.set_xlabel(r'$\beta$')
ax.set_ylabel(r'$\log\ \Delta M/M_\ast$')
ax.legend()
fig.tight_layout()
fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/paper/deltaM_vs_beta.pdf')
