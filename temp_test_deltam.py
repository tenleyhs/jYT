import yt
import numpy as np
execfile('/groups/dark/lawsmith/jYT/my_settings.py')

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


odata = np.loadtxt('/storage/dark/lawsmith/FLASH4.3/runs/43_b2.0_48k/pruned_sinks_evol.dat', dtype='float64')

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

edata = np.loadtxt('/storage/dark/lawsmith/FLASH4.3/runs/43_b2.0_48k/extras.dat', dtype='float64')
gmpt = G*edata[6]

#for chk in p[2]:
for chk in ['0100']:#['0040','0060','0080']:

    f = '/storage/dark/lawsmith/FLASH4.3/runs/43_b2.0_48k/multitidal_hdf5_chk_' + chk

    pf = yt.load(f)

    tindex = abs(time - pf.current_time.v).argmin()

    myvars = var('sb_mass')

    yt.add_field(("gas","sb_mass"), function=myvars.mesh, take_log=False, force_override=True, units="g")

    sp = pf.sphere("center", (1.0, "Mpc"))
    dm = sp.quantities.total_quantity("sb_mass").in_units('Msun')

    print chk, dm
