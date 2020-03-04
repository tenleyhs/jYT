"""
idea is to find last checkpoint at which all the debris is in the box. use like:

python ~/jYT/temp_test_deltam.py --runs m0.3_p1_b0.6_300k m0.3_p1_b0.7_300k m0.3_p1_b0.8_300k m0.3_p1_b0.9_300k m0.3_p1_b1.0_300k

actually it looks like I switched to using this for testing whether delta m changes from chk 40, 60, 80, 100, etc.

TODO I stopped editing this file b/c I was able to just look at slices for m0.3_p1..., should go through it line by line if picking up again.
"""

import yt
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--cluster', help='cluster, e.g., fend/lux/pfe', type=str, default='fend')
parser.add_argument('--runs', help='runs, e.g., m1.0_p16_b3.0_300k', type=str, nargs='+')
parser.add_argument('--chks', help='checkpoints, e.g., 0100', type=str, nargs='+')
args = parser.parse_args()



if args.cluster == 'fend':
    execfile('/groups/dark/lawsmith/jYT/my_settings.py')
    rundir = '/storage/dark/lawsmith/FLASH4.3/runs/'
    savedir = '/groups/dark/lawsmith/results/delta_m/'
    parfile = '/groups/dark/lawsmith/jYT/deltaM_fend.par'
elif args.cluster == 'lux':
    execfile('/home/lawsmith/jYT/my_settings.py')
    rundir = '/data/groups/ramirez-ruiz/lawsmith/FLASH4.3/'
    savedir = '/data/groups/ramirez-ruiz/lawsmith/FLASH4.3/results/delta_m/'
    parfile = '/home/lawsmith/jYT/deltaM_lux.par'
elif args.cluster == 'pfe':
    rundir = '/nobackup/jlawsmit/'

#LOAD_FILES = rundir + args.run + '/bondi_hdf5_chk_*' + args.files
#LOAD_FILES = rundir + args.run + '/bondi_hdf5_plt_cnt_*'
#yt.enable_parallelism()
#ts = yt.DatasetSeries(LOAD_FILES)

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

outputfile = open(savedir + 'sb.dat', 'a')

runs, chks = np.loadtxt(parfile, unpack=True, dtype='str')

#for (run, chk) in zip(args.runs, args.chks):
for (run, chk) in zip(runs, chks):
    f = rundir + run + '/multitidal_hdf5_chk_' + chk
    pf = yt.load(f)

    odata = np.loadtxt(rundir + run + '/pruned_sinks_evol.dat', dtype='float64')

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

    tindex = abs(time - pf.current_time.v).argmin()

    myvars = var('sb_mass')
    yt.add_field(("gas","sb_mass"), function=myvars.mesh, take_log=False, force_override=True, units="g")
    sp = pf.sphere("center", (1.0, "Mpc"))
    sb = sp.quantities.total_quantity("sb_mass").in_units('Msun').v

    print run, chk, sb
    outputfile.write(run + '\t\t' + chk + '\t\t' + str(sb) + '\n')

outputfile.close()
