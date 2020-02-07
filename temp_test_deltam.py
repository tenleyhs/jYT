"""
idea is to find last checkpoint at which all the debris is in the box. use like:

python ~/jYT/temp_test_deltam.py --runs m0.3_p1_b0.6_300k m0.3_p1_b0.7_300k m0.3_p1_b0.8_300k m0.3_p1_b0.9_300k m0.3_p1_b1.0_300k

TODO I stopped editing this file b/c I was able to just look at slices for m0.3_p1..., should go through it line by line if picking up again.
"""

import yt
import numpy as np
execfile('/groups/dark/lawsmith/jYT/my_settings.py')

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--cluster', help='cluster, e.g., fend/lux/pfe', type=str, default='fend')
parser.add_argument('--runs', help='runs, e.g., m1.0_p16_b3.0_300k', type=str, nargs='+')
parser.add_argument('--chks', help='checkpoints, e.g., 0100', type=str)
args = parser.parse_args()

if args.cluster == 'fend':
    clusterdir = '/storage/dark/lawsmith/FLASH4.3/runs/'
    savepath = '/groups/dark/lawsmith/fend_slices/'
elif args.cluster == 'lux':
    clusterdir = '/data/groups/ramirez-ruiz/lawsmith/FLASH4.3/'
    savepath = '/data/groups/ramirez-ruiz/lawsmith/FLASH4.3/lux_slices/'
elif args.cluster == 'pfe':
    clusterdir = '/nobackup/jlawsmit/'
    savepath = '/nobackup/jlawsmit/pfe_slices/'
else:
    print("haven't set up this cluster in script yet")
    exit()

#if not os.path.exists(savepath + args.run + '/' + args.var + '_' + str(args.width) + 'rsun/'):
#	os.makedirs(savepath + args.run + '/' + args.var + '_' + str(args.width) + 'rsun/')

chks = ['0100'] #['0040','0060','0080']

if args.chkplt == 'chk':
    LOAD_FILES = clusterdir + args.run + '/multitidal_hdf5_chk_' + args.files
elif args.chkplt == 'plt':
    LOAD_FILES = clusterdir + args.run + '/multitidal_hdf5_plt_cnt_' + args.files
#LOAD_FILES = clusterdir + args.run + '/bondi_hdf5_chk_*'
#LOAD_FILES = clusterdir + args.run + '/bondi_hdf5_plt_cnt_*'
yt.enable_parallelism()
ts = yt.DatasetSeries(LOAD_FILES)

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

outputfile = open(fn, 'a')


for run in args.runs:
    outputfile.write(run+'\n')

    for chk in args.chks:

        f = clusterdir + run + '/multitidal_hdf5_chk_' + chk

        pf = yt.load(f)

        tindex = abs(time - pf.current_time.v).argmin()

        myvars = var('sb_mass')

        yt.add_field(("gas","sb_mass"), function=myvars.mesh, take_log=False, force_override=True, units="g")

        sp = pf.sphere("center", (1.0, "Mpc"))
        dm = sp.quantities.total_quantity("sb_mass").in_units('Msun')

        print chk, dm

outputfile.close()
