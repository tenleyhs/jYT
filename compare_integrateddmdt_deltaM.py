"""
compared integrated dmdts with delta M's
"""
import numpy as np

dint = '/pfs/lawsmith/FLASH4.3/runs/results/dmdts/integrals/'
ddm = '/pfs/lawsmith/FLASH4.3/runs/results/deltam_vs_beta/'
ps = ['p1', 'p10', 'p16']

for p in ps:
    bs, log_deltam_m = np.loadtxt(ddm + 'deltaM_vs_beta_' + p + '.dat',
                        unpack=True, skiprows=1)
    print
    print p
    for b, ldm in zip(bs, log_deltam_m):
        int, int2 = np.loadtxt(dint + 'm1_0_' + p + '_b' + str(b).replace(".","_") + '_0040_ev.dat',
                        unpack=True, skiprows=1)
        print b, round(int*2, 3), round(10**ldm,3)
