import yt
import numpy as np
from yt.fields.api import ValidateParameter
from yt.config import ytcfg; ytcfg["yt","serialize"] = "False"

g = 6.67428e-8
msun = 1.9889225e33
rsun = 6.955e10

runs = [
    ['m1.0_p1_b2.0', '0208', 'p1_b20_0208'],
    ['m1.0_p10_b2.0_256', '0207', 'p10_b20_256_0207'],
    ['m1.0_p10_b2.0_256', '0800', 'p10_b20_256_0800']
]

class var(object):
    def __init__(self, name):
        self.name = name
    def mesh(self, field, data):
        if (self.name == 'fb_time') :
            pos = np.zeros(data['x'].shape, dtype='float64')
            vel2 = pos.copy()
            for i, ax in enumerate(['x','y','z']) :
                pos += (data[ax].v - ptvec[tindex,i])**2.
            for i, ax in enumerate(['velx', 'vely', 'velz']) :
                vel2 += (data[ax].v - ptvec[tindex,i+3])**2.
            np.sqrt(pos, pos)
            #arr = sign*data['dens']*(gmpt/pos - 0.5*vel2)
            energy = sign*(-gmpt/pos + 0.5*vel2)
            #energyArray.append(energy).ravel()
            #arr[pos < args.minradius*peridist] = float("nan")
            time = (2*np.pi*gmpt)/np.power(-2*energy, 3./2.)
            #timeArray.append(time).ravel()
            return time

sign = 1

for run in runs:
    odata = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/'+run[0]+'/pruned_sinks_evol.dat', dtype='float64')


    # are the part_tags always the same? or can just sort by mass
    # pt is black hole, ob is star
    if odata[0,14] > odata[1,14]:
        part_tag_pt = odata[0,0]
    else:
        part_tag_pt = odata[1,0]
    odata_pt = odata[np.where(odata[:,0]==part_tag_pt)[0]]
    odata_ob = odata[np.where(odata[:,0]!=part_tag_pt)[0]]

    ptvec = odata_pt[:,2:8]
    obvec = odata_ob[:,2:8]
    # just setting self-bound mass vector to object vector. should be very similar
    boundvec = obvec
    time = odata_pt[:,1]

    edata = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/'+run[0]+'/extras.dat', dtype='float64')
    gmpt = g*edata[6]
    peridist = edata[3]*np.power(edata[6]/edata[4]/msun, 1./3.)


    f = '/pfs/lawsmith/FLASH4.3/runs/'+run[0]+'/multitidal_hdf5_plt_cnt_'+run[1]

    pf = yt.load(f)

    tindex = abs(time - pf.current_time.v).argmin()

    myvars = var('fb_time')

    yt.add_field(("gas","fb_time"), function=myvars.mesh, take_log=False, force_override=True)

    ad = pf.all_data()
    dense_ad = ad.cut_region(['obj["dens"] > 1e-8'])

    s = yt.SlicePlot(pf, 'z', 'fb_time', data_source=dense_ad, width=(8,'rsun'))
    s.set_cmap(field="fb_time", cmap='jet')
    s.set_log('fb_time', True)
    s.set_zlim('fb_time',1e5, 1e8)
    #s.annotate_timestamp(text_args={'color':'black'})
    s.save('/pfs/lawsmith/FLASH4.3/runs/results/eas/eas_'+run[2]+'.pdf')


    s = yt.SlicePlot(pf, 'z', 'dens', width=(8,'rsun'))
    s.save('/pfs/lawsmith/FLASH4.3/runs/results/eas/dens_'+run[2]+'.pdf')
