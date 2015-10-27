import matplotlib
matplotlib.use('Agg')
import argparse

def csv(value):
	return map(float, value.split(","))

parser = argparse.ArgumentParser(description='Generate a histogram plot of FLASH data.')
parser.add_argument('filename',                                 help='filename', metavar='filename', nargs='+')
parser.add_argument('--parallel',           dest='parallel',    help='run in parallel', default=False, action='store_true')
parser.add_argument('--vars', '-v',         dest='vars',        help='variables to bin over', default=['Density'], type=str, nargs='+')
parser.add_argument('--excludevars', '-ev', dest='excludevars', help='exclude where these variables is greater than some threshold', default=[], type=str, nargs='+')
parser.add_argument('--excludethr', '-et',  dest='excludethr',  help='Threshold values for exclusions', default=[0.], type=float, nargs='+')
parser.add_argument('--excludetype', '-ey', dest='excludetype', help='Threshold type for exclusions', default=['max'], type=str, nargs='+')
parser.add_argument('--log', '-l',          dest='log',         help='take the log of the variable', action='store_true', default=False)
parser.add_argument('--rotstep', '-rs',     dest='rotstep',     help='angle in degrees between rotations', default=0., type=float)
parser.add_argument('--rotdir', '-rd',      dest='rotdir',      help='normal vector of axis to rotate about', nargs='+', default=[1.,0.,0.], type=float)
parser.add_argument('--res', '-r',          dest='res',         help='image resolution', default=512, type=int)
parser.add_argument('--min',                dest='min',         help='minimum value to plot', type=float)
parser.add_argument('--max',                dest='max',         help='maximum value to plot', type=float)
parser.add_argument('--nbins', '-n',        dest='nbins',       help='number of bins to include', type=int, default=100)
parser.add_argument('--negative',           dest='negative',    help='plots -1.0*value', action='store_true', default=False)
parser.add_argument('--silent',             dest='silent',      help='no plot windows generated', action='store_true', default=False)
parser.add_argument('--fileprefix', '-fp',  dest='fileprefix',  help='prefix for output filename', type=str, default='')
parser.add_argument('--outputfile', '-o',   dest='outputfile',  help='explicitly set output filename', type=str, default='')
parser.add_argument('--minradius',  '-mr',  dest='minradius',   help='minimum distance from point mass to consider, scaled to the pericenter distance', type=float, default=0.)
parser.add_argument('--undersample', '-u',  dest='undersample', help='do not include blocks this many level away from maximum or less', type=int, default=0)
parser.add_argument('--oversample', '-ov',  dest='oversample',  help='interpolate maximally refined blocks by this many levels', type=int, default=0)
parser.add_argument('--subsample', '-s',    dest='subsample',   help='sample lower res blocks within this number of levels of maximum refinement', type=int, default=0)
parser.add_argument('--extras', '-ex',      dest='extras',      help='input extras array when extras.dat does not exist', default='1.,1.e6,1.', type=csv)
parser.add_argument('--version', '-ve',     dest='version',     help='version of the MTP software', default=2, type=int)
args = parser.parse_args()

from yt.config import ytcfg; ytcfg["yt","serialize"] = "False"
from yt.mods import *
import matplotlib.pyplot as plt
import string
import numpy as na
from jYT import congrid
from progressbar import Bar, ETA, Percentage, ProgressBar

class var(object):
	def __init__(self, name):
		self.name = name
	def mesh(self, field, data):
		if (self.name == 'kine') :
			return sign*0.5*data['dens']*(data['velx']**2 + data['vely']**2 + data['velz']**2)
		elif (self.name == 'vel2') :
			return sign*(data['velx']**2 + data['vely']**2 + data['velz']**2)
		elif (self.name == 'angmom') :
			posx = data['x'] - ptvec[tindex,0]
			posy = data['y'] - ptvec[tindex,1]
			posz = data['z'] - ptvec[tindex,2]
			velx = data['velx'] - ptvec[tindex,3]
			vely = data['vely'] - ptvec[tindex,4]
			velz = data['velz'] - ptvec[tindex,5]
			return sign*na.sqrt((posy*velz - posz*vely)**2 + (posz*velx - posx*velz)**2 + (posx*vely - posy*velx)**2)
		elif (self.name == 'rp') :
			posx = data['x'] - ptvec[tindex,0]
			posy = data['y'] - ptvec[tindex,1]
			posz = data['z'] - ptvec[tindex,2]
			velx = data['velx'] - ptvec[tindex,3]
			vely = data['vely'] - ptvec[tindex,4]
			velz = data['velz'] - ptvec[tindex,5]
			angmom2 = (posy*velz - posz*vely)**2 + (posz*velx - posx*velz)**2 + (posx*vely - posy*velx)**2

			pos = na.sqrt(posx**2 + posy**2 + posz**2)
			vel2 = velx**2 + vely**2 + velz**2
			ener = (-gmpt/pos + 0.5*vel2)

			a = -0.5*gmpt/ener
			ecc = na.sqrt(-angmom2/(gmpt*a) + 1.0)
			return sign*(-ecc + 1.0)*a
		elif (self.name == 'argperi'):
			posx = data['x'] - ptvec[tindex,0]
			posy = data['y'] - ptvec[tindex,1]
			posz = data['z'] - ptvec[tindex,2]
			velx = data['velx'] - ptvec[tindex,3]
			vely = data['vely'] - ptvec[tindex,4]
			velz = data['velz'] - ptvec[tindex,5]
			pos = na.sqrt(posx**2 + posy**2 + posz**2)
			#angmomx = posy*velz - posz*vely
			#angmomy = posz*velx - posx*velz
			#angmomz = posx*vely - posy*velx
			#eccx = (vely*angmomz - velz*angmomy)/gmpt - posx/pos
			#eccy = (velz*angmomx - velx*angmomz)/gmpt - posy/pos
			#eccz = (velx*angmomy - vely*angmomx)/gmpt - posz/pos
			#ecc = na.sqrt(eccx**2 + eccy**2 + eccz**2)
			#return sign*na.arccos(eccx/ecc)

			vel2 = velx**2 + vely**2 + velz**2
			rv = posx*velx + posy*vely + posz*velz
			eccx = (vel2*posx - rv*velx)/gmpt - posx/pos
			eccy = (vel2*posy - rv*vely)/gmpt - posy/pos
			eccz = (vel2*posz - rv*velz)/gmpt - posz/pos
			ecc = na.sqrt(eccx**2 + eccy**2 + eccz**2)
			#return sign*np.arctan2(eccy/ecc, eccx/ecc)
			return sign*na.arccos(eccx/ecc)
		elif (self.name == 'selfbound'):
			vel2 = na.zeros(data['x'].shape, dtype='float64')
			for i, ax in enumerate(['velx', 'vely', 'velz']) :
				vel2 += (data[ax].v - boundvec[tindex,i+3])**2.
			arr = 0.5*data['gpot'].v + 0.5*vel2
			return arr
		elif (self.name == 'bhbound') :
			pos = na.zeros(data['x'].shape, dtype='float64')
			vel2 = pos.copy()
			for i, ax in enumerate(['x','y','z']) :
				pos += (data[ax].v - ptvec[tindex,i])**2.
			for i, ax in enumerate(['velx', 'vely', 'velz']) :
				vel2 += (data[ax].v - ptvec[tindex,i+3])**2.
			na.sqrt(pos, pos)	
			#arr = sign*data['dens']*(gmpt/pos - 0.5*vel2)
			arr = sign*(-gmpt/pos + 0.5*vel2)
			#arr[pos < args.minradius*peridist] = float("nan")
			return arr
		elif (self.name == 'ni56dens') :
			return sign*data['ni56']*data['dens']
		elif (self.name == 'allenergies') :
			pos = na.zeros(data['x'].shape, dtype='float64')
			vel2 = pos.copy()
			for i, ax in enumerate(['x','y','z']) :
				pos += (data[ax] - ptvec[tindex,i])**2.
			for i, ax in enumerate(['velx', 'vely', 'velz']) :
				vel2 += (data[ax] - ptvec[tindex,i+3])**2.
			na.sqrt(pos, pos)	
			#arr = sign*data['dens']*(gmpt/pos - 0.5*vel2)

			vel2s = na.zeros(data['x'].shape, dtype='float64')
			for i, ax in enumerate(['velx', 'vely', 'velz']) :
				vel2s += (data[ax] - boundvec[tindex,i+3])**2.
			selfbound = 0.5*data['gpot'] + 0.5*vel2s

			arr = sign*(-gmpt/pos + 0.5*vel2 + selfbound + data['eint'])
			#arr[pos < args.minradius*peridist] = float("nan")
			return arr
		else :
			return sign*data[self.name]

g = 6.67428e-8
msun = 1.9889225e33
rsun = 6.955e10

if (args.negative) :
	sign = -1.0
else :
	sign = 1.0

nghost = 1

if (len(args.excludevars) > 0):
	if len(args.excludevars) != len(args.excludethr) or len(args.excludevars) != len(args.excludetype):
		print 'ERROR: Exclude variable names, threshold values, and types must all be the same length.'
		sys.exit()

if (set(['bhbound','selfbound','angmom']) & set(args.vars + args.excludevars)):
	odata = na.loadtxt('pruned_orbit.dat', dtype='float64')
	ptvec = odata[:,1:7]
	obvec = odata[:,7:13]
	boundvec = odata[:,13:19]
	totvec = odata[:,19:25]
	time = odata[:,0]
	if args.version == 2 :
		ptvec += totvec - obvec
	else :
		mpolevec = odata[:,25:31]
		ptvec += totvec - mpolevec

	import os.path
	if os.path.isfile('extras.dat') :
		edata = na.loadtxt('extras.dat', dtype='float64')
		gmpt = g*edata[6]
		peridist = edata[3]*na.power(edata[6]/edata[4]/msun, 1./3.)
	else :
		extras = na.array(args.extras)
		gmpt = g*msun*extras[1]
		peridist = extras[0]*rsun*na.power(extras[1]/extras[2], 1./3.)

myvars = []
for e, ev in enumerate(list(set(args.vars) | set(args.excludevars))):
	myvars += [var(ev)]
	add_field(ev, function=myvars[e].mesh, take_log=args.log)

#my_storage = {}
for f in args.filename:
#for sto, f in parallel_objects(args.filename, 8, storage = my_storage):
	pf = load(f)

	if (set(['bhbound','selfbound','angmom']) & set(args.vars + args.excludevars)):
		odata = na.loadtxt('pruned_orbit.dat', dtype='float64')
		time = odata[:,0]
		tindex = abs(time - pf.current_time.v).argmin()

	if args.subsample >= 0 and pf.h.max_level - args.undersample < args.subsample:
		print 'ERROR: Subsample must be less than max refine level - undersample.'
		sys.exit()

	maxval = na.empty(len(args.vars))
	minval = na.empty(len(args.vars))
	maxval.fill(-float("inf"))
	minval.fill(float("inf"))
	vals = list()
	pbar = ProgressBar(widgets=['Determining histogram bounds and initial pass of data: ', Percentage(), Bar(), ' ', ETA()], maxval=len(pf.index.grids)).start()
	for cnt, g in enumerate(pf.index.grids):
		if g.Level > pf.h.max_level - args.undersample: continue
		if len(g.Children) != 0 and g.Level != pf.h.max_level - args.undersample: continue

		evals = list()
		vvals = list()
#vvals = g.get_data(args.var).ravel()
		for e, ev in enumerate(args.vars):
			vvals.append(g[ev].ravel())
		dvals = g["dens"].ravel()
		for e, ev in enumerate(args.excludevars):
			evals.append(g[ev].ravel())
		for e in range(len(evals)):
			if len(evals[e]) == 0: continue
			if args.excludetype[e] == 'min':
				dvals = dvals[evals[e] > args.excludethr[e]]
				for e2, ev in enumerate(args.vars):
					vvalstemp = vvals[e2]
					vvals[e2] = vvalstemp[evals[e] > args.excludethr[e]]
				for e2, eval2 in enumerate(evals[e+1:]):
					evals[e+e2+1] = eval2[evals[e] > args.excludethr[e]]
			else:
				dvals = dvals[evals[e] < args.excludethr[e]]
				for e2, ev in enumerate(args.vars):
					vvalstemp = vvals[e2]
					vvals[e2] = vvalstemp[evals[e] < args.excludethr[e]]
				for e2, eval2 in enumerate(evals[e+1:]):
					evals[e+e2+1] = eval2[evals[e] < args.excludethr[e]]

		if args.subsample < 0:
			ss = max(pf.h.max_level + args.subsample - g.Level, 0)
		else:
			ss = min(pf.h.max_level - g.Level, args.subsample)

		for e in range(len(vvals)):
			if len(vvals[e]) > 0:
				lmax = na.nanmax(vvals[e])
				lmin = na.nanmin(vvals[e])
				if lmax > maxval[e]: maxval[e] = lmax
				if lmin < minval[e]: minval[e] = lmin

		if ss == 0:
			vals.append([vvals,dvals,na.product(g.dds)])

		pbar.update(cnt+1)
	pbar.finish()

	histbins = list()
	for e in range(len(args.vars)):
		histbins.append(na.linspace(minval[e], maxval[e], args.nbins+1))

	histdims = na.empty(len(args.vars))
	histdims.fill(args.nbins)
	hist = na.zeros(histdims)

	print minval
	print maxval

#Oversampling needs to account for multiple variables
	if args.oversample == 0:
		if len(vals) > 0:
			pbar = ProgressBar(widgets=['Binning raw data: ', Percentage(), Bar(), ' ', ETA()], maxval=len(vals)).start()
			for cnt, val in enumerate(vals):
				if len(val[1]) == 0: continue
				lhist = na.histogramdd(na.transpose(val[0]), bins=histbins, weights=val[1])[0]
				hist += val[2]*lhist
				pbar.update(cnt+1)
			pbar.finish()
		else:
			print 'No cells satisfy cuts in initial pass!'

	pbar = ProgressBar(widgets=['Secondary pass of data with sub-sampling: ', Percentage(), Bar(), ' ', ETA()], maxval=len(pf.index.grids)).start()
	for cnt, g in enumerate(pf.index.grids):
		if g.Level > pf.h.max_level - args.undersample: continue
		if len(g.Children) != 0 and g.Level != pf.h.max_level - args.undersample: continue

		if args.subsample < 0:
			ss = max(pf.h.max_level + args.subsample - g.Level, 0)
		else:
			ss = min(pf.h.max_level - g.Level, args.subsample) + args.oversample

		evals = list()
		vvals = list()
		if (ss != 0):
			excludelist = list()
			for e, ev in enumerate(args.excludevars):
				excludelist.append(ev)
			varstrings = args.vars + ["dens"] + excludelist
			gz = g.retrieve_ghost_zones(nghost, varstrings, smoothed=False)
			ssg = 2**(ss - 1)*nghost
			cgshape = 2**ss*(na.array(na.shape(gz['x'])) - 1) + 1
			ofs = 2.**(-ss-1)
			for e, ev in enumerate(args.vars):
				vvals.append(congrid(gz[ev], cgshape, minusone=True, offset=ofs, buffer=ssg).ravel())
			dvals = congrid(gz["dens"], cgshape, minusone=True, offset=ofs, buffer=ssg).ravel()
			for e, ev in enumerate(args.excludevars):
				evals.append(congrid(gz[ev], cgshape, minusone=True, offset=ofs, buffer=ssg).ravel())
		for e in range(len(evals)):
			if len(evals[e]) == 0: continue
			if args.excludetype[e] == 'min':
				dvals = dvals[evals[e] > args.excludethr[e]]
				for e2, ev in enumerate(args.vars):
					vvalstemp = vvals[e2]
					vvals[e2] = vvalstemp[evals[e] > args.excludethr[e]]
				for e2, eval2 in enumerate(evals[e+1:]):
					evals[e+e2+1] = eval2[evals[e] > args.excludethr[e]]
			else:
				dvals = dvals[evals[e] < args.excludethr[e]]
				for e2, ev in enumerate(args.vars):
					vvalstemp = vvals[e2]
					vvals[e2] = vvalstemp[evals[e] < args.excludethr[e]]
				for e2, eval2 in enumerate(evals[e+1:]):
					evals[e+e2+1] = eval2[evals[e] < args.excludethr[e]]
		if len(vvals) > 0:
			if len(vvals[0]) > 0:
				vol = na.product(g.dds)/2**(3*ss)
				hist += vol*na.histogramdd(na.transpose(vvals), bins=histbins, weights=dvals)[0]
		pbar.update(cnt+1)
	pbar.finish()

	for e in range(len(args.vars)):
		histbinstemp = histbins[e]
		histbins[e] = histbinstemp[0:-1] + (histbinstemp[1:] - histbinstemp[0:-1])/2.

	if args.outputfile != '':
		fn = args.outputfile
	else:
		fn = ('_'.join(args.vars))+'_histogram_'+f+'.dat'
		if args.fileprefix != '':
			fn = args.fileprefix+'_'+fn

	#output histogram here
	f = open(fn, 'w')
	f.write(str(len(args.vars))+'\n')
	f.write(str(args.nbins)+'\n')
	f.write(str(na.sum(hist))+'\n')	
	f.write(str(time[tindex])+'\n')
	if len(args.vars) == 1:
		f.write(string.join(["% 15.10E" % x for x in histbins[0]])+'\n')
		f.write(string.join(["% 15.10E" % x for x in hist])+'\n')
	else:
		f.write(string.join(["% 15.10E" % x for x in [item for sublist in histbins for item in sublist]])+'\n')
		f.write(string.join(["% 15.10E" % x for x in [item for sublist in hist for item in sublist]])+'\n')
	f.close()

	if not args.silent:
		lhist = na.log10(hist)
		
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.plot(histbins[0], lhist)
		fig.savefig(('_'.join(args.vars))+'_histogram_'+str(f)+'.png')


