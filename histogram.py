import argparse

parser = argparse.ArgumentParser(description='Generate a histogram plot of FLASH data.')
parser.add_argument('filename',                                 help='filename', metavar='filename', nargs='+')
parser.add_argument('--parallel',           dest='parallel',    help='run in parallel', default=False, action='store_true')
parser.add_argument('--var', '-v',          dest='var',         help='variable', default='Density', type=str)
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
args = parser.parse_args()

from yt.config import ytcfg; ytcfg["yt","serialize"] = "False"
from yt.mods import *
import pylab as pl
import string
from jYT import congrid
from progressbar import Bar, ETA, Percentage, ProgressBar

g = 6.67428e-8
msun = 1.9889225e33

if (args.negative) :
	sign = -1.0
else :
	sign = 1.0

nghost = 1

if (len(args.excludevars) > 0):
	if len(args.excludevars) != len(args.excludethr) or len(args.excludevars) != len(args.excludetype):
		print 'ERROR: Exclude variable names, threshold values, and types must all be the same length.'
		sys.exit()

if (set(['bhbound','selfbound']) & set([args.var] + args.excludevars)):
	odata = na.loadtxt('pruned_orbit.dat', dtype='float64')
	ptvec = odata[:,1:7]
	obvec = odata[:,7:13]
	totvec = odata[:,19:25]
	peakvec = odata[:,25:31]
	time = odata[:,0]

	ptvec += totvec - obvec

	edata = na.loadtxt('extras.dat', dtype='float64')
	gmpt = g*edata[6]
	peridist = edata[3]*na.power(edata[6]/edata[4]/msun, 1./3.)

def loadvar(name, fieldname):
	if (name == 'kine') :
		def _var(field, data):
			return sign*0.5*data['dens']*(data['velx']**2 + data['vely']**2 + data['velz']**2)
	elif (name == 'vel2') :
		def _var(field, data):
			return sign*(data['velx']**2 + data['vely']**2 + data['velz']**2)
	elif (name == 'selfbound'):
		def _var(field, data):
			vel2 = na.zeros(data['x'].shape, dtype='float64')
			for i, ax in enumerate(['velx', 'vely', 'velz']) :
				vel2 += (data[ax] - peakvec[tindex,i+3])**2.
			arr = 0.5*data['gpot'] + 0.5*vel2
			return arr
	elif (name == 'bhbound') :

		def _var(field, data):
			pos = na.zeros(data['x'].shape, dtype='float64')
			vel2 = pos.copy()
			for i, ax in enumerate(['x','y','z']) :
				pos += (data[ax] - ptvec[tindex,i])**2.
			for i, ax in enumerate(['velx', 'vely', 'velz']) :
				vel2 += (data[ax] - ptvec[tindex,i+3])**2.
			na.sqrt(pos, pos)	
			#arr = sign*data['dens']*(gmpt/pos - 0.5*vel2)
			arr = sign*(gmpt/pos - 0.5*vel2)
			#arr[pos < args.minradius*peridist] = float("nan")
			return arr
	else :
		def _var(field, data):
			return sign*data[name]
	add_field(fieldname, function=_var, take_log=args.log)

loadvar(args.var, "var")

for e, ev in enumerate(args.excludevars):
	loadvar(ev, "excludevar"+str(e))

#my_storage = {}
for f in args.filename:
#for sto, f in parallel_objects(args.filename, 8, storage = my_storage):
	pf = load(f)

	odata = na.loadtxt('pruned_orbit.dat', dtype='float64')
	time = odata[:,0]

	tindex = abs(time - pf.current_time).argmin()

	if args.subsample >= 0 and pf.h.max_level - args.undersample < args.subsample:
		print 'ERROR: Subsample must be less than max refine level - undersample.'
		sys.exit()

	maxval = -float("inf")
	minval = float("inf")
	vals = list()
	pbar = ProgressBar(widgets=['Determining histogram bounds and initial pass of data: ', Percentage(), Bar(), ' ', ETA()], maxval=len(pf.h.grids)).start()
	for cnt, g in enumerate(pf.h.grids):
		if g.Level > pf.h.max_level - args.undersample: continue
		if len(g.Children) != 0 and g.Level != pf.h.max_level - args.undersample: continue

		evals = list()
		vvals = g.get_data("var").ravel()
		dvals = g.get_data("dens").ravel()
		for e in range(len(args.excludevars)):
			evals.append(g.get_data("excludevar"+str(e)).ravel())
		for e in range(len(evals)):
			if len(evals[e]) == 0: continue
			if args.excludetype[e] == 'min':
				dvals = dvals[evals[e] > args.excludethr[e]]
				vvals = vvals[evals[e] > args.excludethr[e]]
				for e2, eval2 in enumerate(evals[e+1:]):
					evals[e+e2+1] = eval2[evals[e] > args.excludethr[e]]
			else:
				dvals = dvals[evals[e] < args.excludethr[e]]
				vvals = vvals[evals[e] < args.excludethr[e]]
				for e2, eval2 in enumerate(evals[e+1:]):
					evals[e+e2+1] = eval2[evals[e] < args.excludethr[e]]
		if len(vvals) > 0:
			if args.subsample < 0:
				ss = max(pf.h.max_level + args.subsample - g.Level, 0)
			else:
				ss = min(pf.h.max_level - g.Level, args.subsample) + args.oversample

			lmax = na.nanmax(vvals)
			lmin = na.nanmin(vvals)
			if lmax > maxval: maxval = lmax
			if lmin < minval: minval = lmin
			if ss == 0:
				vals.append([vvals,dvals,na.product(g.dds)])
		pbar.update(cnt+1)
	pbar.finish()

	histbins = na.arange(minval, maxval, (maxval-minval)/(args.nbins + 1))

	hist = na.zeros(len(histbins) - 1)
	if len(vals) > 0:
		pbar = ProgressBar(widgets=['Binning raw data: ', Percentage(), Bar(), ' ', ETA()], maxval=len(vals)).start()
		for cnt, val in enumerate(vals):
			hist += val[2]*na.histogram(val[0], bins=histbins, weights=val[1])[0]
			pbar.update(cnt+1)
		pbar.finish()
	else:
		print 'No cells satisfy cuts in initial pass!'

	pbar = ProgressBar(widgets=['Secondary pass of data with sub-sampling: ', Percentage(), Bar(), ' ', ETA()], maxval=len(pf.h.grids)).start()
	for cnt, g in enumerate(pf.h.grids):
		if g.Level > pf.h.max_level - args.undersample: continue
		if len(g.Children) != 0 and g.Level != pf.h.max_level - args.undersample: continue

		if args.subsample < 0:
			ss = max(pf.h.max_level + args.subsample - g.Level, 0)
		else:
			ss = min(pf.h.max_level - g.Level, args.subsample) + args.oversample

		evals = list()
		if (ss != 0):
			excludelist = list()
			for e in range(len(args.excludevars)):
				excludelist.append("excludevar"+str(e))
			varstrings = ["var", "dens"] + excludelist
			gz = g.retrieve_ghost_zones(nghost, varstrings, smoothed=False)
			ssg = 2**(ss - 1)*nghost
			cgshape = 2**ss*(na.array(na.shape(gz["var"])) - 1) + 1
			ofs = 2.**(-ss-1)
			vvals = congrid(gz["var"], cgshape, minusone=True, offset=ofs, buffer=ssg).ravel()
			dvals = congrid(gz["dens"], cgshape, minusone=True, offset=ofs, buffer=ssg).ravel()
			for e in range(len(args.excludevars)):
				evals.append(congrid(gz["excludevar"+str(e)], cgshape, minusone=True, offset=ofs, buffer=ssg).ravel())
		for e in range(len(evals)):
			if len(evals[e]) == 0: continue
			if args.excludetype[e] == 'min':
				dvals = dvals[evals[e] > args.excludethr[e]]
				vvals = vvals[evals[e] > args.excludethr[e]]
				for e2, eval2 in enumerate(evals[e+1:]):
					evals[e+e2+1] = eval2[evals[e] > args.excludethr[e]]
			else:
				dvals = dvals[evals[e] < args.excludethr[e]]
				vvals = vvals[evals[e] < args.excludethr[e]]
				for e2, eval2 in enumerate(evals[e+1:]):
					evals[e+e2+1] = eval2[evals[e] < args.excludethr[e]]
		if len(vvals) > 0:
			vol = na.product(g.dds)/2**(3*ss)
			hist += vol*na.histogram(vvals, bins=histbins, weights=dvals)[0]
		pbar.update(cnt+1)
	pbar.finish()

	histbins = histbins[0:-1] + (histbins[1:] - histbins[0:-1])/2.

	if args.outputfile != '':
		fn = args.outputfile
	else:
		fn = args.var+'_histogram_'+f+'.dat'
		if args.fileprefix != '':
			fn = args.fileprefix+'_'+fn

	#output histogram here
	f = open(fn, 'w')
	f.write(str(args.nbins)+'\n')
	f.write(str(na.sum(hist))+'\n')	
	f.write(str(time[tindex])+'\n')
	f.write(string.join(["% 15.10E" % x for x in -histbins])+'\n')
	f.write(string.join(["% 15.10E" % x for x in hist])+'\n')
	f.close()

	if not args.silent:
		lhist = na.log10(hist)
		pl.plot(histbins, lhist)
		pl.draw()
		pl.show()

