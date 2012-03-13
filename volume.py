import argparse

parser = argparse.ArgumentParser(description='Generate a volumetric plot of FLASH data.')
parser.add_argument('filename',                                help='filename', metavar='filename', nargs='+')
parser.add_argument('--parallel',          dest='parallel',    help='run in parallel', default=False, action='store_true')
parser.add_argument('--var', '-v',         dest='var',         help='variable', default='Density', type=str)
parser.add_argument('--log', '-l',         dest='log',         help='take the log of the variable', action='store_true', default=False)
parser.add_argument('--rotstep', '-rs',    dest='rotstep',     help='angle in degrees between rotations', default=0., type=float)
parser.add_argument('--rotdir', '-rd',     dest='rotdir',      help='normal vector of axis to rotate about', nargs='+', default=[1.,0.,0.], type=float)
parser.add_argument('--res', '-r',         dest='res',         help='image resolution', default=512, type=int)
parser.add_argument('--min',               dest='min',         help='minimum value to plot', type=float)
parser.add_argument('--max',               dest='max',         help='maximum value to plot', type=float)
parser.add_argument('--center', '-c',      dest='center',      help='which point to center on', choices=('min', 'max', 'geometric', 'custom'), default='geometric')
parser.add_argument('--centervar', '-cv',  dest='centervar',   help='variable to use when centering', type=str, default='var')
parser.add_argument('--width', '-w',       dest='width',       help='width of region to plot', type=float)
parser.add_argument('--sub_samples', '-s', dest='sub_samples', help='sub_samples parameter for camera', type=int, default=5)
parser.add_argument('--no_ghost', '-ng',   dest='no_ghost',    help='include ghost cells (slower)', action='store_true', default=False)
parser.add_argument('--colormap', '-cm',   dest='colormap',    help='colormap used for plot', type=str, default='jet')
parser.add_argument('--ncontours', '-n',   dest='ncontours',   help='number of contours to draw', type=int, default=6)
parser.add_argument('--negative',          dest='negative',    help='plots -1.0*value', action='store_true', default=False)
parser.add_argument('--fileprefix', '-fp', dest='fileprefix',  help='prefix for output filename', type=str, default='')
args = parser.parse_args()

from yt.mods import *
from yt.visualization.volume_rendering.camera import PerspectiveCamera

up = [0.,1.,0.]
L = na.array([0.,0.,1.])

if (args.negative) :
	sign = -1.0
else :
	sign = 1.0

fieldname = "var"
if (args.var == 'kine') :
	def _var(field, data):
		return sign*0.5*data['dens']*(data['velx']**2 + data['vely']**2 + data['velz']**2)
elif (args.var == 'vel2') :
	def _var(field, data):
		return sign*(data['velx']**2 + data['vely']**2 + data['velz']**2)
elif (args.var == 'bhbound') :
	g = 6.67428e-8

	odata = na.loadtxt('pruned_orbit.dat', dtype='float64')
	ptvec = odata[:,1:7]
	obvec = odata[:,7:13]
	totvec = odata[:,19:25]
	time = odata[:,0]

	ptvec += totvec - obvec

	gmpt = g*na.loadtxt('extras.dat', dtype='float64')[6]

	def _var(field, data):
		tindex = abs(time - pf.current_time).argmin()
		pos = na.zeros(data['x'].shape, dtype='float64')
		vel2 = pos.copy()
		for i, ax in enumerate('xyz') :
			r = na.abs(data[ax] - ptvec[tindex,i])
			pos += r**2.0
		for i, ax in enumerate(['velx', 'vely', 'velz']) :
			v = na.abs(data[ax] - ptvec[tindex,i+3])
			vel2 += v**2.0
		na.sqrt(pos, pos)	
		arr = sign*data['dens']*(gmpt/pos - 0.5*vel2)
		#arr = sign*(gmpt/pos - 0.5*vel2)
		if (args.log == True) :
			if (args.min != None) :
				arr[arr < 0] = args.min/100.
			else :
				arr[arr < 0] = na.finfo(na.float32).tiny
		return arr
else :
	def _var(field, data):
		return sign*data[args.var]

add_field(fieldname, function=_var, take_log=args.log)

# Load up your dataset (CHANGE THIS)
for f in args.filename:
	pf = load(f)

	if (args.width == None) :
		W = sqrt(sum((pf.domain_right_edge - pf.domain_left_edge)**2))
	else :
		W = args.width
		
	if args.center == 'geometric' :
		c = na.array((pf.domain_right_edge - pf.domain_left_edge) / 2.)
	elif args.center == 'max' :
		v,cmax = pf.h.find_max(args.centervar)
		c = na.array(cmax)
	elif args.center == 'min' :
		v,cmin = pf.h.find_min(args.centervar)
		c = na.array(cmin)

	sp = pf.h.sphere(c, W)
	mi, ma = sp.quantities['Extrema']('var')[0]
	print mi, ma

	if (args.min != None) :
		mi = args.min
	if (args.max != None) :
		ma = args.max

	if (args.log) : ma = na.log10(ma)
	if (args.log) :
		if mi < 0. :
			print 'Warning: Minimum less than 0 for log variable! Setting minimum to 0.'
			mi = 0
		else :
			mi = na.log10(mi)
	
	# Construct transfer function, pad the TF space by a bit so that
	# gaussians sampling the data range don't hit the edge.
	tf = ColorTransferFunction((mi-2, ma+2), nbins=1024)
	
	# Sample transfer function with several gaussians.  Use col_bounds keyword
	# to restrict color mapping to true data range.
	tf.add_layers(args.ncontours, w=0.001, alpha=na.logspace(-0.5,0,args.ncontours), col_bounds = (mi,ma), colormap=args.colormap)
	
	# Create the camera object
	#cam = Camera(c, L, W, (args.res,args.res), fields=[fieldname], log_fields=[True],
	#	             transfer_function=tf,pf=pf ,north_vector=up, no_ghost=True)
	cam = Camera(c, L, W, (args.res,args.res), fields=[fieldname], log_fields=[args.log],
		transfer_function=tf, pf=pf, no_ghost=args.no_ghost, north_vector=up, sub_samples=args.sub_samples)
	
	# Take a "snapshot".  This is where the rendering actually occurs.  
	if args.fileprefix == '':
		fn = f
	else :
		fn = args.fileprefix+'_'+f

	image = cam.snapshot(fn+'_'+args.var+'.png')
	
	if args.rotstep != 0.:
		for i, snapshot in enumerate(cam.rotation(2.*na.pi, int(360./args.rotstep), args.rotdir)):
			if cam.comm.rank == 0:
				write_bitmap(snapshot, fn+'_'+args.var+('_rot_%04i' % i) + '.png')

