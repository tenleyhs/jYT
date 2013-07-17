import argparse

parser = argparse.ArgumentParser(description='Generate a volumetric plot of FLASH data.')
parser.add_argument('filename',                                help='filename', metavar='filename', nargs='+')
parser.add_argument('--parallel',          dest='parallel',    help='run in parallel', default=False, action='store_true')
parser.add_argument('--var', '-v',         dest='var',         help='variable', default='Density', type=str)
parser.add_argument('--log', '-l',         dest='log',         help='take the log of the variable', action='store_true', default=False)
parser.add_argument('--res', '-r',         dest='res',         help='image resolution', default=512, type=int)
parser.add_argument('--min',               dest='min',         help='minimum value to plot', type=float)
parser.add_argument('--max',               dest='max',         help='maximum value to plot', type=float)
parser.add_argument('--center', '-c',      dest='center',      help='which point to center on', choices=('min', 'max', 'geometric', 'custom', 'pttrack'), default='geometric')
parser.add_argument('--centervar', '-cv',  dest='centervar',   help='variable to use when centering', type=str, default='var')
parser.add_argument('--width', '-w',       dest='width',       help='width of region to plot', type=float)
parser.add_argument('--sub_samples', '-s', dest='sub_samples', help='sub_samples parameter for camera', type=int, default=5)
parser.add_argument('--colormap', '-cm',   dest='colormap',    help='colormap used for plot', type=str, default='jet')
parser.add_argument('--ncontours', '-n',   dest='ncontours',   help='number of contours to draw', type=int, default=6)
parser.add_argument('--negative',          dest='negative',    help='plots -1.0*value', action='store_true', default=False)
parser.add_argument('--axis', '-a',        dest='axis',        help='set sliceplane axis [x,y,z]', type=str, default='x')
args = parser.parse_args()

from yt.mods import *
from math import *
from yt.visualization.volume_rendering.camera import PerspectiveCamera

up = [0.,-1.,0.]
L = na.array([0.,1.,1.])

if (args.negative) :
	sign = -1.0
else :
	sign = 1.0

fieldname = args.var
if (args.var == 'kine') :
	def _var(field, data):
		return sign*0.5*data['dens']*(data['velx']**2 + data['vely']**2 + data['velz']**2)
elif (args.var == 'vel2') :
	def _var(field, data):
		return sign*(data['velx']**2 + data['vely']**2 + data['velz']**2)
elif (args.var == 'eratio') :
	def _var(field, data):
		return sign*(data['velx']**2 + data['vely']**2 + data['velz']**2)/data['eint']
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
	elif args.center == 'pttrack' :
		odata = na.loadtxt('pruned_orbit.dat', dtype='float64')
		ptvec = odata[:,1:7]
		obvec = odata[:,7:13]
		totvec = odata[:,19:25]
		ptvec += totvec - obvec
		time = odata[:,0]
		tindex = abs(time - pf.current_time).argmin()
		c = na.array([ptvec[tindex,0],ptvec[tindex,1]+8e12,ptvec[tindex,2]])

	sp = pf.h.sphere(c, W)
	mi, ma = sp.quantities['Extrema'](args.var)[0]
	print mi, ma

	if (args.min != None) :
		mi = args.min
	if (args.max != None) :
		ma = args.max

	sp = SlicePlot(pf, args.axis, args.var, center = c, width = W)
	sp.set_zlim('all', mi, ma)
	sp.save()
