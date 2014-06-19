import argparse

parser = argparse.ArgumentParser(description='Generate a projection of FLASH data.')
parser.add_argument('filename',                                help='filename', metavar='filename', nargs='+')
parser.add_argument('--var', '-v',         dest='var',         help='variable', default='Density', type=str)
parser.add_argument('--min',               dest='min',         help='minimum value to plot', type=float)
parser.add_argument('--max',               dest='max',         help='maximum value to plot', type=float)
parser.add_argument('--width',             dest='width',       help='size of cube to project through', type=float)
parser.add_argument('--center',            dest='center',      help='center of cube', type=float, nargs=3)
args = parser.parse_args()

from yt.mods import *
from scipy.ndimage.filters import gaussian_filter
import h5py

fieldname = args.var

def _ecoodens(field, data):
	return (-data["ecoo"]*data["dens"])

fieldname = fieldname.strip()

add_field(fieldname, function=_ecoodens)

for f in args.filename:
	pf = load(f)
	#arr = off_axis_projection(pf, normal_vector=[0.297421, 0.863774, -0.406737], field=fieldname, center=args.center, width=args.width, north_vector=[-0.181724, -0.367011, -0.912293], resolution=2000)
	arr = off_axis_projection(pf, normal_vector=[0.297421, 0.863774, -0.406737], field=fieldname, center=args.center, width=[args.width, args.width, args.width/10.],
							  north_vector=[-0.181724, -0.367011, -0.912293], resolution=2000)

	if (args.min != None) :
		mi = args.min
	else:
		mi = 'min'
	if (args.max != None) :
		ma = args.max
	else:
		ma = 'max'

	#p.set_zlim('all', mi, ma)
	#p.set_cmap('all', 'Greys_r')
	outname = args.var + '_projection_' + f
	fh = h5py.File(outname, 'w')
	#np.savetxt(fh, np.array([args.width]))
	header = np.concatenate((np.array([args.width, 2000]), args.center))
	fh.create_dataset("header", data=header)
	#fh.create_dataset(args.var, data=arr, compression='gzip')
	fh.create_dataset(args.var, data=arr)
	fh.close()
	#p.save()
