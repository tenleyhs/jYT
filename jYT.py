import numpy as na
import scipy.interpolate
import scipy.ndimage

def congrid(a, newdims, method='linear', centre=False, minusone=False, offset=0., buffer=0):
	'''Arbitrary resampling of source array to new dimension sizes.
	Currently only supports maintaining the same number of dimensions.
	To use 1-D arrays, first promote them to shape (x,1).

	Uses the same parameters and creates the same co-ordinate lookup points
	as IDL''s congrid routine, which apparently originally came from a VAX/VMS
	routine of the same name.

	method:
	neighbour - closest value from original data
	nearest and linear - uses n x 1-D interpolations using
						 scipy.interpolate.interp1d
	(see Numerical Recipes for validity of use of n 1-D interpolations)
	spline - uses ndimage.map_coordinates

	centre:
	True - interpolation points are at the centres of the bins
	False - points are at the front edge of the bin

	minusone:
	For example- inarray.shape = (i,j) & new dimensions = (x,y)
	False - inarray is resampled by factors of (i/x) * (j/y)
	True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
	This prevents extrapolation one element beyond bounds of input array.
	'''
	if not a.dtype in [na.float64, na.float32]:
		a = na.cast[float](a)

	m1 = na.cast[int](minusone)
	ofs = na.cast[int](centre) * 0.5
	old = na.array( a.shape )
	ndims = len( a.shape )
	if len( newdims ) != ndims:
		print "[congrid] dimensions error. " \
			  "This routine currently only support " \
			  "rebinning to the same number of dimensions."
		return None
	newdims = na.asarray( newdims, dtype=float )
	dimlist = []

	if method == 'neighbour':
		for i in range( ndims ):
			base = na.indices(newdims)[i]
			dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
							* (base + ofs) - ofs - offset)
		cd = na.array( dimlist ).round().astype(int)
		newa = a[list( cd )]
		return newa

	elif method in ['nearest','linear']:
		# specify old dims
		olddims = [na.arange(i, dtype = na.float) for i in list( a.shape )]
		# calculate new dims
		for i in range( ndims ):
			base = na.arange( buffer+1, newdims[i]-buffer )
			newval = (old[i] - m1) / (newdims[i] - m1) \
							* (base + ofs) - ofs - offset
			newval = [t for t in newval if t >= olddims[-1][0] and t <= olddims[-1][-1]]
			dimlist.append( newval )

		# first interpolation - for ndims = any
		mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
		newa = mint( dimlist[-1] )

		trorder = [ndims - 1] + range( ndims - 1 )
		for i in range( ndims - 2, -1, -1 ):
			newa = newa.transpose( trorder )

			mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
			newa = mint( dimlist[i] )

		if ndims > 1:
			# need one more transpose to return to original dimensions
			newa = newa.transpose( trorder )

		return newa
	elif method in ['spline']:
		oslices = [ slice(0,j) for j in old ]
		oldcoords = na.ogrid[oslices]
		nslices = [ slice(0,j) for j in list(newdims) ]
		newcoords = na.mgrid[nslices]

		newcoords_dims = range(na.rank(newcoords))
		#make first index last
		newcoords_dims.append(newcoords_dims.pop(0))
		newcoords_tr = newcoords.transpose(newcoords_dims)
		# makes a view that affects newcoords

		newcoords_tr += ofs

		deltas = (na.asarray(old) - m1) / (newdims - m1)
		newcoords_tr *= deltas

		newcoords_tr -= ofs

		newa = scipy.ndimage.map_coordinates(a, newcoords)
		return newa
	else:
		print "Congrid error: Unrecognized interpolation type.\n", \
			  "Currently only \'neighbour\', \'nearest\',\'linear\',", \
			  "and \'spline\' are supported."
		return None
