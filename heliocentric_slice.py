# Description: this program reads the orbit.dat to find the position of the point mass and make slices centered on the point mass. The point mass will be represented as a filled circle in the slices.
# Written by Shangfei Liu for visualizaing disrupted giant planets.
from yt.mods import *

# Make sure the following directory exists!
outputdir = 'slice/'

begin = 0
end = 55
step = 1

for i in range(begin, end+1, step):
	if i < 10:
		pf = load('multibipoly_hdf5_plt_cnt_0000' + str(i))
	elif i < 100:
		pf = load('multibipoly_hdf5_plt_cnt_000' + str(i))
	else:
		pf = load('multibipoly_hdf5_plt_cnt_00' + str(i))

	time = pf.current_time
#	print time

	index = 0

# Could be more efficient if it only needs to read orbit.dat once.
	for line in open('orbit.dat'):
		if index > 0 :
			parts = line.split()
			if float(parts[0]) >= time :
				"""
				if (float(parts[0]) - time) > (time - float(lastparts[0])) :
					print parts[0], lastparts[0]
				index -= 1
				"""
				"""
				x1=parts[7]
				x2=parts[1]
				y1=parts[8]
				y2=parts[2]
				"""
				print parts[0]
				break
		index += 1

	pc = PlotCollection(pf, [5.25e12-float(parts[7])+float(parts[1]),5.75e12-float(parts[8])+float(parts[2]),5e12])


	p=pc.add_slice("Density",2)
#	pc.add_projection("Density",2)

	pc.set_width(0.195, 'au')
	pc.set_zlim(1e-12,1e1)
#	pc.set_zlim(1e-5,1e9)
	pc.set_cmap('gist_heat')

	p.modify["sphere"]((5e12-float(parts[7])+float(parts[1]), 5e12-float(parts[8])+float(parts[2])), 6.955e10, {'ec': '#e67817', 'fc': '#f8c301','fill': True})
#	p.modify["marker"]([5.25e12-float(parts[7])+float(parts[1]), 5.75e12-float(parts[8])+float(parts[2])])
#	p.modify["contour"]("Density")

	pc.save(outputdir+"%s" % pf)
