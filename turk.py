from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from yt.mods import *
mylog.setLevel(40)
 
pf = load("fiducial_1to3_b0.273d_hdf5_plt_cnt_0175")
sphere = pf.h.sphere(pf.domain_center, 0.25*pf.domain_width.max())
mi, ma = sphere["temp"].min(), sphere["temp"].max()
surface = pf.h.surface(sphere, "temp", (mi + ma) / 2.0)
colors = apply_colormap(np.log10(surface["dens"]), cmap_name="hot")
