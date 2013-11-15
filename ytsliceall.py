from yt.mods import * # set up our namespace
import matplotlib.pyplot as plt

# pick appropriate filenames, and range that you want... eg:
#filenames = glob.glob('multibipoly_hdf5_chk_[0-9][0-9][0-9][0-9][0-9]')
filenames = glob.glob('multibipoly_hdf5_plt_cnt_00[2-3][0-9][0-9]')
filenames = sorted(filenames)


# loop through files
for i,fn in enumerate(filenames):
    pf = load(fn) # load data
    pc = PlotCollection(pf) # defaults to center at most dense point
    p = pc.add_slice("Density", 2,use_colorbar=False) # 2 = z-axis
    #p.modify["contour"]("Density",6,4,True,[1.e-8,1.e-3],{'colors':'cyan'})
    #p.modify["contour"]("Temperature",5,4,True,[1.e4,1.e9],{'colors':'red'})
    #p.modify["velocity"]()
    #p.modify["coord_axes"]()
    pc.set_width(4.e13, 'cm') # change width of all plots in cm
    pc.set_cmap("afmhot")
    pc.set_zlim(1.e-7,5)
    # save plot
    pc.save(fn) 
