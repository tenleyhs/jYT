from yt.mods import * # set up our namespace
import matplotlib.colorbar as cb
from matplotlib.colors import LogNorm
from matplotlib import figure
import numpy as na
from yt.funcs import *
from matplotlib.backends.backend_agg import FigureCanvasAgg

filenames = ["multibipoly_hdf5_plt_cnt_00038",
             "multibipoly_hdf5_plt_cnt_00053",
             "multibipoly_hdf5_plt_cnt_00068",
             "multibipoly_hdf5_plt_cnt_00083",
             "multibipoly_hdf5_plt_cnt_00113",
             "multibipoly_hdf5_plt_cnt_00172"]

timelabels = [r'$t=-1 t_{\rm dyn}$',
              r'$t=0$',
              r'$t= 1 t_{\rm dyn}$',
              r'$t= 2 t_{\rm dyn}$',
              r'$t= 4 t_{\rm dyn}$',
              r'$t= 8 t_{\rm dyn}$']


# set up our Fixed Resolution Buffer parameters: a width, resolution, and center
width = [4*8.75e11,60*8.75e11] 
res = [500,500]

# set up to do one row of length the same as # of filenames... 
orient = 'vertical'
#fig, axes, colorbars = get_multi_plot( 6, 2, colorbar=orient, bw = 3,dpi=300)
nx = 6
ny = 2
bw = 3
dpi = 300
hf,wf = 1.0/ny, 1.0/nx
cbar_padding = 0.4
fudge_x = nx/(cbar_padding+nx)
fudge_y = 1.0
print "nx,ny=",nx,ny
print "dpi=",dpi
print "hf,wf=",hf,wf
fig = figure.Figure((bw*nx/fudge_x, bw*ny/fudge_y), dpi=dpi)

fig.set_canvas(FigureCanvasAgg(fig))
fig.subplots_adjust(wspace=0.0, hspace=0.0,
                    top=1.0, bottom=0.0,
                    left=0.0, right=1.0)
tr = []
print "fudge_x, fudge_y=",fudge_x, fudge_y
for j in range(ny):
    tr.append([])
    for i in range(nx):
        left = i*wf*fudge_x
        bottom = fudge_y*(1.0-(j+1)*hf) + (1.0-fudge_y)
        ax = fig.add_axes([left, bottom, wf*fudge_x, hf*fudge_y])
        tr[-1].append(ax)
cbars = []


#for j in range(ny):
#    ax = fig.add_axes([wf*(nx+0.05)*fudge_x, hf*fudge_y*(ny-(j+0.95)),
#                       wf*fudge_x*0.05, hf*fudge_y*0.90])
ax = fig.add_axes([wf*(nx+0.05)*fudge_x, 0.075,
                   wf*fudge_x*0.05, 0.85])
ax.clear()
cbars.append(ax)

fig = fig
axes = tr
colorbars = cbars


plots = []

for i,fn in enumerate(filenames):
    # load data
    pf = load(fn) 
    #c = pf.h.find_max('Density')[1]
    c = [5.e14,5.e14,5.e14]
    
    # make a slice + frb
    for j in range(2):
        sli = pf.h.slice(2,c[2]) #z=const slice
        frb = sli.to_frb(width[j], res, center=c)
        thisframe = axes[j][i]
        # get rid of axes
        thisframe.xaxis.set_visible(False)
        thisframe.yaxis.set_visible(False)
        # annotate?
        # with an arrow
        #thisframe.annotate('text',xy=(250,250),xycoords='data',xytext=(-40,30),textcoords='offset points',arrowprops=dict(arrowstyle="->"))

        #time labels on lower panel
        if j == 1:
            thisframe.annotate(timelabels[i],xy=(0.1,0.05),xycoords='axes fraction',color='white',size=20)

        # add the plots to the figure
        plots.append(thisframe.imshow(frb['Density'], norm=LogNorm(),origin='lower'))
        plots[-1].set_clim((1.e-8, 10))
        if j==0:
            plots[-1].set_cmap("afmhot")
        elif j==1:
            plots[-1].set_cmap("afmhot")
            



    titles = [r'$\mathrm{Density}\ (\mathrm{g\ cm^{-3}})$']
    for p,cb,t in zip(plots,colorbars,titles):
        cbar = fig.colorbar(p,cax=cb,orientation=orient)
        cbar.set_label(t)
        

# And now we're done!  
fig.savefig("timeseries.pdf",format='pdf')
