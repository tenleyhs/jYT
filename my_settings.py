'''
some plotting settings, and some constants
'''

import matplotlib as mpl

mpl.rcParams['figure.figsize'] = (6,5)

mpl.rcParams['legend.frameon'] = False
mpl.rcParams['legend.borderpad'] = 0.1
mpl.rcParams['legend.labelspacing'] = 0.1
mpl.rcParams['legend.handletextpad'] = 0.1

mpl.rcParams['xtick.top'] = True
mpl.rcParams['xtick.major.size'] = 5 #3.5
mpl.rcParams['xtick.minor.size'] = 3 #2
mpl.rcParams['xtick.major.width'] = 1.0 #0.8
mpl.rcParams['xtick.minor.width'] = 1.0 #0.6
mpl.rcParams['xtick.minor.visible'] = True

mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.major.size'] = 5 #3.5
mpl.rcParams['ytick.minor.size'] = 3 #2
mpl.rcParams['ytick.major.width'] = 1.0 #0.8
mpl.rcParams['ytick.minor.width'] = 1.0 #0.6
mpl.rcParams['ytick.minor.visible'] = True

#mpl.rcParams['legend.fontsize'] = 14
#mpl.rcParams['font.size'] = 16
#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['font.family'] = 'stixgeneral'

#xtick.top            : False   # draw ticks on the top side
#xtick.bottom         : True   # draw ticks on the bottom side
#xtick.major.size     : 3.5      # major tick size in points
#xtick.minor.size     : 2      # minor tick size in points
#xtick.major.width    : 0.8    # major tick width in points
#xtick.minor.width    : 0.6    # minor tick width in points
#xtick.major.pad      : 3.5      # distance to major tick label in points
#xtick.minor.pad      : 3.4      # distance to the minor tick label in points
#xtick.color          : k      # color of the tick labels
#xtick.labelsize      : medium # fontsize of the tick labels
#xtick.direction      : out    # direction: in, out, or inout
#xtick.minor.visible  : False  # visibility of minor ticks on x-axis
#xtick.major.top      : True   # draw x axis top major ticks
#xtick.major.bottom   : True   # draw x axis bottom major ticks
#xtick.minor.top      : True   # draw x axis top minor ticks
#xtick.minor.bottom   : True   # draw x axis bottom minor ticks

#ytick.left           : True   # draw ticks on the left side
#ytick.right          : False  # draw ticks on the right side
#ytick.major.size     : 3.5      # major tick size in points
#ytick.minor.size     : 2      # minor tick size in points
#ytick.major.width    : 0.8    # major tick width in points
#ytick.minor.width    : 0.6    # minor tick width in points
#ytick.major.pad      : 3.5      # distance to major tick label in points
#ytick.minor.pad      : 3.4      # distance to the minor tick label in points
#ytick.color          : k      # color of the tick labels
#ytick.labelsize      : medium # fontsize of the tick labels
#ytick.direction      : out    # direction: in, out, or inout
#ytick.minor.visible  : False  # visibility of minor ticks on y-axis
#ytick.major.left     : True   # draw y axis left major ticks
#ytick.major.right    : True   # draw y axis right major ticks
#ytick.minor.left     : True   # draw y axis left minor ticks
#ytick.minor.right    : True   # draw y axis right minor ticks


# CONSTANTS [cgs]
# todo should just use astropy probably
day = 86400.
yr = 3.154e7
c = 2.993e10
G = 6.674e-8
pc = 3.086e18
AU =  1.496e13
M_sun = 1.988435e33
R_sun = 6.955e10
M_jup = 1.898e30
R_jup = 6.9173e9
L_sun = 3.828e33
