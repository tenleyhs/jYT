'''
some plotting settings, and some constants
'''

import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (6,5)
mpl.rcParams['legend.fontsize'] = 16
mpl.rcParams['font.size'] = 20
#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['font.family'] = 'stixgeneral'

mpl.rcParams['legend.frameon'] = False
mpl.rcParams['legend.borderpad'] = 0.1
mpl.rcParams['legend.labelspacing'] = 0.1
mpl.rcParams['legend.handletextpad'] = 0.1

mpl.rcParams['xtick.top'] = True
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['xtick.major.size'] = 15
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.minor.width'] = 2
mpl.rcParams['xtick.minor.visible'] = True

mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['ytick.major.size'] = 15
mpl.rcParams['ytick.minor.size'] = 5
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.width'] = 2
mpl.rcParams['ytick.minor.visible'] = True

mpl.rcParams['axes.linewidth'] = 1

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
