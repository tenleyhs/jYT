'''
some plotting settings and colors, and some constants
'''

import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = (6,5)
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['legend.borderpad'] = 0.1
plt.rcParams['legend.labelspacing'] = 0.1
plt.rcParams['legend.handletextpad'] = 0.1
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'stixgeneral'
plt.rcParams['font.size'] = 16

plt.rcParams['xtick.major.size'] = 6.0
plt.rcParams['ytick.major.size'] = 6.0
plt.rcParams['xtick.minor.size'] = 2.0
plt.rcParams['ytick.minor.size'] = 2.0
plt.rcParams['xtick.major.width'] = 1.0
plt.rcParams['ytick.major.width'] = 1.0
plt.rcParams['xtick.minor.width'] = 0.6
plt.rcParams['ytick.minor.width'] = 0.6
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True

plt.rcParams['xtick.major.top'] = True
plt.rcParams['xtick.minor.top'] = True
plt.rcParams['ytick.major.right'] = True
plt.rcParams['ytick.minor.right'] = True


# http://oobrien.com/2012/01/tube-colours/
brown = '#B36305'
red = '#E32017'
yellow = '#FFD300'
green = '#00782A'
pink = '#F3A9BB'
gray = '#A0A5A9'
mauve = '#9B0056'
blue = '#003688'
brightblue = '#0098D4'
softblue = '#95CDBA'
otherblue = '#00A4A7'
orange = '#EE7C0E'
softgreen = '#84B817'
brightred = '#E21836'
purple = '#7156A5'

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
