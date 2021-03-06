import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import UnivariateSpline
#execfile('/nobackup/jlawsmit/jYT/my_settings.py')
execfile('/pfs/lawsmith/jYT/my_settings.py')
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['font.size'] = 18

# PARAMS
M_bh = 1e6*M_sun
lw = 1.5
do_smoothing = True

ds = [
	['m1.0_p1_b1.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0100.dat', 300, 0., 'p1_b1.0', 'age=0Gyr, '+r'$\beta=1.0$'],
	['m1.0_p1_b2.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0080.dat', 250, -0.15, 'p1_b2.0', 'age=0Gyr, '+r'$\beta=2.0$'],
	['m1.0_p1_b3.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0060.dat', 250, 0.25, 'p1_b3.0', 'age=0Gyr, '+r'$\beta=3.0$'],
	['m1.0_p10_b1.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0100.dat', 500, 0.2, 'p10_b1.0', 'age=4.8Gyr, '+r'$\beta=1.0$'],
	['m1.0_p10_b1.0_256/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0100.dat', 500, 0.45, 'p10_b1.0_256', 'age=4.8Gyr, '+r'$\beta=1.0$'+' (256)'],
	['m1.0_p10_b1.5/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0090.dat', 200, -0.2, 'p10_b1.5', 'age=4.8Gyr, '+r'$\beta=1.5$'],
	['m1.0_p10_b2.0_128/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0075.dat', 200, 0.25, 'p10_b2.0', 'age=4.8Gyr, '+r'$\beta=2.0$'],
	['m1.0_p10_b2.0_256/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0080.dat', 200, 0.8, 'p10_b2.0_256', 'age=4.8Gyr, '+r'$\beta=2.0$'+' (256)'],
	['m1.0_p10_b3.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0060.dat', 200, 0.45, 'p10_b3.0', 'age=4.8Gyr, '+r'$\beta=3.0$'],
	['m1.0_p10_b4.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0050.dat', 200, 0.3, 'p10_b4.0', 'age=4.8Gyr, '+r'$\beta=4.0$'],
	['m1.0_p10_b5.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0055.dat', 200, 0.35, 'p10_b5.0', 'age=4.8Gyr, '+r'$\beta=5.0$'],
	['m1.0_p16_b2.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0075.dat', 300, 0.6, 'p16_b2.0', 'age=8.4Gyr, '+r'$\beta=2.0$'],
	['m1.0_p16_b3.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0060.dat', 500, 0., 'p16_b3.0', 'age=8.4Gyr, '+r'$\beta=3.0$'],
	['m1.0_p16_b4.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0050.dat', 500, -0.25, 'p16_b4.0', 'age=8.4Gyr, '+r'$\beta=4.0$'],
	['m1.0_p16_b5.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0050.dat', 500, -0.25, 'p16_b5.0', 'age=8.4Gyr, '+r'$\beta=5.0$']
	]

els = ['h1', 'he4', 'o16', 'c12', 'ne20', 'n14']

# for dmde's
def smooth(x,window_len=11,window='hanning'):
	if x.ndim != 1:
		raise ValueError, "smooth only accepts 1 dimension arrays."
	if x.size < window_len:
		raise ValueError, "Input vector needs to be bigger than window size."
	if window_len<3:
		return x
	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
	s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
	if window == 'flat': #moving average
		w=np.ones(window_len,'d')
	else:
		w=eval('np.'+window+'(window_len)')
	y=np.convolve(w/w.sum(),s,mode='same')
	return y[window_len:-window_len+1]



for d in ds:
    #fig, ax = plt.subplots()
    fig2, ax2 = plt.subplots()

    for el in els:
        e, dm = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/'+d[0]+el+d[1], skiprows=4)
        de = e[1]-e[0]
        dm_de = dm/de

        # calculate mdot
        e_bound = e[np.where(e<0.)]
        t = 2.*np.pi*G*M_bh/((2*np.abs(e_bound))**(3./2))
        de_dt = (1./3)*((2*np.pi*G*M_bh)**(2./3))*t**(-5./3)
        dm_de_bound = dm_de[np.where(e<0)]
        mdot = dm_de_bound*de_dt

        # take logs. better to smooth what I'm actually plotting (log-linear or log-log)
        log_t = np.log10(t)
        log_t_yr = np.log10(t/yr)
        log_mdot = np.log10(mdot)
        log_mdot_moyr = np.log10(mdot*yr/M_sun)
        log_dm_de = np.log10(dm_de)

        # mdot eddington
        #eta = 0.1
        #mdot_edd = (2.2e-8)*(eta/0.1)*(M_bh/M_sun)*M_sun/yr
        #log_mdot_edd = np.log10(mdot_edd)
        
        if do_smoothing:
            # smooth hanning dmde
            #slog_dm_de = smooth(log_dm_de, window_len=d[2])
            slog_mdot_moyr = smooth(log_mdot_moyr, window_len=d[2])
            #ax.plot(e/1e17, slog_dm_de, lw=lw, label=el)
            sel = np.where(log_t_yr<d[3])[0]
            ax2.plot(log_t_yr[sel], slog_mdot_moyr[sel], lw=lw, label=el)

        else:
            #ax.plot(e/1e17, log_dm_de, lw=lw, rasterized=True, label=el)
            ax2.plot(log_t_yr, log_mdot_moyr, lw=lw, rasterized=True, label=el)


    ax2.set_xlim(-2.25, 1)
    ax2.set_ylim(-6, 1)
    ax2.set_ylabel(r'$\log\ \dot M\ {\rm [M_\odot/yr]}$')
    ax2.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
    ax2.legend()
    ax2.text(-2.2,0.5, d[5])
    fig2.tight_layout(pad=0.3)
    fig2.savefig('/pfs/lawsmith/FLASH4.3/runs/results/mdot_comp/mdot_comp_'+ d[4] +'.png')
    plt.close()

