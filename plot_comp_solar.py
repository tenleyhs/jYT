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

f_h1_ZAMS = 0.7153
f_he4_ZAMS = 0.2704
f_c12_ZAMS = 2.435e-3
f_n14_ZAMS = 7.509e-4
f_o16_ZAMS = 6.072e-3
f_ne20_ZAMS = 1.232e-3

do_smoothing = True
plot_ttpeak = True

ds = [ 
    ['m1.0_p1_b1.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0100.dat', 300, 0., 'p1_b1.0', 'age=0Gyr, '+r'$\beta=1.0$'],
    ['m1.0_p1_b2.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0080.dat', 250, -0.15, 'p1_b2.0', 'age=0Gyr, '+r'$\beta=2.0$'],
    ['m1.0_p1_b3.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0060.dat', 250, 0.25, 'p1_b3.0', 'age=0Gyr, '+r'$\beta=3.0$'],
    ['m1.0_p10_b1.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0100.dat', 500, 0.2, 'p10_b1.0', 'age=4.8Gyr, '+r'$\beta=1.0$'],
    ['m1.0_p10_b1.0_256/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0100.dat', 500, 0.45, 'p10_b1.0_256', 'age=4.8Gyr, '+r'$\beta=1.0$'+' (256)'],
    ['m1.0_p10_b1.5/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0090.dat', 200, -0.2, 'p10_b1.5', 'age=4.8Gyr, '+r'$\beta=1.5$'],
    ['m1.0_p10_b2.0_128/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0075.dat', 200, 0.15, 'p10_b2.0', 'age=4.8Gyr, '+r'$\beta=2.0$'],
    ['m1.0_p10_b2.0_256/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0080.dat', 200, 0.8, 'p10_b2.0_256', 'age=4.8Gyr, '+r'$\beta=2.0$'],
    ['m1.0_p10_b3.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0060.dat', 200, 0.45, 'p10_b3.0', 'age=4.8Gyr, '+r'$\beta=3.0$'],
    ['m1.0_p10_b4.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0050.dat', 200, 0.3, 'p10_b4.0', 'age=4.8Gyr, '+r'$\beta=4.0$'],
    ['m1.0_p10_b5.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0055.dat', 200, 0.35, 'p10_b5.0', 'age=4.8Gyr, '+r'$\beta=5.0$'],
    ['m1.0_p16_b2.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0075.dat', 300, 0.6, 'p16_b2.0', 'age=8.4Gyr, '+r'$\beta=2.0$'],
    ['m1.0_p16_b3.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0060.dat', 500, 0., 'p16_b3.0', 'age=8.4Gyr, '+r'$\beta=3.0$'],
    ['m1.0_p16_b4.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0050.dat', 500, -0.25, 'p16_b4.0', 'age=8.4Gyr, '+r'$\beta=4.0$'],
    ['m1.0_p16_b5.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0050.dat', 500, -0.25, 'p16_b5.0', 'age=8.4Gyr, '+r'$\beta=5.0$']
]  

#els = ['ev','h1', 'he4', 'c12', 'n14', 'o16', 'ne20']
# monica order
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
    fig, ax = plt.subplots()
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

        if do_smoothing:
            # smooth hanning dmde
            slog_dm_de = smooth(log_dm_de, window_len=d[2])
            slog_mdot_moyr = smooth(log_mdot_moyr, window_len=d[2])
            sel = np.where((log_t_yr<d[3]) & (np.isfinite(slog_mdot_moyr)))[0]
            tpeak = 10**log_t_yr[sel][np.argmax(slog_mdot_moyr[sel])]
            if el == 'h1':
                h1_mdot = 10**slog_mdot_moyr[sel]
                continue
            elif el == 'o16':
                el_lr = (10**slog_mdot_moyr[sel]/h1_mdot)/(f_o16_ZAMS/f_h1_ZAMS)
            elif el == 'he4':
                el_lr = (10**slog_mdot_moyr[sel]/h1_mdot)/(f_he4_ZAMS/f_h1_ZAMS)
            elif el == 'c12':
                el_lr = (10**slog_mdot_moyr[sel]/h1_mdot)/(f_c12_ZAMS/f_h1_ZAMS)
            elif el == 'n14':
                el_lr = (10**slog_mdot_moyr[sel]/h1_mdot)/(f_n14_ZAMS/f_h1_ZAMS)

            if plot_ttpeak:
                ax2.plot(np.log10(10**log_t_yr[sel]/tpeak), el_lr, lw=lw, label=el)
            else:
                ax2.plot(log_t_yr[sel], el_lr, lw=lw, label=el)

        else:
            ax2.plot(log_t_yr, log_mdot_moyr, lw=lw, rasterized=True, label=el)

    
    if plot_ttpeak:
        ax2.set_xlim(-0.5, 1.0)
        ax2.text(-0.48,2.75,d[5])
        ax2.set_xlabel(r'$\log\ t/t_{\rm peak}$')
    else:
        ax2.set_xlim(-1.5, 0.0)
        ax2.text(-1.45,2.8,d[5])
        #ax2.text(-1.95,3.75,d[5])
        ax2.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
    ax2.set_ylim(0, 3)
    ax2.axhline(1, ls=':', lw=1.5, c='k', alpha=0.5)
    ax2.set_ylabel(r'$X/X_\odot$')
    ax2.legend(loc=1)
    fig2.tight_layout()
    if plot_ttpeak:
        fig2.savefig('/pfs/lawsmith/FLASH4.3/runs/results/mdot_comp_solar/mdot_comp_solar_'+d[4]+'_ttpeak.png')
    else:
        fig2.savefig('/pfs/lawsmith/FLASH4.3/runs/results/mdot_comp_solar/mdot_comp_solar_'+d[4]+'.png')
