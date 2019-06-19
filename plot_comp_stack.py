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
plot_ttpeak = False

ds = [
    ['m1.0_p1_b1.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0100.dat', 300, 0., 'p1_b1.0', 'age=0Gyr, '+r'$\beta=1.0$'],
    ['m1.0_p1_b2.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0080.dat', 250, -0.15, 'p1_b2.0', 'age=0Gyr, '+r'$\beta=2.0$'],
    ['m1.0_p1_b3.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0060.dat', 250, 0.25, 'p1_b3.0', 'age=0Gyr, '+r'$\beta=3.0$'],
    ['m1.0_p10_b1.0/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0100.dat', 500, 0.2, 'p10_b1.0', 'age=4.8Gyr, '+r'$\beta=1.0$'],
    ['m1.0_p10_b1.0_256/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0100.dat', 500, 0.45, 'p10_b1.0_256', 'age=4.8Gyr, '+r'$\beta=1.0$'+' (256)'],
    ['m1.0_p10_b1.5/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0090.dat', 200, -0.2, 'p10_b1.5', 'age=4.8Gyr, '+r'$\beta=1.5$'],
    ['m1.0_p10_b2.0_128/b10000_','_bhbound_histogram_multitidal_hdf5_chk_0075.dat', 200, 0.25, 'p10_b2.0', 'age=4.8Gyr, '+r'$\beta=2.0$'],
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
els = ['ev','h1', 'he4', 'o16', 'c12', 'ne20', 'n14']

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

    tpeak = 0.
    last_mdot = 0.
    cum_mdot = 0.
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
            #slog_dm_de = smooth(log_dm_de, window_len=d[2])
            slog_mdot_moyr = smooth(log_mdot_moyr, window_len=d[2])
            sel = np.where((log_t_yr<d[3]) & (np.isfinite(slog_mdot_moyr)))[0]
            if el == 'ev':
                total_mdot = 10**slog_mdot_moyr[sel]
                tpeak = 10**log_t_yr[sel][np.argmax(total_mdot)]
            else:
                if plot_ttpeak:
                    ax2.fill_between(np.log10(10**log_t_yr[sel]/tpeak), np.log10(1.-cum_mdot/total_mdot), 
                        np.log10(1.-(10**slog_mdot_moyr[sel]+cum_mdot)/total_mdot), label=el, alpha=0.75)
                else:
                    ax2.fill_between(log_t_yr[sel], np.log10(1.-cum_mdot/total_mdot), 
                        np.log10(1.-(10**slog_mdot_moyr[sel]+cum_mdot)/total_mdot), label=el, alpha=0.75)
                last_mdot = 10**slog_mdot_moyr[sel]
                if el == 'h1':
                    cum_mdot = last_mdot
                else:
                    cum_mdot += last_mdot

        else:
            ax2.plot(log_t_yr, log_mdot_moyr, lw=lw, rasterized=True, label=el)


    ax2.set_ylim(-2.5, 0)
    ax2.set_ylabel(r'$\log\ \dot M_x/\dot M_{\rm total}$')
    if plot_ttpeak:
        ax2.fill_between(np.log10(10**log_t_yr[sel]/tpeak), np.log10(1.-cum_mdot/total_mdot), 
            np.log10(0.00001), label='other', color='gray', alpha=0.5)
        ax2.set_xlabel(r'$\log\ t/t_{\rm peak}$')
        ax2.text(-0.95,-0.2, d[5])
        ax2.set_xlim(-1, 1.7)
    else:
        ax2.fill_between(log_t_yr[sel], np.log10(1.-cum_mdot/total_mdot), 
            np.log10(0.00001), label='other', color='gray', alpha=0.5)
        ax2.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
        #ax2.text(-1.95,-0.2, d[5])
        ax2.text(-1.45,-0.2, d[5])
        ax2.set_xlim(-1.5, 0)
    ax2.legend()
    fig2.tight_layout(pad=0.3)
    if plot_ttpeak:
        fig2.savefig('/pfs/lawsmith/FLASH4.3/runs/results/mdot_comp_stack/mdot_comp_stack_'+ d[4] +'_ttpeak.png')
    else:
        fig2.savefig('/pfs/lawsmith/FLASH4.3/runs/results/mdot_comp_stack/mdot_comp_stack_'+ d[4] +'.png')
