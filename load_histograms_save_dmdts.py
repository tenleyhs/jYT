"""
read in all the histograms and save dmdts
optionally make plots. plot raw with alpha to compare. go through these plots to
decide what window lens, tmin, tmax to choose
comment out some ds if just ran a few new simulations and don't want to do everything

run like:
python load_histograms_save_dmdts.py --cluster=fend
"""
## TODO need to shift dmdes


import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.io import ascii
from scipy.ndimage.filters import gaussian_filter
from scipy.integrate import simps
from scipy.integrate import cumtrapz
from scipy import interpolate
from scipy import stats
import os
import argparse, sys

parser = argparse.ArgumentParser()
parser.add_argument('--cluster', help='cluster, e.g., hyades/fend/pfe', type=str, default='fend')
args = parser.parse_args()

if args.cluster == 'hyades': execfile('/pfs/lawsmith/jYT/my_settings.py')
elif args.cluster == 'fend': execfile('/groups/dark/lawsmith/jYT/my_settings.py')
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['font.size'] = 18

M_bh = 1e6*M_sun

# if choosing what smoothing to use
#LOOP_THRU_SIGMAS = True
#sigmas = [1,5,10,15,20,25,30,35,40,45,50]

PLOT_DMDES = True

NBINS = 800
#SMOOTHS = [10,50,75,100,110,120,130,140,150,160,170,180,190,200]
SMOOTHS = [100]

"""
ds = [
    # e1 and e2 are for the simple minimum-shifting-to-0 of dmde
    # todo not sure if to use M*, R* from multitidal.log or what, or if this makes a difference at all
    # run,                  chk,        sigma,  tmin,   tmax,   slope,  split,  sigma2, e1,     e2,     M*,     R*,     beta
    # d[0],                 d[1],       d[2],   d[3],   d[4],   d[5],   d[6],   d[7],   d[8],   d[9],   d[10],  d[11],  d[12]
    ['m0.3_p1_b0.6_300k',   '0100',     20,     -2,     2,      0,      0,      0,      -1,     1,      0.3,    0.2814, '0.600'],
    ['m0.3_p1_b0.7_300k',   '0100',     20,     -2,     2,      0,      0,      0,      -1,     1,      0.3,    0.2814, '0.700'],
    ['m0.3_p1_b0.8_300k',   '0100',     20,     -2,     2,      0,      0,      0,      -1,     1,      0.3,    0.2814, '0.800'],
    ['m0.3_p1_b0.9_300k',   '0100',     20,     -2,     2,      0,      0,      0,      -1,     1,      0.3,    0.2814, '0.900'],
    ['m0.3_p1_b1.0_300k',   '0100',     20,     -2,     2,      0,      0,      0,      -1,     1,      0.3,    0.2814, '1.000'],

    ['m0.7_p50_b0.8_300k',  '0100',     20,     -2,     2,      0,      0,      0,      -1,     1,      0.7,    0.6793, '0.800'],
    ['m0.7_p50_b1.0_300k',  '0100',     20,     -2,     2,      0,      0,      0,      -1,     1,      0.7,    0.6793, '1.000'],
    ['m0.7_p50_b1.15_300k',  '0100',     20,     -2,     2,      0,      0,      0,      -1,     1,      0.7,    0.6793, '1.150'],
    ['m0.7_p50_b1.3_300k',  '0100',     20,     -2,     2,      0,      0,      0,      -1,     1,      0.7,    0.6793, '1.300'],
    ['m0.7_p50_b1.5_300k',  '0100',     20,     -2,     2,      0,      0,      0,      -1,     1,      0.7,    0.6793, '1.500'],
    ['m0.7_p50_b1.7_300k',  '0104',     20,     -2,     2,      0,      0,      0,      -1,     1,      0.7,    0.6793, '1.700'],

    ['m1.0_p10_b1.0_300k',  '0100',     20,     -2,     2,      0,      0,      0,      -1,     1,      1.0,    1.0455, '1.000'],
    ['m1.0_p10_b1.5_300k',  '0090',     20,     -2,     2,      0,      0,      0,      -1,     1,      1.0,    1.0455, '1.500'],
    ['m1.0_p10_b2.0_300k',  '0075',     20,     -2,     2,      0,      0,      0,      -1,     1,      1.0,    1.0455, '2.000'],

    ['m1.0_p16_b1.5_300k',  '0090',     20,     -2,     2,      0,      0,      0,      -1,     1,      1.0,    1.2872, '1.500'],
    ['m1.0_p16_b2.0_300k',  '0075',     20,     -2,     2,      0,      0,      0,      -1,     1,      1.0,    1.2872, '2.000'],
    ['m1.0_p16_b3.0_300k',  '0060',     20,     -2,     2,      0,      0,      0,      -1,     1,      1.0,    1.2872, '3.000'],
    ['m1.0_p16_b4.0_300k',  '0050',     20,     -2,     2,      0,      0,      0,      -1,     1,      1.0,    1.2872, '4.000'],

    ['m3.0_p16_b1.0_300k',  '0080',     20,     -2,     2,      0,      0,      0,      -1,     1,      3.0,    3.3192, '1.000'],
    ['m3.0_p16_b2.0_300k',  '0075',     20,     -2,     2,      0,      0,      0,      -1,     1,      3.0,    3.3192, '2.000'],
    ['m3.0_p16_b3.0_300k',  '0060',     20,     -2,     2,      0,      0,      0,      -1,     1,      3.0,    3.3192, '3.000'],
    ['m3.0_p16_b4.0_300k',  '0050',     20,     -2,     2,      0,      0,      0,      -1,     1,      3.0,    3.3192, '4.000'],
    ['m3.0_p16_b4.5_300k',  '0051',     20,     -2,     2,      0,      0,      0,      -1,     1,      3.0,    3.3192, '4.500'],
    ['m3.0_p16_b4.5_300k_lr15',  '0052',     20,     -2,     2,      0,      0,      0,      -1,     1,      3.0,    3.3192, '4.500'],
]
"""
ds = [
    # e1 and e2 are for the simple minimum-shifting-to-0 of dmde
    # todo not sure if to use M*, R* from multitidal.log or what, or if this makes a difference at all
    # run,                  chk,        sigma,  tmin,   tmax,   slope,  split,  sigma2, e1,     e2,     M*,     R*,     beta        fpp
    # d[0],                 d[1],       d[2],   d[3],   d[4],   d[5],   d[6],   d[7],   d[8],   d[9],   d[10],  d[11],  d[12],      d[13]
    #['m3.0_p16_b1.0_300k',  '0080',     20,     -2,     2,      0,      0,      0,      -1,     1,      3.0,    3.3192, '1.000',    ''],
    #['m3.0_p16_b1.0_300k',  '0080',     20,     -2,     2,      0,      0,      0,      -1,     1,      3.0,    3.3192, '1.000',    'yt26_'],
    #['m3.0_p16_b1.0_300k',  '0080',     20,     -2,     2,      0,      0,      0,      -1,     1,      3.0,    3.3192, '1.000',    'yt26_1em11_'],
    ['m3.0_p16_b1.0_300k',  '0080',     20,     -2,     2,      0,      0,      0,      -1,     1,      3.0,    3.3192, '1.000',    'yt26_1em7_'],
    ['m3.0_p16_b1.0_300k',  '0080',     20,     -2,     2,      0,      0,      0,      -1,     1,      3.0,    3.3192, '1.000',    'yt26_1em6_'],
    ['m3.0_p16_b1.0_300k',  '0080',     20,     -2,     2,      0,      0,      0,      -1,     1,      3.0,    3.3192, '1.000',    'yt26_1em5_'],
]


#els = ['ev', 'h1', 'he4', 'o16', 'c12', 'ne20', 'n14']
#els = ['ev', 'h1', 'he3', 'he4', 'c12', 'c13', 'n14', 'o16', 'ne20', 'na23', 'mg24', 'al27', 'si28', 's34']
els = ['ev']

#ax directory = '/pfs/lawsmith/FLASH4.3/runs/results/dmdes/gaussian_wrap/' + d[0].replace(".","_")
#ax2 directory = '/pfs/lawsmith/FLASH4.3/runs/results/dmdts/dmde_first/' + d[0].replace(".","_")

# commenting above out because of some weird behavior where it wouldn't loop
# over els; would only do 6

# test make binning array for (log t, log mdot)
###a = np.array([0])

# plot final versions
#else:
for SMOOTH in SMOOTHS:
    for d in ds:
        for i, el in enumerate(els):
            fig2, ax2 = plt.subplots()
            fig3, ax3 = plt.subplots()
            fig4, ax4 = plt.subplots()
            fig5, ax5 = plt.subplots()
            if args.cluster == 'hyades': directory = '/pfs/lawsmith/FLASH4.3/runs/'
            elif args.cluster == 'fend': directory = '/storage/dark/lawsmith/FLASH4.3/runs/'
            e, dm = np.loadtxt(directory+d[0]+'/' + d[13] + 'b10000_'+el+'_bhbound_histogram_multitidal_hdf5_chk_'+d[1]+'.dat', skiprows=4)
            de = e[1]-e[0]
            
            # prune where dm=0, for b-spline fitting later
            e = e[dm!=0.0]
            dm = dm[dm!=0.0]

            dm_de = dm/de
            log_dm_de = np.log10(dm_de)
            # smooth linear value of dmde first
            s_dm_de = gaussian_filter(dm_de, d[2], mode='wrap')
            slog_dm_de = np.log10(s_dm_de)

            # todo should put the binning and b-spline up here, then center

            # shift dmde histogram to be centered at 0
            if d[0] == 'm1.0_p1_b1.0':
                e = e - e[(-1<e/1e17) & (e/1e17<1)][np.argmax(slog_dm_de[(-1<e/1e17) & (e/1e17<1)])]
            if d[0] == 'm1.0_p1_b2.0':
                e = e + 0.1e17
            if d[0] == 'm1.0_p1_b3.0':
                e = e - e[(-1<e/1e17) & (e/1e17<1)][np.argmin(slog_dm_de[(-1<e/1e17) & (e/1e17<1)])]
            if d[0] == 'm1.0_p10_b1.0':
                e = e - e[(-1<e/1e17) & (e/1e17<1)][np.argmin(slog_dm_de[(-1<e/1e17) & (e/1e17<1)])]
            if d[0] == 'm1.0_p10_b2.0':
                e = e - e[np.argmax(slog_dm_de)]
            if d[0] == 'm1.0_p10_b3.0':
                pass
            if d[0] == 'm1.0_p16_b1.0':
                e = e - e[(-1<e/1e17) & (e/1e17<1)][np.argmin(slog_dm_de[(-1<e/1e17) & (e/1e17<1)])]
            if d[0] == 'm1.0_p16_b2.0':
                e = e - e[(-1<e/1e17) & (e/1e17<1)][np.argmax(slog_dm_de[(-1<e/1e17) & (e/1e17<1)])]
            if d[0] == 'm1.0_p16_b3.0':
                e = e - e[(-1<e/1e17) & (e/1e17<1)][np.argmax(slog_dm_de[(-1<e/1e17) & (e/1e17<1)])]
            if d[0] == 'm1.0_p16_b3.0_96k_1000':
                e = e - e[(0<e/1e17) & (e/1e17<1)][np.argmax(slog_dm_de[(0<e/1e17) & (e/1e17<1)])]

            if d[0] == 'm1.0_p16_b4.0':
                e = e - e[(0<e/1e17) & (e/1e17<0.3)][np.argmin(slog_dm_de[(0<e/1e17) & (e/1e17<0.3)])]
            if d[0] == 'm3.0_p1_b2.0_48k':
                e = e - e[(-1<e/1e17) & (e/1e17<1)][np.argmin(slog_dm_de[(-1<e/1e17) & (e/1e17<1)])]
            if d[0] == 'm3.0_p16_b4.0_48k':
                e = e - e[np.argmax(slog_dm_de)]

            # todo fast shift dmde for m0.3, all minimums. this isn't correct for the strange double-minimums for some
            e = e - e[(d[8]<e/1e17) & (e/1e17<d[9])][np.argmin(slog_dm_de[(d[8]<e/1e17) & (e/1e17<d[9])])]

            # try different binning rather than smoothing, trying to copy Gafton & Rosswog 2019
            dm_de_binned = stats.binned_statistic(e, dm_de, statistic='mean', bins=NBINS)
            dm_de_binned_std = stats.binned_statistic(e, dm_de, statistic='std', bins=NBINS)

            # try different binning rather than smoothing, trying to copy Gafton & Rosswog 2019
            bin_width = (dm_de_binned[1][1] - dm_de_binned[1][0])
            bin_centers = dm_de_binned[1][1:] - bin_width/2

            # try without weights..
            sp1 = interpolate.splrep(bin_centers, dm_de_binned[0], w=1./dm_de_binned_std[0], s=SMOOTH)
            #sp1 = interpolate.splrep(bin_centers, dm_de_binned[0], s=SMOOTH) <-- doesn't seem to work
            #If the weights represent the inverse of the standard-deviation of y, then a good s value should be found 
            #in the range (m-sqrt(2*m),m+sqrt(2*m)) where m is the number of datapoints in x, y, and w. 
            #default : s=m-sqrt(2*m) if weights are supplied. s = 0.0 (interpolating) if no weights are supplied.
            # ^ this gives range of 760 to 840
            x2 = np.linspace(min(bin_centers), max(bin_centers), NBINS)
            y2 = interpolate.splev(x2, sp1)

            # todo move / make general. fast overplot of GRR2013 5/3.
            #g13_53 = np.loadtxt('/groups/dark/lawsmith/Guillochon2013_dmdts/5-3/' + d[12] + '.dat')
            # adjust to M~0.3 msun, R~0.28 rsun
            # M* = d[10],     R* = d[11]
            #ax2.plot(g13_53[:,0] - np.log10(d[10]) + 1.5*np.log10(d[11]), g13_53[:,1] + 2.0*np.log10(d[10]) - 1.5*np.log10(d[11]), c='g', label='GRR13')


            e_bound = e[np.where(e<0.)]
            t = 2.*np.pi*G*M_bh/((2*np.abs(e_bound))**(3./2))
            de_dt = (1./3)*((2.*np.pi*G*M_bh)**(2./3))*t**(-5./3)
            dm_de_bound = 10**slog_dm_de[np.where(e<0.)]
            dm_de_bound_orig = dm_de[np.where(e<0.)]
            mdot = dm_de_bound*de_dt
            mdot_orig = dm_de_bound_orig*de_dt

            log_t_yr = np.log10(t/yr)
            log_mdot_moyr = np.log10(mdot*yr/M_sun)
            log_mdot_moyr_orig = np.log10(mdot_orig*yr/M_sun)

            # plot unsmoothed
            ax2.scatter(log_t_yr, log_mdot_moyr_orig, c='C0', rasterized=True, alpha=0.5, s=1, edgecolors='none')

            # old gaussian mdots
            sel = np.where((log_t_yr>d[3]) & (log_t_yr<d[4]))[0]
            x = log_t_yr[sel]
            y = log_mdot_moyr[sel]
            ax2.plot(x, y, c='r', label='gaussian_filter')
            ax3.plot(x, (y-np.min(y[np.isfinite(y)]))/(np.max(y)-np.min(y)), c='C1')

            # plot B-spline mdots
            x2_bound = x2[np.where(x2<0.)]
            t2 = 2.*np.pi*G*M_bh/((2*np.abs(x2_bound))**(3./2))
            de_dt2 = (1./3)*((2.*np.pi*G*M_bh)**(2./3))*t2**(-5./3)
            dm_de_bound2 = y2[np.where(x2<0.)]
            mdot2 = dm_de_bound2*de_dt2
            ax2.plot(np.log10(t2/yr), np.log10(mdot2*yr/M_sun), c='k', label='B-spline')

            # todo commented out below 2/15/2020
            """
            # hacky thing to fix extended slope for a few dmdts
            # TODO change this
            if d[6] != 0:
                slog_dm_de1 = gaussian_filter(log_dm_de, d[2], mode='wrap')
                dm_de_bound1 = 10**slog_dm_de1[np.where(e<0)]
                mdot1 = dm_de_bound1*de_dt
                log_mdot_moyr1 = np.log10(mdot1*yr/M_sun)

                slog_dm_de2 = gaussian_filter(log_dm_de, d[7], mode='wrap')
                dm_de_bound2 = 10**slog_dm_de2[np.where(e<0)]
                mdot2 = dm_de_bound2*de_dt
                log_mdot_moyr2 = np.log10(mdot2*yr/M_sun)

                # join them
                sel1 = np.where(log_t_yr < d[6])[0]
                sel2 = np.where(log_t_yr > d[6])[0]
                log_mdot_moyr = np.append(log_mdot_moyr1[sel1], log_mdot_moyr2[sel2])
            """

            #linear time
            #ax5.plot(10**(g13_53[:,0] - np.log10(d[10]) + 1.5*np.log10(d[11]))*yr/day, g13_53[:,1] + 2.0*np.log10(d[10]) - 1.5*np.log10(d[11]), c='g', label='GRR13')
            ax5.scatter(t/day, log_mdot_moyr_orig, c='C0', rasterized=True, alpha=0.5, s=1, edgecolors='none')
            ax5.plot((10**x)*yr/day, y, c='r', label='gaussian_filter')
            ax5.plot(t2/day, np.log10(mdot2*yr/M_sun), c='k', label='B-spline')

            # try gaussian filtering mdot instead of dmde..idk. the number of points in t bins dereases as some power law, so want something, idk
            # gives same
            #s_mdot = gaussian_filter(mdot_orig, d[2], mode='wrap')
            #ax2.plot(np.log10(t/yr), np.log10(s_mdot*yr/M_sun), c='C1', label='gaussian_filter of mdot')

            # todo try binning with increasing bin widths
            ###mdot_binned = stats.binned_statistic(np.log10(t/yr), mdot_orig, statistic='mean', bins=)

            # todo commented out below 2/15/2020
            """
    		# extend dmdt from slope near end. could maybe improve slope function
            # this should just apply to the else statement above. could combine?
            if d[4] < 2.0:
                ninf = ((y[-10] - y[-1]) / (x[-10] - x[-1])) + d[5]
                yext = np.linspace(y[-1], y[-1] + (2.0 - x[-1])*ninf, num=100)[1:]
                xext = np.linspace(x[-1], 2.0, num=100)[1:]
                ax2.plot(xext, yext, c='C2')
                ax3.plot(xext, (yext-np.min(yext))/(np.max(yext)-np.min(yext)), c='C2')
                x = np.append(x, xext)
                y = np.append(y, yext)
                ax3.plot(x, (y-np.min(y[np.isfinite(y)]))/(np.max(y)-np.min(y)), c='C3')
            """

            # interpolate for smoother integrating/interpolating
            # doesn't change integration answer though.
            f = interpolate.interp1d(x, y)
            xnew = np.linspace(min(x), max(x), num=10000)
            ynew = f(xnew)
            x = xnew
            y = ynew

            # write for later plotting
            if args.cluster == 'hyades': directory = '/pfs/lawsmith/results/dmdts/data/'
            elif args.cluster == 'fend': directory = '/groups/dark/lawsmith/results/'+d[0].replace(".","_")+'/'
            if not os.path.exists(directory): os.makedirs(directory)
            ascii.write([x, y], directory + d[13] + d[0].replace(".","_")+'_'+d[1]+'_'+el+'.dat', overwrite=True, names=['log_t_yr','log_mdot_moyr'])

            """
            # integrate mdot curves, to compare with delta m, and save for later
            if el == 'ev':
                int_trapz = np.trapz(10**y, 10**x)
                int_simps = simps(10**y, 10**x)
                if args.cluster == 'hyades':
                    ascii.write([[int_trapz], [int_simps]],
                        '/pfs/lawsmith/results/dmdts/integrals/'+d[0].replace(".","_")+
                        '_'+d[1]+'_'+el+'_'+str(d[2])+'.dat', overwrite=True, names=['int_trapz','int_simps'])
                elif args.cluster == 'fend':
                    ascii.write([[int_trapz], [int_simps]],
                        '/groups/dark/lawsmith/results/'+d[0].replace(".","_")+'/'+
                        'ints_'+d[0].replace(".","_")+
                        '_'+d[1]+'_'+el+'_'+str(d[2])+'.dat', overwrite=True, names=['int_trapz','int_simps'])

                # cumulative plots
                int_cumtrapz = cumtrapz(10**y, 10**x, initial=0)
                ax3.plot(x, int_cumtrapz)
                ax3.set_xlim(-2, 2)
                ax3.set_ylim(0, 1)
                #ax3.set_ylabel(r'$\log\ \dot M\ {\rm [M_\odot/yr]}$')
                ax3.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
                fig3.tight_layout(pad=0.3)
                if args.cluster == 'hyades':
                    directory = '/pfs/lawsmith/results/dmdts/integrals/'
                elif args.cluster == 'fend':
                    directory = '/groups/dark/lawsmith/results/'+d[0].replace(".","_")+'/'+'ints_'
                if not os.path.exists(directory): os.makedirs(directory)
                fig3.savefig(directory \
                	+ d[0].replace(".","_") + '_' + d[1] + '_' + el + '_' + str(d[2]) + '.pdf')
            """

            if PLOT_DMDES:
                fig, ax = plt.subplots()
                #try normalizing to \Delta e, as in Ryu+2020b, Rees, Phinney
                #delta_e = (G * M_star / R_star) * (M_bh / M_star)**(1./3)
                #M*, R*
                #d[10], d[11]
                delta_e =  G * (d[10]*M_sun)**(2./3) * M_bh**(1./3) / (d[11]*R_sun)

                #ax.scatter(e/1e17, log_dm_de, c='C0', rasterized=True, alpha=0.5, s=1, edgecolors='none')
                #ax.plot(e/1e17, slog_dm_de, c='C1')
                ax.scatter(e/delta_e, np.log10(dm_de * delta_e/(d[10]*M_sun)), c='C0', rasterized=True, alpha=0.5, s=1, edgecolors='none')
                ax.plot(e/delta_e, np.log10(s_dm_de * delta_e/(d[10]*M_sun)), c='r', label='gaussian_filter')
                
                # or could do by [10^-4 c^2] like Gafton & Rosswog 2019

                #ax.hlines(np.log10(dm_de_binned[0] * delta_e/(d[10]*M_sun)), dm_de_binned[1][:-1]/delta_e, dm_de_binned[1][1:]/delta_e, colors='k', label='binned_statistic')
                ax.scatter(bin_centers/delta_e, np.log10(dm_de_binned[0] * delta_e/(d[10]*M_sun)), c='r', s=5, label='binned_statistic', edgecolors='none', rasterized=True)

                # b spline from binned dm_de ^
                ax.plot(x2/delta_e, np.log10(y2 * delta_e/(d[10]*M_sun)), c='k', label='B-spline')

                #ax.set_xlabel(r'$E\ \mathrm{[10^{17}\ erg\ g^{-1}]}$')
                #ax.set_ylabel(r'$\log\ dM/dE\ \mathrm{[g^2\ erg^{-1}]}$')
                ax.set_xlabel(r'$e/\Delta e$')
                ax.set_ylabel(r'$\log\ dM/de\ [M_\star / \Delta e]$')

                ax.legend()
                ax.set_title(d[0])

                #ax.set_xlim(-10, 10)
                #ax.set_ylim(6, 16)

                ax.grid(which='both')
                ax.grid(which='minor', lw=0.1, ls=':')
                fig.tight_layout(pad=0.3)
                if args.cluster == 'hyades': 
                    directory = '/pfs/lawsmith/results/dmdes/gaussian_wrap/' + d[0].replace(".","_")
                elif args.cluster == 'fend': 
                    #todo temp
                    directory = '/groups/dark/lawsmith/results/'+d[0].replace(".","_")
                    #directory = '/groups/dark/lawsmith/results/temp'
                if not os.path.exists(directory): os.makedirs(directory)
                #fig.savefig(directory+'/dmde_'+d[0].replace(".","_")+'_'+d[1]+'_'+el+'_'+str(d[2])+'.pdf')
                fig.savefig(directory+'/dmde_'+d[13] +d[0].replace(".","_")+'_'+d[1]+'_'+el+'_'+str(SMOOTH)+'.pdf')

                # make histogram of e, to try out other smoothing algorithms
                """
                ax4.hist(np.log10(t/yr), 100)
                ax4.set_xlabel(r'$\log t\ \mathrm{[yr]}$')
                ax4.set_ylabel('N')
                ax4.grid()
                fig4.tight_layout(pad=0.3)
                fig4.savefig(directory+'/hist_t_'+d[0].replace(".","_")+'_'+d[1]+'_'+el+'_'+str(d[2])+'.pdf')            
                """

            #ax2.set_xlim(-2, 2)
            #ax2.set_ylim(-8, 1)
            ax2.set_xlim(right=2)
            ax2.set_ylim(bottom=-8)

            ax2.legend()
            ax2.set_title(d[0])

            ax2.set_ylabel(r'$\log\ \dot M\ {\rm [M_\odot/yr]}$')
            ax2.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
            ax2.grid(which='both')
            ax2.grid(which='minor', lw=0.1, ls=':')
            fig2.tight_layout(pad=0.3)
            if args.cluster == 'hyades': 
                directory = '/pfs/lawsmith/results/dmdts/final'
            elif args.cluster == 'fend': 
                #todo temp
                directory = '/groups/dark/lawsmith/results/'+d[0].replace(".","_")
                #directory = '/groups/dark/lawsmith/results/temp'
            #fig2.savefig(directory+'/dmdt_'+d[0].replace(".","_")+'_'+d[1]+'_'+el+'_'+str(d[2])+'.pdf')
            fig2.savefig(directory+'/dmdt_'+d[13] +d[0].replace(".","_")+'_'+d[1]+'_'+el+'_'+str(SMOOTH)+'.pdf')

            ax5.set_xlim(0, 100)
            ax5.set_ylim(bottom=-8)
            ax5.legend()
            ax5.set_title(d[0])
            ax5.set_ylabel(r'$\log\ \dot M\ {\rm [M_\odot/yr]}$')
            ax5.set_xlabel(r'$t\ \mathrm{[day]}$')
            ax5.grid(which='both')
            ax5.grid(which='minor', lw=0.1, ls=':')
            fig5.tight_layout(pad=0.3)
            #fig5.savefig(directory+'/dmdt_linear_'+d[0].replace(".","_")+'_'+d[1]+'_'+el+'_'+str(d[2])+'.pdf')
            fig5.savefig(directory+'/dmdt_linear_'+d[13] +d[0].replace(".","_")+'_'+d[1]+'_'+el+'_'+str(SMOOTH)+'.pdf')

            plt.close('all')
