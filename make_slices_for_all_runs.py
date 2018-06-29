"""
makes all the slices for a bunch of runs at once, and ideally movies too,
so that I can just rsync it all at once to my local machine
"""
import os

cluster = '/groups/dark/lawsmith/'
#cluster = '/pfs/lawsmith/'

dirs = [
    #'m1.0_p10_b1.5',
    #'m1.0_p10_b2.5',
    #'m1.0_p16_b1.5',
    #'m1.0_p16_b2.0',
    #'m1.0_p16_b2.5',
    #'m1.0_p16_b3.0',
    #'m1.0_p16_b3.5',
    #'m1.0_p16_b4.0',
    'new_m1.0_p16_b4.0',
    #'43_b2.0_24k',
    #'43_b2.0_48k',
    #'43_b3.0_24k',
    #'43_b3.0_48k',
]

var = 'dens'
width = 1000

for dir in dirs:
    SAVE_PATH = cluster + 'fend_slices/' + dir + '/' + var + '_' + str(width) + 'rsun/'
    if not os.path.exists(SAVE_PATH): os.makedirs(SAVE_PATH)
    os.system('mpirun -np 8 python ~/jYT/plot_movie.py ' + dir + ' ' + var + ' ' + str(width) + ' ' + SAVE_PATH)
    movie_name = SAVE_PATH + dir.replace(".","_")
    slices_names = SAVE_PATH + 'multitidal_hdf5_chk_%04d_Slice_z_dens.png'
    os.system('ffmpeg -y -r 24 -i ' + slices_names + ' -b:v 4M ' \
        + movie_name + '.mp4')
    #  ' >> outputfile &'
