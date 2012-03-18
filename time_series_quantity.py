"""
TRY TO LOOP OVER FILES AND PRINT ENCLOSED MASS
"""
import asciitable
from yt.mods import * # set up our namespace

# First set up our times and quantities lists
data = {
    'a_fn':[],
    'b_t':[],
    'c_menc0':[],
    'c_menc1':[],
    'c_menc2':[],
    'c_menc3':[],
    'c_menc4':[],
    'c_menc5':[],
    'c_menc6':[],
    'c_menc7':[],
    'c_menc8':[],
    'c_menc9':[]
    }

# Read filenames
filenames = glob.glob('multibipoly_hdf5_plt_cnt_00[0-2][0-9][0-9]')
filenames = sorted(filenames)

#define sampling radii
radii = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]


mass =  [0,0,0,0,0,0,0,0,0,0] #just to initialize array
#loop over files
for i,filename in enumerate(filenames):
    print "FILE {}, {}".format(i,filename)
    pf = load(filename)
    v, c = pf.h.find_max("Density")
    for j,rad in enumerate(radii):
        rad = rad*8.75e11
        sp = pf.h.sphere(c, rad)
        #print "j {} rad {}".format(j,rad)
        mlocal = sp.quantities["TotalQuantity"](
            ["CellMassMsun"], lazy_reader=True)
        mass[j] = mlocal[0]
        #print "mass {} mlocal {}".format(mass[j],mlocal[0])
        

    # now construct the output table
    data['a_fn'].append(filename)
    data['b_t'].append(pf.current_time)
    data['c_menc0'].append(mass[0])
    data['c_menc1'].append(mass[1])
    data['c_menc2'].append(mass[2])
    data['c_menc3'].append(mass[3])
    data['c_menc4'].append(mass[4])
    data['c_menc5'].append(mass[5])
    data['c_menc6'].append(mass[6])
    data['c_menc7'].append(mass[7])
    data['c_menc8'].append(mass[8])
    data['c_menc9'].append(mass[9])

    # print to screen
    print "Time {}, mass within radius {} cm is {} Msun".format(data['b_t'][i],radii[3],data['c_menc3'][i])

# write the data table
asciitable.write(data,"menc.dat")




