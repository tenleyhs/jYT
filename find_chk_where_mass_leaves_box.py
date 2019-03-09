import yt

for i in range(101):


    ds = yt.load('multitidal_hdf5_chk_' + str(i).zfill(4))
    sp = ds.sphere("center", (1.0, "Mpc"))
    mass = sp.quantities.total_quantity(["cell_mass"]).in_units('Msun')
    print i, round(mass, 5)
