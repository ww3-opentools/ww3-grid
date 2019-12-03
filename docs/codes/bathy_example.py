import xarray
import matplotlib.pyplot as plt
from SMCPy.bathy.bathy_interp import bathy_interp

filename_in = "../examples/data/GEBCO_2014_2D.nc"
filename_out = "/tmp/NZ.nc"

res = bathy_interp(filename_in, filename_out, 'elevation', 'lon', 'lat', 155, 178, 0.5, -60, -30, 0.5)
dset = xr.open_dataset('/tmp/output.nc')
dset.z.plot()
plt.show()
