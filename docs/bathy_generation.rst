Bathymetry Generation
======================

SMC grids require input bathymetry. The library requires a bathymetry file at the same resolution as the highest tear SMC grid. For example, if you are creating a global grid with a base resolution of 1 degree and three SMC tier levels, the required input resolution is 0.25 degree. 


There are two methods included in this library for making bathymetry, the first is using the provided grid interpolation routines included in this library, the second is using NOAA's gridgen software. Each has it advantages and disadvantages, these are discussed below. 



Included interpolation library
--------------------------------


These approximate the techniques used in the NOAA gridgen software, but are based on raw input bathymetry only, and do not include any accounting of coastal polygons. The underlying code is written in C, so its fast and requires no additional software installs. For course grids, this is likely adequate. For finer grids, and where calculation of subgridscale obstacles are required, consider using the second method below. 

This method is relatively straight forward. A demonstration of the generation of a 0.5 degree grid over Australia is given below: 

.. plot:: 
   :include-source:

   import xarray as xr
   import matplotlib.pyplot as plt
   from SMCPy.bathy.bathy_interp import bathy_interp
   import tempfile

   filename_in = "../examples/data/etopo2.nc"
   filename_out = tempfile.mktemp()

   res = bathy_interp(filename_in, filename_out, 'z', 'x', 'y', 100, 165, 0.5, -50, -5, 0.5)
   dset = xr.open_dataset(filename_out)
   dset.z.plot()
   plt.show()


In the provided examples, configuration files are provided to easily construct this file, using the CLI tool. For example, a full configuration file for New Zealand is shown below:

.. literalinclude:: ../examples/nz/nc_bath.yml
   :language: yaml

For more examples, see the examples `here <examples.html#Examples>`__ 

Gridgen
--------

This method uses the standard NOAA gridgen library. This code is written in MATLAB. It is wrapped here in python, and the user need not be familiar with MATLAB in order to use this approach, however, you will need to have MATLAB installed, and in your path. You will also need to have gridgen installed, and you may need change the paths in SMCPy/matlab/create_grid_smcbase.m to point to the gridgen routines in your particular installation.

TODO - need more documentation here

An example using gridgen is given `here <examples.html#New Zealand>`__ 
