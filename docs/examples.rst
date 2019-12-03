Examples
===========

New Zealand
------------

All configuration for this example can be found in `examples/nz`. In this
directory are three files


NetCDF Bathmetry
~~~~~~~~~~~~~~~~~~

.. Produce bathymetry for plots

.. plot:: 

   import os
   from SMCPy.utils import load_object, load_config 
   path = '../examples/nz'
   os.chdir(path)
   os.system('rm -f bathy_nz.nc')
   config = load_config('nc_bath.yml')
   load_object(config)

First, lets create the bathymetry::

   cd examples/glob1
   bathyfromnc nc_bath.yml

This config is shown below:

.. literalinclude:: ../examples/nz/nc_bath.yml
   :language: yaml

Now let's run the smc grid generation with the following configuration:

.. literalinclude:: ../examples/nz/nz_smc_nc.yml
   :language: yaml

::

   smcgridgen nz_smc_nc.yml

This will produce all the required inputs for WAVEWATCHIII in the `output`
directory. It will also produce the followings plots of cells and filled cells. 


.. plot:: 
   :include-source:

   import os
   from SMCPy.utils import load_object, load_config 
   path = '../examples/nz'
   os.chdir(path)
   config = load_config('nz_smc_nc.yml')
   smc = load_object(config)
   smc.run()
   plt.show()


Gridgen Bathmetry
~~~~~~~~~~~~~~~~~~

This example is also replicated here using the second method of creating
bathymetry. The matlab file nz.mat is already included here, but can be
recreated by running the following script from within the nz examples
directory::

   python run_gridgen.py

TODO - This needs to be documented better

Then, we can generate SMC grids in a similar way pointing to the mat config.
This is exactly the same as the one above, except its points to a different
bathymetry source, and uses the GRIDGEN class rather than the NC class. 

.. literalinclude:: ../examples/nz/nz_smc_mat.yml
   :language: yaml

.. plot:: 

   import os
   from SMCPy.utils import load_object, load_config 
   path = '../examples/nz'
   os.chdir(path)
   config = load_config('nz_smc_mat.yml')
   smc = load_object(config)
   smc.run()
   plt.show()


Global 1 degree
------------------

.. Produce bathymetry for plots

.. plot:: 

   import os
   from SMCPy.utils import load_object, load_config 
   path = '../examples/glob1'
   os.chdir(path)
   os.system('rm -f glob1.nc')
   config = load_config('nc_bathy.yml')
   load_object(config)

Global examples from `examples/glob1`


First, lets create the bathymetry

.. literalinclude:: ../examples/glob1/nc_bathy.yml
   :language: yaml


::

   cd examples/glob1
   bathyfromnc nc_bathy.yml


3-tier land
~~~~~~~~~~~~

A config to create a 3 tier global SMC grid would look something like this:


.. literalinclude:: ../examples/glob1/glob1-3_smc_land.yaml
   :language: yaml


Running that, produces the following::

   smcgridgen glob1-3_smc.yaml

.. plot:: 

   import os
   from SMCPy.utils import load_object, load_config 
   path = '../examples/glob1'
   os.chdir(path)
   config = load_config('glob1-3_smc_land.yaml')
   smc = load_object(config)
   smc.run()
   plt.show()


3-tier depth
~~~~~~~~~~~~

Refining on depth, rather than just proximity to land, here we refine on all
point in a depth of less than 150 m. 


.. literalinclude:: ../examples/glob1/glob1-3_smc.yaml
   :language: yaml

.. plot:: 

   import os
   from SMCPy.utils import load_object, load_config 
   path = '../examples/glob1'
   os.chdir(path)
   config = load_config('glob1-3_smc.yaml')
   smc = load_object(config)
   smc.run()
   plt.show()


High latitude merging
~~~~~~~~~~~~~~~~~~~~~~~

A single merge at 60 degrees

.. literalinclude:: ../examples/glob1/antactic_merge-1.yaml
   :language: yaml

.. plot:: 

   import os
   from SMCPy.utils import load_object, load_config 
   path = '../examples/glob1'
   os.chdir(path)
   config = load_config('antactic_merge-1.yaml')
   smc = load_object(config)
   smc.run()
   plt.show()

Two merges at 55 and again at 65 degrees

.. literalinclude:: ../examples/glob1/antactic_merge-2.yaml
   :language: yaml

.. plot:: 

   import os
   from SMCPy.utils import load_object, load_config 
   path = '../examples/glob1'
   os.chdir(path)
   config = load_config('antactic_merge-1.yaml')
   smc = load_object(config)
   smc.run()
   plt.show()


Exclusions
~~~~~~~~~~~~~~~~~~~~~~~~~

Points can be completely removed from consideration. At the moment, this is
done based on a number of pre-defined polygons for enclosed bays and seas that
are either turned or off based on a user-defined flag file. This is similar to
the approach taken by the NOAA gridgen software, and indeed regions and flag
file formats are identical. The example shown below is the same as the 3-teir
global example above, but with all of these flags activated. 

.. literalinclude:: ../examples/glob1/user_polygons.flag
   :language: text

.. literalinclude:: ../examples/glob1/glob1-3_smc_exclude.yaml
   :language: yaml

.. plot:: 

   import os
   from SMCPy.utils import load_object, load_config 
   path = '../examples/glob1'
   os.chdir(path)
   config = load_config('glob1-3_smc_exclude.yaml')
   smc = load_object(config)
   smc.run()
   plt.show()


Manual Regional Refining
~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, some manual regional refining is demonstrated. Specifically,
we refine over a large rectangular area over the the U.K. and turn off all
automatic refineing over North America. This example also included at high
latitude merge at 60 deg.

.. literalinclude:: ../examples/glob1/atlantic_manual.yaml
   :language: yaml

.. plot:: 

   import os
   from SMCPy.utils import load_object, load_config 
   path = '../examples/glob1'
   os.chdir(path)
   config = load_config('atlantic_manual.yaml')
   smc = load_object(config)
   smc.run()
   plt.show()



