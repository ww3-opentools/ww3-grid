
How to install
--------------

Library install
~~~~~~~~~~~~~~~~

Download this repository::

      git clone git@gitlab.com:oceanum/pysmc.git


The package should install with a standard pip install::

    pip install --user -r requirements.txt
    pip install --user .

A user install is recommended so as not to interfere with the system python
install, but things can of course be installed however you want. If you plan to
edit the code, replace the second line with a developer install e.g.::

    pip install --user -e .


This install compiles the required C and Fortran modules, and as you will need
to have C and Fortran compilers installed. The library also requires a netCDF
library and gtk library. For a fresh install of Ubuntu 19.04, the following
should fulfil these requirements::

   apt-get -y install python3-pip python3-cartopy gcc g++ gfortran \ 
              automake pkg-config build-essential libgtk2.0-dev \
              libnetcdf-dev 

Note that python3-cartopy is installed here through the system package manager
rather than relying in the pip install through the requirements file. This is
done as it provides an easy way to make sure the system requirements are
installed, but the install doesn't have to be carried out in this way. 

Command line tools
~~~~~~~~~~~~~~~~~~~~

There are two small helper scripts in the `bin` directory. To make it easier to
run the examples, put these in your path, e.g.:: 

   export PATH=$PATH:`pwd`/bin

This adds the following tools to your path for generating bathymetry and smc
grids from config files::

   bathyfromnc <bathy_config.yaml>
   smcgridgen <smc_config.yaml>

Next steps
~~~~~~~~~~~

Create a bathmetry file (e.g. `here <bathy_generation.html#Bathymetry Generation>`__), 
create an smc config file, and then use the smcgridgen command to build the grid. 
Grid options are shown in the config file below, serval examples are shown 
`here <examples.html#Examples>`__ 


.. literalinclude:: ./smc-config-full.yaml
   :language: yaml


Development notes
~~~~~~~~~~~~~~~~~~

The above installation compiles some C and Fortran code. 

Generate python modules from Fortran code (requires f2py)::

    cd SMCPY/fortran
    ./generate_fpys.sh

TODO - Clean up this info on the dev environment here
