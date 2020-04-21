README FOR gridgen AND genSides FUNCTIONS

#------------------------
USING THE gridgen LIBRARY

PURPOSE

The purpose of this library is to generate a WAVEWATCH III format regular and 
SMC grid files for use by ww3_grid


RUNNING THE CODE

The library can be directly imported into a python environment, or run via the
'run_gridgen.py' script using a given set of actions. These need to be associated
with a user defined category in the namelist supplied to 'run_gridgen.py', since
the action could be re-used multiple times(see worked examples):


For SMC models:

   - smcbase: creates an SMC base grid in a new netCDF file
   - tiergen: creates an SMC tier in a new netCDF file 
   - tiercombine: combines SMC base/parent and SMC tier data in a new netCDF file
   - writeWW3: writes data out in WAVEWATCH III format for ww3_grid plus smcGrid.nml
               file for cell face generation and propagation testing

Default options for the reduce/correct functions are encoded, however it will usually
be necessary to set these via a configuration file (default name './gebco_reduce.cfg').
From a raw copy of this library a new file can be set up by copying the file 

'gridgen_defaults.cfg'

to a new file 

USERPATH/NEWNAME.cfg


#-------------------------
USING THE genSides LIBRARY

PURPOSE

The purpose of this library is to generate a the face cells arrays needed by
WAVEWATCH III format SMC grid files for use in ww3_grid


RUNNING THE CODE

Before running the code, 'genSides.f90' needs to be built using the
make_genSides.sh (or other local) bash script

The genSides functions ('genSides.exe' and 'countijsdnew') are then
run using the bash script 'run_genSides.sh' as

./run_genSides.sh WORKINGDIRECTORYPATH smcGrid.nml

where WORKINGDIRECTORYPATH is the path to a directory containing the 'ww3Cels.dat'
('ww3ArcCels.dat') and 'smcGrid.nml' ('arcGrid.nml') files generated using 'gridgen.py'.
The .nml namelists provide the necessary control information for genSides and can
be adapted if needed (e.g. different input Cels.dat filename).


