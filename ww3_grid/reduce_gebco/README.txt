USING THE gebco_reduce LIBRARY

PURPOSE

The purpose of this library is to generate a reduced version of the GEBCO
bathymetry product, based on a simple interpolation over an integer number
of cells. The reduction makes subsequent work with gridgen quicker and easier,
since the cell finding loops in that code then get to work with smaller arrays.

Running gebco_reduce also allows the user to correct sea levels/remove large
inland lakes and areas of land below mean sea level, and provides a percentage
land value for each of the coarsened cells.


RUNNING THE CODE

The library can be directly imported into a python environment, or run via the
'run_gebco_reduce.py' script using a given set of actions:
   - reduce: creates a reduced GEBCO bathymetry in a new netCDF file
   - correct: modifies a bathymetry file to correct/remove lakes
   - plot: visualize depth and/or land-sea mask information from a modified bathymetry file

Default options for the reduce/correct functions are encoded, however it is recommended
to set these via a configuration file (default name './gebco_reduce.cfg'). From a raw
copy of this library a new file can be set up by copying the file 

'gebco_reduce_defaults.cfg'

to a new file 

USERPATH/NEWNAME.cfg

The options can then be run as:

python run_gebco_reduce.py reduce USERPATH/NEWNAME.cfg

or
 
python run_gebco_reduce.py correct USERPATH/NEWNAME.cfg

or
 
python run_gebco_reduce.py plot USERPATH/NEWNAME.cfg
