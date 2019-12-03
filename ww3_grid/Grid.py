#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import netCDF4 as nc
import xarray
import os
from tqdm import tqdm
import logging
from scipy import ndimage
import yaml
import pandas as pd
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep
from matplotlib.path import Path  # Path is much faster than Shapely
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader

from SMCPy.plotting import UnSMC, CartopyMap
from SMCPy.fortran.GenCellSide import gencellside
from SMCPy import Grid

logging.basicConfig(level=logging.INFO)

class NC2SMC(object):

    """
    This class will calculate Sperical-Multi-Cell (SMC) grids from regular
    netCDF file based on user defined criteria. It will perform all the
    required steps to produce the input files for running SMC grids in
    WAVEWATCH III

        - Calculate SMC grid cell based on user defined criteria
        - Calculate face arrays

    Note that the resoluiton of the input netcdf grid is the lowest tier or the
    SMC grid

    Args:
        bathymetry (str):   path to netcdf bathymetry.
        workdir    (str):   output directory.
        gid        (str):   grid id.
        grid       (dict):
            - Rotated       (bool):    Rotated grid (True, False).
            - urlat         (float):   Upper right latitude (otherwise, take from file).
            - lllat         (float):   Lower left latitude (otherwise, take from file).
            - urlon         (float):   Upper right longitude (otherwise, take from file).
            - lllon         (float):   Lower left longitude (otherwise, take from file)
            - mergelat      (float):   First merge latitude
            - mergelat2     (float):   Second merge latitude
            - zoom          (int):     Zoom grid
            - polygon_check (bool):    Exclude based on coastal polygon test (True, False)
            - exclude_flags (bool):    path to exclude flags
            - refine_region (dict):    force refining of a region defined by the supplied polygon
                - all: # e.g. Refine over oceania
                    - lats: [-50, -50, 20, 20, 20]
                    - lons: [100, 190, 190, 100, 100]
                - none: # e.g. ignore the Americas
                    - lats: [-50, -50, 20, 20, 20]
                    - lons: [250, 330, 330, 250, 250]
        conventions (dict):
            - yname         (str):     y variable name
            - xname         (str):     x variable name
            - zname         (str):     z variable name
            - zscale        (str):     z scale
            - xyorder       (bool):    xyorder (True, False)
        smc (dict):
            - smctiers      (int):     Number of smc tiers
            - depthlim      (int):     Depth to refine on
            - depthvar      (int):     Depth variation to refine on
        ww3meta (dict):
            - latlonscale   (int):     latlonscale
            - llcrnrscale   (int):     latlonscale
            - lsmdepth:     (int):     TODO
            - mindepth:     (int):     Model minimum depth
            - blockscale    (int):     Scale factor for blocking information - TODO - not sure what this is doing
        plotting (dict):
            - projection (str):        Projection, default PlateCarree
            - clon (float):            Central latitude, default 180.0
            - clat (float):            Central longitude, default 0.0
            - fillplot (bool):         Plot with filled cells, default True
            - cellplot (bool):         Plot with unfilled cells, default True
            - interactive (bool):      True
        actions (dict):
            - cell          (bool):    True
            - face          (bool):    True
            - plot          (bool):    True

    Returns:

        None

    Outputs:

        The following files are produced

          - <outputdir>/<gid>BPlist.txt   - boundary point list
          - <outputdir>/<gid>Cell.dat     - SMC Cell file
          - <outputdir>/<gid>_grid.inp    - WW3 grid input dimensions
          - <outputdir>/<gid>ISide.dat    - east-west face arrays
          - <outputdir>/<gid>JSide.dat    - north-south face arrays
          - <outputdir>/<gid>_META.txt    - grid metadata
          - <outputdir>/<gid>Obs.dat      - obstruction grid
          - <outputdir>/<gid>SMCSide.txt - grid statistics


    """

    def __init__(
        self,
        bathymetry,
        gid="grid",
        workdir="./",
        grid={},
        conventions={},
        smc={},
        ww3meta={},
        actions={},
        plotting={},
        logger=logging,
    ):


        self.workdir = workdir
        self.bathymetry = bathymetry
        self.gid = gid
        self.grid = grid
        self.conventions = {"xyorder": False}
        self.conventions.update(conventions)
        self.smc = smc
        self.ww3meta = ww3meta
        self.logger = logger
        # set input rotation and grid values
        self.rotated = False
        if "Rotated" in self.grid:
            self.rotated = self.grid.get("Rotated", False)
            if self.rotated:
                self.rlat = np.float(self.grid["rlat"])
                self.rlon = np.float(self.grid["rlon"])
                self.ingrid_lllon = np.float(self.grid["lllon"])
                self.ingrid_lllat = np.float(self.grid["lllat"])
                self.dx = np.float(self.grid["dx"])
                self.dy = np.float(self.grid["dy"])

        self.logger.info("==============================================")
        self.logger.info("Creating SMC Grid")
        self.logger.info("\t input bathymetry %s" % bathymetry)
        self.logger.info("\t output bathymetry %s" % workdir)
        self.logger.info("Grid ---> {0}".format(grid))
        self.logger.info("Conventions ---> {0}".format(conventions))
        self.logger.info("SMC ---> {0}".format(smc))
        self.logger.info("==============================================")

        # grid constraints for limited areas
        if "llx" in self.grid and "lllon" in self.grid:
            raise Exception("Can only respect llx or lllon, remove one from config")
        if "llx" in self.grid:
            self.llx = np.int(self.grid["llx"])
        else:
            self.llx = 0
        if "lly" in self.grid and "lllat" in self.grid:
            raise Exception("Can only respect lly or llat, remove one from config")
        if "lly" in self.grid:
            self.lly = np.int(self.grid["lly"])
        else:
            self.lly = 0
        if "urx" in self.grid and "urlon" in self.grid:
            raise Exception("Can only respect urx or urlon, remove one from config")
        if "urx" in self.grid:
            self.urx = np.int(self.grid["urx"]) + 1
        else:
            self.urx = None
        if "ury" in self.grid and "urlat" in self.grid:
            raise Exception("Can only respect ury or urlat, remove one from config")
        if "ury" in self.grid:
            self.ury = np.int(self.grid["ury"]) + 1
        else:
            self.ury = None
        if "zoom" in self.grid:
            self.zoom = np.int(self.grid["zoom"])
        else:
            self.zoom = False
        self.polygon_check = self.grid.get("polygon_check", False)

        # conventions for reading file
        self.xname = self.conventions["xname"]
        self.yname = self.conventions["yname"]
        self.zname = self.conventions["zname"]
        self.zscale = np.float(
            self.conventions["zscale"]
        )  # combined scale and pos-neg depth convention
        self.bathyscale = np.abs(self.zscale)
        self.xyorder = self.conventions[
            "xyorder"
        ]  # true if bathy variable is x-y array, false if y-x array

        if self.smc:
            self.smctiers = np.int(self.smc["smctiers"])
            self.smcscale = 2.0 ** (self.smctiers - 1.0)

            # smc output file names and info for WW3Meta
            self.WW3Cels = "%sCell.dat" % self.gid.upper()
            self.cellfile = self.workdir + "/" + self.WW3Cels

            # Note that SMC not support subgridscale blocking for multilevel
            # SMC grids. However, a file is needed for consistency. Hence a
            # dummy file with all zero blocking is created here.
            self.WW3Obs = "%sObs.dat" % self.gid.upper()
            self.WW3BPs = "%sBPlist.txt" % self.gid.upper()
            self.WW3Meta = "%s_META.txt" % self.gid.upper()
            self.WW3GDef = "%s_grid.inp" % self.gid.upper()
            self.unitbathy = 30
            self.idlabathy = 3
            self.idfmbathy = 1
            self.smc_dcheck = False
            if "depthlim" in self.smc:
                self.smc_dcheck = True
                self.smc_dlim = np.float(self.smc["depthlim"])
                self.smc_dvar = np.float(self.smc["depthvar"])

        # scale factor used for lat-lon dx-dy - best to hardwire this; using
        # 1.0 will remove a lot of confusion in set-up??
        self.latlonscale = 1.0
        self.llscale = 1.0
        self.depthlim = np.float(self.ww3meta["lsmdepth"])
        self.moddepthmin = np.float(self.ww3meta["mindepth"])
        self.blockscale = np.float(self.ww3meta["blockscale"])
        self.mindepth_switch = self.ww3meta.get("setmindepth", False)

        self.bplist = []
        self.tiers = {}
        self.tx = {}
        self.ty = {}
        self.ntc = {}
        self.actions = {"cell": True, "face": True, "plot": True}
        self.actions.update(actions)
        self.plotting = {"projection": "Robinson",
                         "clon": 180.0,
                         "clat": 0,
                         "fillplot": True,
                         "cellplot": True,
                         "interactive": True}
        self.plotting.update(plotting)

    def read_bathy(self):
        """
        read netCDF file
        """
        d = nc.Dataset(self.bathymetry)
        lat = d.variables[self.yname]
        lon = d.variables[self.xname]
        depth = d.variables[self.zname]
        if self.zoom:
            self.lat = lat[:: self.zoom]
            self.lon = lon[:: self.zoom]
            # self.depth = map_interp(depth[:,:], lon, lat, self.lon, self.lat)
            self.depth = map_interp(depth[:, :], lon, lat, self.lon, self.lat)
            if not ma.is_masked(self.depth):
                self.depth = ma.masked_array(self.depth)
            self.depth.mask = [self.depth > 0 | self.depth.mask]
        else:
            self.lat = lat
            self.lon = lon
            self.depth = depth
        self.lons, self.lats = np.meshgrid(self.lon, self.lat)
        self.lats, self.lons = np.meshgrid(self.lat, self.lon)
        self._set_grid_paramaters()

    def _set_grid_paramaters(self):

        if "lllon" in self.grid:
            self.llx = np.argmin(np.abs(self.lon[:] - float(self.grid["lllon"])))
        if "lllat" in self.grid:
            self.lly = np.argmin(np.abs(self.lat[:] - float(self.grid["lllat"])))
        if "urlon" in self.grid:
            self.urx = np.argmin(np.abs(self.lon[:] - float(self.grid["urlon"])))
        if "urlat" in self.grid:
            self.ury = np.argmin(np.abs(self.lat[:] - float(self.grid["urlat"])))
        if "mergelat" in self.grid:
            mergelat = float(self.grid["mergelat"])
            ury = self.ury or len(self.lat)
            lly = self.lly or 0
            self.mergeny = ury - np.argmin(np.abs(self.lat[:] - mergelat))
            self.mergesy = np.argmin(np.abs(self.lat[:] + mergelat)) - lly
        else:
            self.mergeny = None
            self.mergesy = None
        if "mergelat2" in self.grid:
            mergelat2 = float(self.grid["mergelat2"])
            if mergelat2 <= mergelat:
                raise Exception(
                    "Raise mergelat2 (%s) is less than mergelat (%s)"
                    % (mergelat2, mergelat)
                )
            self.mergeny2 = ury - np.argmin(np.abs(self.lat[:] - mergelat2))
            self.mergesy2 = np.argmin(np.abs(self.lat[:] + mergelat2)) - lly
        else:
            self.mergeny2 = None
            self.mergesy2 = None

    def extract_region(self):
        """
        Extract region of interest
        """

        if self.xyorder == "True":
            self.ingrid_xpts = np.shape(
                self.depth[self.llx:self.urx, self.lly:self.ury])[0]
            self.ingrid_ypts = np.shape(
                self.depth[self.llx:self.urx, self.lly:self.ury])[1]
        else:
            self.ingrid_xpts = np.shape(
                self.depth[self.lly:self.ury, self.llx:self.urx])[1]
            self.ingrid_ypts = np.shape(
                self.depth[self.lly:self.ury, self.llx:self.urx])[0]

        # get grid ll corner and dx,dy values from regular grid file
        if not self.rotated:
            if self.xyorder == "True":
                self.ingrid_lllon = self.lons[self.lly, self.llx]
                self.ingrid_lllat = self.lats[self.lly, self.llx]
                self.dx = lons[self.lly, self.llx + 1] - self.ingrid_lllon
                self.dy = lons[self.lly + 1, self.llx] - self.ingrid_lllat
            else:
                self.ingrid_lllon = self.lons[self.llx, self.lly]
                self.ingrid_lllat = self.lats[self.llx, self.lly]
                self.ingrid_urlon = self.lons[self.urx, self.ury]
                self.ingrid_urlat = self.lats[self.urx, self.ury]
                self.dx = self.lons[self.llx + 1, self.lly] - self.ingrid_lllon
                self.dy = self.lats[self.llx, self.lly + 1] - self.ingrid_lllat
        # sort out the number of x and y cells to actually use - factor of smcscale
        self.use_xpts = np.int(self.ingrid_xpts - self.ingrid_xpts % self.smcscale)
        self.use_ypts = np.int(self.ingrid_ypts - self.ingrid_ypts % self.smcscale)
        self.logger.info(
            "Read in %s xpts; using %s xpts" % (self.ingrid_xpts, self.use_xpts)
        )
        self.logger.info(
            "Read in %s xpts; using %s xpts" % (self.ingrid_ypts, self.use_ypts)
        )

        # calculate the output grid sw corner cell centre
        self.lllon = self.ingrid_lllon - self.dx / 2.0 + (self.smcscale / 2.0) * self.dx
        self.lllat = self.ingrid_lllat - self.dy / 2.0 + (self.smcscale / 2.0) * self.dy

        # calculate the output grid ne corner cell centre
        self.urlon = (
            self.ingrid_lllon
            - self.dx / 2.0
            + (self.use_xpts - self.smcscale / 2.0) * self.dx
        )
        self.urlat = (
            self.ingrid_lllat
            - self.dy / 2.0
            + (self.use_ypts - self.smcscale / 2.0) * self.dy
        )

        self.glob = (((self.lllon + self.urlon + self.dx) % 360) < self.dx / 10.0) or (
            ((self.lllon + self.urlon) % 360) < self.dx / 10.0
        )
        self.logger.info(
            "Located input grid sw corner cell center at: %s,%s"
            % (self.ingrid_lllon, self.ingrid_lllat)
        )
        self.logger.info("\t this gets used for the ww3_grid.inp metadata")
        self.logger.info(
            "Located SMC grid ne corner cell center at: %s,%s"
            % (self.urlon, self.urlat)
        )
        self.logger.info("\t this gets used for the boundary point metadata")
        self.logger.info("Global grid: %s" % self.glob)
        self.logger.info("Approximate resolutions")
        for level in range(self.smctiers):
            res = self.dx * float(2 ** (level))
            self.logger.info("\tTier %s: %s deg (~%s km)" % (level + 1, res, res * 111))
        if self.glob:
            self._check_global_bounds()

    def _check_global_bounds(self):
        """
        Check if grid resolution is compatible with requested SMC leves and
        regions
        """
        self.maxcellwidth = self._determine_max_gridwidth()
        residual = 360.0 % (self.maxcellwidth * self.dx)
        if residual > self.dx / 10.0:
            raise Exception(
                "Grid definitions invalid: Level %s merge of cellwidth %s not divisable by 360"
                % (self.maxcellwidth / self.smcscale / 2, self.maxcellwidth)
            )

    def _determine_max_gridwidth(self, lllat=None, urlat=None):
        """
        Determine maximum gridwidth
        """
        lllat = lllat or self.lllat
        urlat = urlat or self.urlat
        ret = self.smcscale
        for lat in (lllat, urlat):
            if self.mergeny2:
                if abs(lat) > self.mergeny2:
                    ret = max(ret, 4 * self.smcscale)
                    self.merge2_active = True
            if self.mergeny:
                if abs(lat) > self.mergeny:
                    ret = max(ret, 2 * self.smcscale)
                    self.merge_active = True
        return int(ret)

    def define_writedepths(self):
        """
        Establish writedepth array using (lat,lon) convention
        """
        if self.xyorder == True:
            self.writedepths = (
                np.rot90(
                    self.depth[
                        self.llx:self.llx + self.use_xpts,
                        self.lly:self.lly + self.use_ypts,
                    ]
                )
                * self.zscale
            )  # need to check if this line works properly!
        else:
            self.writedepths = (
                self.depth[
                    self.lly:self.lly + self.use_ypts,
                    self.llx:self.llx + self.use_xpts,
                ]
                * self.zscale
            ).filled(999)
        self.logger.info(
            "Dimensions of grid for analysis: %s %s" % (np.shape(self.writedepths))
        )

        self.writelats = self.lats[self.llx, self.lly:self.lly + self.use_ypts]
        self.writelons = self.lons[self.llx:self.llx + self.use_xpts, self.lly]

        if self.mindepth_switch:
            self.writedepths[
                (self.writedepths >= -1.0 * mindepth) & (self.writedepths < 0.0)
            ] = (-1.0 * mindepth)

        if "exclude_flags" in self.grid:
            fll_kwds = dict(lon1d=self.writelons, lat1d=self.writelats)
            user_polys = np.loadtxt(self.grid["exclude_flags"], usecols=(0, 1))
            self.logger.debug("Removing polygons specified in %s" %
                    self.grid["exclude_flags"])
            fnumbers = user_polys[:, 0][user_polys[:, 1] == 1]
            for f in fnumbers:
                p = np.loadtxt(
                    os.path.join(
                        os.path.dirname(Grid.__file__),
                        "user_polygons/user_polygon-%s.txt" % int(f),
                    )
                )
                self.logger.debug("Removing polygon %s" % int(f))
                loc = FindLevelLoc(p, **fll_kwds)
                self.writedepths[loc] = 999  # exclude points inside polygon
        self._define_writedepth_properties()

    def _define_writedepth_properties(self):
        # self.nx = np.shape(self.writedepths)[1] - self.maxcellwidth
        self.nx = np.shape(self.writedepths)[1]
        self.ny = np.shape(self.writedepths)[0]
        self.writemask = np.ones_like(self.writedepths) * 3.0
        self.writeblock = np.zeros([self.ny, self.nx])

    def tier_grid(self):
        """
        Tier up the requested region
        """
        self.logger.info("Analysing tiers")
        self.smcscli = np.int(self.smcscale) * 2
        self._adjacent_land_tier(1)
        self._depth_tier(2)
        for ii in range(3, self.smctiers + 1):
            self._deepwater_tier(ii)

    def _adjacent_land_tier(self, tier):
        """
        1st tier establish locations next to land
        """
        self.tier = tier
        self.logger.info("   ...Analysing tier %s" % self.tier)
        self._check_fixed_tier()
        for self.lpy in tqdm(range(self.ny)):
            for self.lpx in range(0, self.nx, self.cellwidth):
                if not self.writeblock[self.lpy, self.wrap(self.lpx)]:
                    if self.is_land(self.lpy, self.wrap(self.lpx)):
                        self.writemask[self.lpy, self.wrap(self.lpx)] = 0
                        if self.glob:
                            if self.is_land(
                                self.lpy, self.wrap(self.lpx - self.cellwidth)
                            ):
                                self.writemask[
                                    self.lpy, self.wrap(self.lpx - self.cellwidth)
                                ] = 1
                        else:
                            if self.lpx - self.cellwidth >= 0:
                                if self.is_land(
                                    self.lpy, self.wrap(self.lpx - self.cellwidth)
                                ):
                                    self.writemask[
                                        self.lpy, self.wrap(self.lpx - self.cellwidth)
                                    ] = 1
                            if self.lpx + self.cellwidth < self.nx:
                                if self.is_land(
                                    self.lpy, self.wrap(self.lpx + self.cellwidth)
                                ):
                                    self.writemask[
                                        self.lpy, self.wrap(self.lpx + self.cellwidth)
                                    ] = 1
                        if self.lpy - 1 >= 0:
                            if self.is_land(self.lpy - 1, self.wrap(self.lpx)):
                                self.writemask[self.lpy - 1, self.wrap(self.lpx)] = 1
                        if self.lpy + 1 < self.ny:
                            if self.is_land(self.lpy + 1, self.wrap(self.lpx)):
                                self.writemask[self.lpy + 1, self.wrap(self.lpx)] = 1

    def _depth_tier(self, tier):
        """
        Establish second tier.
        If switched on (smc_dcheck) this will retain type 1 cells for
        water below a cut-off depth and where cell-cell depth variability is
        above a threshold
        """
        self.tier = tier
        self.logger.info("   ...Analysing Tier %s" % self.tier)
        self._check_fixed_tier()
        for self.lpy in tqdm(range(0, self.ny, self.cellsize)):
            for self.lpx in range(0, self.nx, self.cellwidth):
                if not np.all(
                    self.writemask[self.cell_indeces(self.lpy, self.lpx)] == 0
                ):
                    if np.any(
                        self.writemask[self.cell_indeces(self.lpy, self.lpx)] == 1
                    ):
                        for sly in range(2):
                            for slx in range(self.cellwidth):
                                if (
                                    self.writemask[
                                        self.lpy + sly, self.wrap(self.lpx + slx)
                                    ]
                                    > 1
                                ):
                                    self.writemask[
                                        self.lpy + sly, self.wrap(self.lpx + slx)
                                    ] = 1
                    elif self.smc_dcheck:
                        # depth based variability
                        if np.any(
                            self.writedepths[self.cell_indeces(self.lpy, self.lpx)]
                            >= -1 * self.smc_dlim
                        ):
                            dmax = np.max(
                                np.abs(
                                    self.writedepths[
                                        self.cell_indeces(self.lpy, self.lpx)
                                    ]
                                )
                            )
                            dmin = np.min(
                                np.abs(
                                    self.writedepths[
                                        self.cell_indeces(self.lpy, self.lpx)
                                    ]
                                )
                            )
                            dmean = np.mean(
                                np.abs(
                                    self.writedepths[
                                        self.cell_indeces(self.lpy, self.lpx)
                                    ]
                                )
                            )
                            ddep = (dmax - dmin) / dmean
                            if (
                                ddep > self.smc_dvar
                                and not self.writeblock[self.lpy, self.wrap(self.lpx)]
                            ):
                                self.writemask[
                                    self.cell_indeces(self.lpy, self.lpx)
                                ] = 1
                            else:
                                self.writemask[
                                    self.cell_indeces(self.lpy, self.lpx)
                                ] = 2
                        else:
                            self.writemask[self.cell_indeces(self.lpy, self.lpx)] = 2
                    else:
                        self.writemask[self.cell_indeces(self.lpy, self.lpx)] = 2

    def _deepwater_tier(self, tier):
        """
        Establish tier above 3
        """
        self.tier = tier
        self.logger.info("Analysing Tier %s" % self.tier)
        for self.lpy in tqdm(range(0, self.ny, self.cellsize)):
            for self.lpx in range(0, self.nx, self.cellwidth):
                # ensure we never go straight from tier 3 to tier 1 by searching over a +/-1 box
                # if self.lpy-1 >= 0 and self.lpy+self.cellsize+1 < self.ny and self.lpx-1 >=0 and self.lpx+self.gridwidth+1 < self.nx:
                if (
                    self.lpy - 1 >= 0
                    and self.lpy + self.cellsize + 1 < self.ny
                    and self.lpx >= 0
                    and self.lpx + self.cellwidth + 1 < self.nx
                ):
                    if np.all(
                        self.writemask[self.cell_indeces(self.lpy, self.lpx, buffer=1)]
                        >= (tier - 1)
                    ):
                        self.writemask[
                            self.cell_indeces(self.lpy, self.lpx)
                        ] = self.tier

    def _check_fixed_tier(self):
        self.writeblock[:] == 0
        if "refine_region" in self.grid:
            for config in self.grid["refine_region"]:
                for method, coords in config.items():
                    fll_kwds = dict(lon1d=self.writelons, lat1d=self.writelats)
                    maxcellwidth = self._determine_max_gridwidth(
                        min(coords["lats"]), max(coords["lats"])
                    )
                    p = np.array(
                        [
                            myround(coords["lons"], maxcellwidth),
                            myround(coords["lats"], self.cellsize),
                        ]
                    )
                    loc = FindLevelLoc(p.T, **fll_kwds)
                    if method == "all":
                        self.writemask[loc] = 1
                    if method == "none":
                        self.writeblock[loc] = 1

    def define_border(self):
        """
        Set border cells at highest tier value
        """
        self.logger.info("Applying highest tier to border cells")

        lpy = 0
        for lpx in range(self.smcscli, self.nx, self.smcscli):
            if np.all(
                self.writemask[
                    lpy:lpy + (self.smcscli + 1), lpx - 1:lpx + (self.smcscli + 1)
                ]
                >= self.smctiers - 1
            ):
                self.writemask[
                    lpy:lpy + self.smcscli, lpx:lpx + self.smcscli
                ] = self.smctiers
                bcy = self.lllat
                bcx = self.lllon + np.float(lpx) * self.dx
                if [bcx, bcy] not in self.bplist:
                    self.bplist.append([bcx, bcy])
        lpy = self.ny
        for lpx in range(self.smcscli, self.nx, self.smcscli):
            if np.all(
                self.writemask[
                    lpy - (self.smcscli + 1):lpy, lpx - 1:lpx + (self.smcscli + 1)
                ]
                >= self.smctiers - 1
            ):
                self.writemask[
                    lpy - self.smcscli:lpy, lpx:lpx + self.smcscli
                ] = self.smctiers
                bcy = self.urlat
                bcx = self.lllon + np.float(lpx) * self.dx
                if [bcx, bcy] not in self.bplist:
                    self.bplist.append([bcx, bcy])
        lpx = 0
        for lpy in range(self.smcscli, self.ny, self.smcscli):
            if np.all(
                self.writemask[
                    lpy - 1:lpy + (self.smcscli + 1), lpx:lpx + (self.smcscli + 1)
                ]
                >= self.smctiers - 1
            ):
                self.writemask[
                    lpy:lpy + self.smcscli, lpx:lpx + self.smcscli
                ] = self.smctiers
                bcy = self.lllat + np.float(lpy) * self.dy
                bcx = self.lllon
                if [bcx, bcy] not in self.bplist:
                    self.bplist.append([bcx, bcy])
        lpx = self.nx
        for lpy in range(self.smcscli, self.ny, self.smcscli):
            if np.all(
                self.writemask[
                    lpy - 1:lpy + (self.smcscli + 1), lpx - (self.smcscli + 1):lpx
                ]
                >= self.smctiers - 1
            ):
                self.writemask[
                    lpy:lpy + self.smcscli, lpx - self.smcscli:lpx
                ] = self.smctiers
                bcy = self.lllat + np.float(lpy) * self.dy
                bcx = self.urlon
                if [bcx, bcy] not in self.bplist:
                    self.bplist.append([bcx, bcy])

        # set corner cells at highest tier value
        if np.all(
            self.writemask[0:(self.smcscli + 1), 0:(self.smcscli + 1)]
            >= self.smctiers - 1
        ):
            self.writemask[0:self.smcscli, 0:self.smcscli] = self.smctiers
            bcy = self.lllat
            bcx = self.lllon
            if [bcx, bcy] not in self.bplist:
                self.bplist.append([bcx, bcy])
        if np.all(
            self.writemask[0:(self.smcscli + 1), -1 * (self.smcscli + 1) :]
            >= self.smctiers - 1
        ):
            self.writemask[0:self.smcscli, -1 * self.smcscli :] = self.smctiers
            bcy = self.lllat
            bcx = self.urlon
            if [bcx, bcy] not in self.bplist:
                self.bplist.append([bcx, bcy])
        if np.all(
            self.writemask[-1 * (self.smcscli + 1) :, 0:(self.smcscli + 1)]
            >= self.smctiers - 1
        ):
            self.writemask[-1 * self.smcscli :, 0:self.smcscli] = self.smctiers
            bcy = self.urlat
            bcx = self.lllon
            if [bcx, bcy] not in self.bplist:
                self.bplist.append([bcx, bcy])
        if np.all(
            self.writemask[-1 * (self.smcscli + 1) :, -1 * (self.smcscli + 1) :]
            >= self.smctiers - 1
        ):
            self.writemask[-1 * self.smcscli :, -1 * self.smcscli :] = self.smctiers
            bcy = self.urlat
            bcx = self.urlon
            if [bcx, bcy] not in self.bplist:
                self.bplist.append([bcx, bcy])

    def _check_writemask(self, iy, ix):
        """
        Check for valid data before adding to tier list
        """
        if not self.tier in self.tiers:
            self.tiers[self.tier] = []
        if not self.tier in self.ty:
            self.ty[self.tier] = []
        if not self.tier in self.tx:
            self.tx[self.tier] = []
        mydepth = np.mean(self.writedepths[self.cell_indeces(iy, ix)])
        if mydepth != "masked" and not np.isnan(mydepth) and not self.is_land(iy, ix):
            self.tiers[self.tier].append(
                [ix, iy, self.cellwidth, self.cellsize, mydepth * -1]
            )
            self.tx[self.tier].append(ix)
            self.ty[self.tier].append(iy)

    def create_cells_lists(self):
        """
        Establish tier lists from masked data
        """

        self.logger.info("Creating cells list for tiers")

        self.tier = self.smctiers

        def loop_quad(ii, refy, refx):
            for lpy in range(0, 2):
                self.tier = ii
                for lpx in range(0, 2):
                    self.tier = ii
                    myx = refx + (lpx * self.cellwidth)
                    myy = refy + (lpy * self.cellsize)
                    if np.all(self.writemask[self.cell_indeces(myy, myx)] == self.tier):
                        self._check_writemask(myy, myx)
                    else:
                        if ii > 1:
                            loop_quad(ii - 1, myy, myx)
                        else:
                            return

        for lpy in tqdm(range(0, self.ny, self.cellsize)):
            self.tier = self.smctiers
            self.lpy = lpy
            if self.glob:
                nxmax = self.nx
            else:
                nxmax = self.nx - self.cellwidth
            for lpx in range(0, nxmax, self.cellwidth):
                self.tier = self.smctiers
                if np.all(self.writemask[self.cell_indeces(lpy, lpx)] == self.tier):
                    self._check_writemask(lpy, lpx)
                else:
                    loop_quad(self.smctiers - 1, lpy, lpx)

    def sort(self):
        """
        Sort lists by y location (2nd value in lists).
        This is the needed order for processing cell face arrays
        """
        self.logger.info("Sorting tier lists")
        self.ttotc = 0
        for ii in range(1, self.smctiers + 1):
            self.tiers[ii].sort(key=lambda x: int(x[1]))
            self.ntc[ii] = len(self.tiers[ii])
            self.ttotc += self.ntc[ii]
            self.logger.info("\tNumber of tier %i cells %i: " % (ii, self.ntc[ii]))

    def plot_scatter(self, proj=ccrs.Robinson(central_longitude=180.0)):
        """
        Plot tier mask
        """
        fig, ax = CartopyMap(proj, coast=True, gridbase=45, figsize=(20, 15))
        lats = self.writelons - self.dy
        lons = self.writelats - self.dx
        # data = self.writemask[:, 0:-self.maxcellwidth]
        data = self.writemask
        m = ax.pcolormesh(lons, lats, data, transform=ccrs.PlateCarree())
        cols = ["k", "b", "g", "r"]
        for ii in range(1, self.smctiers + 1):
            ax.scatter(
                lons[self.tx[ii]],
                lats[self.ty[ii]],
                s=15,
                marker=".",
                color=cols[ii - 1],
                transform=ccrs.PlateCarree(),
            )
        plt.colorbar(m)
        return ax

    def plot_patches(self, filled=False, center=True,
                     proj=ccrs.Robinson(central_longitude=180.0)):
        """
        Plot cell patches
        """
        return plot_patches(
            self.cellfile,
            self.dy,
            self.dx,
            [self.ingrid_lllon, self.ingrid_lllat],
            filled=filled,
            proj=proj,
            center=center,
        )

    def write_cell(self):
        """
        Write out the cells file

        """
        if not os.path.isdir(self.workdir):
            os.makedirs(self.workdir)
        self.logger.info("Writing cell info to " + self.cellfile)
        with open(self.cellfile, "w") as inp:
            if self.smctiers == 2:
                inp.write(
                    " %8d %8d %8d" % tuple([self.ttotc, self.ntc[1], self.ntc[2]])
                    + "  0  0\r\n"
                )
            elif self.smctiers == 3:
                inp.write(
                    " %8d %8d %8d %8d"
                    % tuple([self.ttotc, self.ntc[1], self.ntc[2], self.ntc[3]])
                    + "  0\r\n"
                )
            elif self.smctiers == 4:
                inp.write(
                    " %8d %8d %8d %8d %8d"
                    % tuple(
                        [self.ttotc, self.ntc[1], self.ntc[2], self.ntc[3], self.ntc[4]]
                    )
                    + "  0\r\n"
                )

            for ii in range(1, self.smctiers + 1):
                for lp in range(len(self.tiers[ii])):
                    inp.write(
                        " %5d %5d %2d %2d %5d" % tuple(self.tiers[ii][lp]) + "\r\n"
                    )

        self.logger.info(
            "Writing dummy obstuction info to " + self.workdir + "/" + self.WW3Obs
        )
        self.obsfile = self.workdir + "/" + self.WW3Obs

        with open(self.obsfile, "w") as inp:
            inp.write(" %8d 2\r\n" % self.ttotc)
            for lp in range(self.ttotc):
                inp.write("   0   0 \r\n")

    def write_meta(self):
        """
        Write info to metadata file
        """

        # calculating output metadata here
        gdx = (
            self.smcscale * self.dx / self.latlonscale
        )  # values for grid_def file - based on largest smc cell size
        gdy = (
            self.smcscale * self.dy / self.latlonscale
        )  # values for grid_def file - based on largest smc cell size

        # calculate limits on CFL and 2nd order swell age - based on largest smc cell size
        # maxlat  = (self.lly / self.llscale) + np.float(self.ny) * (self.dy / self.latlonscale)
        ury = (
            self.ury or self.lat.shape[0] - 1
        )  # TODO Fix this to work with xyorder options
        maxlat = max(abs(self.lat[self.lly]), abs(self.lat[ury]))
        minlon = (
            1853.0
            * 60.0
            * (self.smcscale * self.dx / self.latlonscale)
            * np.cos(np.pi * maxlat / 180.0)
        )
        maxcg = 1.4 * 9.81 * 25.0 / (4.0 * np.pi)
        # cflstep = minlon / maxcg
        cflstep = minlon / maxcg

        sagemax = (
            0.5 * minlon ** 2.0 * 12.0 / ((2.0 * np.pi * maxcg / 24.0) ** 2.0 * cflstep)
        )

        self.logger.info("Approximate resolutions")
        for level in range(self.smctiers):
            res = self.dx / float(2 ** (level))
            self.logger.info("\tTier %s: %s deg (~%s km)" % (level + 1, res, res * 111))

        # write grid data to grid.inp metadata file
        # note grid parameters are defined by the samllest cell size
        # and use the small cell centre for the sw corner
        self.logger.info("Writing WW3 metadata to " + self.workdir + "/" + self.WW3Meta)
        with open(self.workdir + "/" + self.WW3Meta, "w") as inp:

            inp.write(
                "$ Grid minimum cell dx: %.2f" % minlon
                + "m at latitude %.3f" % maxlat
                + " degrees\r\n"
            )
            inp.write(
                "$ CFL minimum timestep (needs rounding down): %i" % cflstep
                + " seconds\r\n"
            )
            inp.write(
                "$ Estimated maximum swell age for 24 direction spectrum: %i" % sagemax
                + " seconds\r\n"
            )
            if self.mindepth_switch:
                inp.write("$ Minimum depth set for model at %f" % mindepth + "m\r\n")
            inp.write("$\r\n")
            inp.write(
                "$ Define grid rules -------------------------------------------------- $\r\n"
            )
            inp.write("$ Four records containing :\r\n")
            inp.write(
                "$  1 NX, NY. As the outer grid lines are always defined as land\r\n"
            )
            inp.write("$    points, the minimum size is 3x3.\r\n")
            inp.write(
                "$  2 Grid increments SX, SY (degr.or m) and scaling (division) factor.\r\n"
            )
            inp.write("$    If NX*SX = 360., latitudinal closure is applied.\r\n")
            inp.write(
                "$  3 Coordinates of (1,1) (degr.) and scaling (division) factor.\r\n"
            )
            inp.write(
                "$  4 Limiting bottom depth (m) to discriminate between land and sea\r\n"
            )
            inp.write(
                "$    points, minimum water depth (m) as allowed in model, unit number\r\n"
            )
            inp.write(
                "$    of file with bottom depths, scale factor for bottom depths (mult.),\r\n"
            )
            inp.write(
                "$    IDLA, IDFM, format for formatted read, FROM and filename.\r\n"
            )
            inp.write("$\r\n")
            inp.write(
                "$ Define grid -------------------------------------------------------- $\r\n"
            )
            inp.write("$\r\n")
            inp.write(" %i" % self.nx + " %i" % self.ny + "\r\n")
            inp.write(
                " %8.6f" % self.dx
                + " %8.6f" % self.dy
                + " %5.1f" % self.latlonscale
                + "\r\n"
            )
            inp.write(
                " %9.6f" % self.ingrid_lllon
                + " %9.6f" % self.ingrid_lllat
                + " %5.1f" % self.llscale
                + "\r\n"
            )
            inp.write(
                " %5.2f" % self.depthlim
                + " %5.1f" % self.moddepthmin
                + " %i" % self.unitbathy
                + " %5.1f" % self.bathyscale
                + " %i" % self.idlabathy
                + " %i" % self.idfmbathy
                + " '(....)' 'NAME' '%s" % self.WW3Cels
                + "'\r\n"
            )
            inp.write("$\r\n")

            inp.close()

        # write grid data to grid_def file
        # note grid_def parameters for pre-procesing are defined by the largest cell size
        # and use a largest cell centre for the sw corner
        nxdef = np.int(self.nx / self.smcscale)
        nydef = np.int(self.ny / self.smcscale)
        self.logger.info(
            "Writing grid_def metadata to " + self.workdir + "/" + self.WW3GDef
        )
        with open(self.workdir + "/" + self.WW3GDef, "w") as inp:

            inp.write(" %i" % nxdef + " %i" % nydef + "\r\n")
            inp.write(
                " %9.6f" % self.lllon
                + " %9.6f" % self.lllat
                + " %8.6f" % gdx
                + " %8.6f" % gdy
                + "\r\n"
            )
            # inp.write(' %6.2f' %rlon +' %6.2f' %rlat +'\r\n') # for standard grid set-up

            inp.close()

    def write_bnd(self):
        """
        write out the boundary points file
        """
        self.logger.info(
            "Writing boundary point info (real world lat lons) to "
            + self.workdir
            + "/"
            + self.WW3BPs
        )
        with open(self.workdir + "/" + self.WW3BPs, "w") as inpbp:
            for lp in range(len(self.bplist)):
                writeBP(self.bplist[lp], inpbp, rotated=self.rotated)
            inpbp.close()

    def generate_face_arrays(self):
        """
        Generate face arrays
        """
        self.logger.info('Generating face arrays')
        genCellSides(self.gid.upper(), self.ny, self.nx, self.dy, self.dx,
                     outdir=self.workdir, logger=self.logger)

    @property
    def gridwidth(self):
        """
        Set the gridwidth based on the current latitude and the mergesy and
        mergesy2 values
        """

        if self.mergeny != None:
            if self.lpy < self.mergesy:
                gw = 2
                if self.mergeny2 != None:
                    if self.lpy < self.mergesy2:
                        gw = 4
            elif self.lpy >= self.ny - self.mergeny:
                gw = 2
                if self.mergeny2 != None:
                    if self.lpy >= self.ny - self.mergeny2:
                        gw = 4
            else:
                gw = 1
        else:
            gw = 1
        return gw

    @property
    def cellsize(self):
        """
        Calculate cellsize based on current tier being processed
        """
        return 2 ** (self.tier - 1)

    @property
    def cellwidth(self):
        """
        Calculate current gridwith based on current cellsize and gridwidth
        """
        return self.gridwidth * self.cellsize

    def is_land(self, lpy, lpx, resolution="50m"):
        """
        Determine of cell is on land from coastal polygon dataset
        """
        lat = self.writelats[lpy]
        lon = self.writelons[lpx]
        if self.writedepths[lpy, self.wrap(lpx)] == 999:
            if self.polygon_check:
                if not hasattr(self, "land"):
                    self.logger.info("Creating land shapefile")
                    land_shp_fname = shpreader.natural_earth(
                        resolution=resolution, category="physical", name="land"
                    )
                    land_geom = unary_union(
                        list(shpreader.Reader(land_shp_fname).geometries())
                    )
                    self.land = prep(land_geom)
                if self.land.contains(sgeom.Point(lon, lat)):
                    return True
                else:
                    self.logger.debug(
                        "%s %s excluded based on polygon test" % (lon, lat)
                    )
                    return False
            return True
        else:
            return False

    def cell_indeces(self, lpy, lpx, buffer=0):
        """
        Calculate cell indices for current location based on cellwidth and
        cellsize
        Args:
            lpy, lpy (int): x, y coordinates in regular grid
            buffer (int):   extra buffer around calculated cells that is also
                            requested
        Returns:
            Cell indices from regular grid for SMC cell currently under
            consideration
        """
        return tuple(
            np.meshgrid(
                np.arange(lpy - buffer, lpy + self.cellsize + buffer),
                self.wrap(np.arange(lpx - buffer, lpx + self.cellwidth + buffer)),
            )
        )

    def wrap(self, x):
        """
        Wrap around 360 if required
        """
        if self.glob:
            xwrap = np.arange(self.nx).take(x, mode="wrap")
            if np.array(x).any() > self.nx:
                self.logger.debug("Wrapping %s to %s" % (x, xwrap))
            return xwrap
        else:
            return x

    def run(self):
        self.read_bathy()
        self.extract_region()
        self.define_writedepths()
        if self.actions["cell"]:
            self.tier_grid()
            self.define_border()
            self.create_cells_lists()
            self.sort()
            self.write_cell()
            self.write_meta()
            self.write_bnd()
        if self.actions["plot"]:
            # self.plot_scatter()
            if self.plotting['projection'] == 'Orthographic':
                proj = getattr(ccrs,
                        self.plotting['projection'])(
                            central_longitude=self.plotting['clon'],
                            central_latitude=self.plotting['clat']
                            )
            else:
                proj = getattr(ccrs,
                        self.plotting['projection'])(
                            central_longitude=self.plotting['clon'],
                            )
            if self.plotting["fillplot"] == True:
                self.plot_patches(filled=True, proj=proj, center=False)
                plt.savefig(os.path.join(self.workdir, self.gid + "_filled.png"))
            if self.plotting["cellplot"] == True:
                self.plot_patches(filled=False, proj=proj, center=True)
                plt.savefig(os.path.join(self.workdir, self.gid + "_cells.png"))
        if self.actions["face"]:
            self.generate_face_arrays()
        if self.plotting["interactive"] == True:
            plt.show()


class GRIDGEN2SMC(NC2SMC):

    def __init(self, *args, **kwargs):
        super(GRIDGEN2SMC, self).__init__(*args, **kwargs)

    def read_bathy(self):
        from scipy.io import loadmat

        matDict = loadmat(self.bathymetry, squeeze_me=True)
        keys = [
            "dlon",
            "dlat",
            "lon",
            "lat",
            "depth",
            "m3",
            "m4",
            "mask_map",
            "sx1",
            "sy1",
        ]
        for key in keys:
            setattr(self, key, matDict[key])
        del keys, matDict
        self.lons = self.lon.transpose()
        self.lats = self.lat.transpose()
        self.lon = self.lons[:, 0]
        self.lat = self.lats[0, :]
        self.depth = ma.masked_array(self.depth)
        self._set_grid_paramaters()


def plot_patches( cellfile, dlat, dlon, refp, filled=False,
        proj=ccrs.PlateCarree(central_longitude=180.0), center=True,
        linewidth=0.25, coast=True, land=True):
    plot_kws = dict(txtloc=(0.25, 0.88),
                    dotsize=0.15,
                    cax_kws=dict(width="2.5%",
                                 height="90%",
                                 borderpad=4.5,
                                 loc=6,
                                 bbox_to_anchor=(0.975, 0.0, 1, 1),),
                                 txtSize=5,
                                 cb_kws=dict(orientation="vertical"),
                                 cbtxtSize=5,
                                 linewidth=linewidth)
    glbSMC = UnSMC(cellfile, dlon=dlon, dlat=dlat, refp=refp)
    fig, ax = CartopyMap(proj, coast=coast, land=land, gridbase=45, figsize=(20, 15))
    glbSMC.genPlot(filled=filled, ax=ax, plot_var="depth", center=center,
            **plot_kws)
    return ax


def map_interp(datain, xin, yin, xout, yout, interpolation="Bilinear"):

    """
       Interpolates a 2D array onto a new grid (only works for linear grids),
       with the Lat/Lon inputs of the old and new grid. Can perfom nearest
       neighbour interpolation or bilinear interpolation (of order 1)'

       This is an extract from the basemap module (truncated)
    """

    # Mesh Coordinates so that they are both 2D arrays
    xout, yout = np.meshgrid(xout, yout)

    # compute grid coordinates of output grid.
    delx = xin[1:] - xin[0:-1]
    dely = yin[1:] - yin[0:-1]

    xcoords = (len(xin) - 1) * (xout - xin[0]) / (xin[-1] - xin[0])
    ycoords = (len(yin) - 1) * (yout - yin[0]) / (yin[-1] - yin[0])

    xcoords = np.clip(xcoords, 0, len(xin) - 1)
    ycoords = np.clip(ycoords, 0, len(yin) - 1)

    # Interpolate to output grid using nearest neighbour
    if interpolation == "NearestNeighbour":
        xcoordsi = np.around(xcoords).astype(np.int32)
        ycoordsi = np.around(ycoords).astype(np.int32)
        dataout = datain[ycoordsi, xcoordsi]

        # Interpolate to output grid using bilinear interpolation.
    elif interpolation == "Bilinear":
        xi = xcoords.astype(np.int32)
        yi = ycoords.astype(np.int32)
        xip1 = xi + 1
        yip1 = yi + 1
        xip1 = np.clip(xip1, 0, len(xin) - 1)
        yip1 = np.clip(yip1, 0, len(yin) - 1)
        delx = xcoords - xi.astype(np.float32)
        dely = ycoords - yi.astype(np.float32)
        dataout = (
            (1.0 - delx) * (1.0 - dely) * datain[yi, xi]
            + delx * dely * datain[yip1, xip1]
            + (1.0 - delx) * dely * datain[yip1, xi]
            + delx * (1.0 - dely) * datain[yi, xip1]
        )

    return dataout


def FindLevelLoc(size_bbox=None, lon1d=None, lat1d=None):
    """
    Locate all the points inside the `size_bbox`.

    Input args:
        size_bbox -- a rectangle or irregular polygon defining the extent of
                     one specific size/level

        lon1d & lat1d -- the x&y axis of bathy array (nlat x nlon)
    """
    size_bbox = np.asarray(size_bbox)  # [lon, lat]
    size_path = Path(size_bbox)

    lon2d, lat2d = np.meshgrid(lon1d, lat1d)
    bathy_points = np.vstack((lon2d.ravel(), lat2d.ravel())).T

    inside = size_path.contains_points(bathy_points).reshape(lon2d.shape)

    return inside


def myround(x, base=5):
    return base * np.round(np.array(x).astype(float) / float(base))


def genCellSides(gid, nlat, nlon, dlat, dlon, outdir='./', logger=logging):
    """
    genCellSides

    Generate cell sides using wrapped fortran code

    Args:
        gid (str):          Grid id
        nlat, nlon (int):   lat/lon dimensions of grid
        dlat, dlon (float): lat/lon resolution of grid
        outdir (str):       output directory
        logger (str):       logger
    """
    gencellside.adapgrid(os.path.join(outdir, gid), nlat, nlon, dlat, dlon)
    SortFaceArray(os.path.join(outdir, gid + "ISide.d"),
                  os.path.join(outdir, gid + "JSide.d"), logger=logger)


def SortFaceArray(ufnm, vfnm, logger=logging):
    """
    SortFaceArray

    Sort U face array and V face array

    args:
        ufnm & vfnm -- fnm of unsorted face/side array
    """
    logger.info("Sorting face arrays...")
    uSide = pd.read_csv(
        ufnm, header=None, sep="\s+", names=["i", "j", "ysize", "ll", "l", "r", "rr"]
    )

    vSide = pd.read_csv(
        vfnm,
        header=None,
        sep="\s+",
        names=["i", "j", "xsize", "bb", "b", "t", "tt", "ysize"],
    )

    for fnm, Side in zip([ufnm, vfnm], [uSide, vSide]):
        Side.sort_values(by=["ysize", "j", "i"], inplace=True)

        NS = Side["ysize"].size
        N1 = np.equal(Side["ysize"], 1).sum()
        N2 = np.equal(Side["ysize"], 2).sum()
        N4 = np.equal(Side["ysize"], 4).sum()
        N8 = np.equal(Side["ysize"], 8).sum()

        fnm = fnm.replace(".d", ".dat")
        fmt = "{:10d} " * 5
        header = fmt[:-1].format(NS, N1, N2, N4, N8)
        np.savetxt(fnm, np.array(Side), fmt="%10d", header=header, comments="")

def writeBP(bcs, inpbp, rotated=False):

    bcx = bcs[0]
    bcy = bcs[1]

    if rotated:
        bplon, bplat = iris.analysis.cartography.unrotate_pole(
            np.array([bcx]), np.array([bcy]), rlon, rlat
        )
        bcx = bplon[0]
        bcy = bplat[0]

    if bcx < 0.0:
        bcx = 360.0 + bcx

    # inpbp.write( '%8.3f' %bcx + ',%8.3f' %bcy + '%8.3f' %bcxr + ',%8.3f' %bcyr + ' 0.0 0.0 1\r\n' )
    inpbp.write("%8.3f" % bcx + ",%8.3f" % bcy + " 0.0 0.0 1\r\n")

    return



if __name__ == "__main__":
    config = load_config("test_global.yml")
    smc = NC2SMC(**config)
    smc.run()
    proj = ccrs.Orthographic(central_longitude=0, central_latitude=60)
    ax = smc.plot_patches(proj=proj)
    plt.show()
