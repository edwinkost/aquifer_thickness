#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pcraster.framework import *
import pcraster as pcr

import logging
logger = logging.getLogger(__name__)

import virtualOS as vos

class MonteCarloAquiferThickness(DynamicModel, MonteCarloModel):

    def __init__(self, clone_map_file, \
                       dem_average_netcdf, dem_floodplain_netcdf, ldd_netcdf, \
                       lookup_table_average_thickness, lookup_table_zscore, \
                       number_of_samples, include_percentile = True,\
                       threshold_sedimentary_basin = 50.0, elevation_F_min = 0.0, elevation_F_max = 50.0):  # values defined in de Graaf et al. (2014)

        DynamicModel.__init__(self)
        MonteCarloModel.__init__(self)
        
        msg  = "\n"
        msg += "\n"
        msg += 'For each step, please refer to de Graaf et al. (2014).'
        msg += "\n"
        logger.info(msg)

        # set clone
        self.clone_map_file = clone_map_file
        pcr.setclone(clone_map_file)
        
        # an option to include_percentile or not
        self.include_percentile = include_percentile
        
        # number of samples
        self.number_of_samples = pcr.scalar(number_of_samples)
        
        logger.info("Step 1: Identify the cells belonging to the sedimentary basin region.")
        dem_average    = vos.netcdf2PCRobjCloneWithoutTime(dem_average_netcdf['file_name'],\
                                                           dem_average_netcdf['variable_name'],\
                                                           clone_map_file)
        dem_average    = pcr.max(0.0, dem_average)
        self.dem_average = pcr.cover(dem_average, 0.0)

        dem_floodplain = vos.netcdf2PCRobjCloneWithoutTime(dem_floodplain_netcdf['file_name'],
                                                           dem_floodplain_netcdf['variable_name'],
                                                           clone_map_file)
        dem_floodplain = pcr.max(0.0, dem_floodplain)
        
        lddMap = vos.netcdf2PCRobjCloneWithoutTime(ldd_netcdf['file_name'],
                                                   ldd_netcdf['variable_name'],
                                                   clone_map_file)
        self.lddMap   = pcr.lddrepair(pcr.lddrepair(pcr.ldd(lddMap)))
        self.landmask = pcr.defined(self.lddMap)

        elevation_F    = dem_average - dem_floodplain
        
        sedimentary_basin_extent = pcr.ifthen(elevation_F < pcr.scalar(threshold_sedimentary_basin), pcr.boolean(1))
        
        # include the continuity along the river network
        sedimentary_basin_extent = pcr.windowmajority(sedimentary_basin_extent, 3.00*vos.getMapAttributes(clone_map_file,"cellsize"))
        sedimentary_basin_extent = pcr.cover(sedimentary_basin_extent, \
                                   pcr.path(self.lddMap, pcr.defined(sedimentary_basin_extent)))  

        # TODO: We should also include the extent of major aquifer basins and unconsolidated sediments in the GLiM map.

        elevation_F    = pcr.ifthen(sedimentary_basin_extent, elevation_F)
        elevation_F    = pcr.max(0.0, elevation_F)
        elevation_F    = pcr.min(50., elevation_F)
        
        logger.info("Step 2: Calculate relative difference and associate z_score.")
        relative_elevation_F = pcr.scalar(1.0) - (elevation_F - elevation_F_min)/(elevation_F_max - elevation_F_min)
        
        z_score_relat_elev_F = pcr.lookupscalar(lookup_table_zscore, relative_elevation_F)  
        self.F = z_score_relat_elev_F   # zscore (varying over the map)
        
        # maximum and minimum z_score
        self.F = pcr.min(  3.75, self.F)
        self.F = pcr.max(-10.00, self.F)

        logger.info("Step 3: Assign average and variation of aquifer thickness.")
        self.lookup_table_average_thickness = lookup_table_average_thickness
        self.lnCV = pcr.scalar(0.1)                                          # According to Inge, this lnCV value corresponds to the table "lookup_table_average_thickness".

    def premcloop(self):
        pass 

    def initial(self):
        pass 
      
    def dynamic(self):

        logger.info("Step 4: Monte Carlo simulation")

        # draw a random value (uniform for the entire map)
        z = pcr.mapnormal() ; #~ self.report(z,"z")
        
        # constraints, in order to make sure that random values are in the table of "lookup_table_average_thickness" 
        z = pcr.max(-5.0, z)
        z = pcr.min( 5.0, z)
        
        # assign average thickness (also uniform for the entire map) based on z
        self.Davg = pcr.lookupscalar(self.lookup_table_average_thickness, z)
        #
        self.report(self.Davg,"davg")
        self.lnDavg = pcr.ln(self.Davg)
      	
        # sedimentary basin thickness (varying over cells and samples)
        lnD = self.F * (self.lnCV * self.lnDavg) + self.lnDavg
        
        # set the minimum depth (must be bigger than zero)        
        minimum_depth = 0.005
        lnD = pcr.max( pcr.ln(minimum_depth), lnD)
        
        # extrapolation 
        lnD = pcr.cover(lnD, \
              pcr.windowaverage(lnD, 1.50*vos.getMapAttributes(self.clone_map_file,"cellsize")))
        lnD = pcr.cover(lnD, \
              pcr.windowaverage(pcr.cover(lnD, pcr.ln(minimum_depth)), 3.00*vos.getMapAttributes(self.clone_map_file,"cellsize")))
        lnD = pcr.cover(lnD, \
              pcr.windowaverage(pcr.cover(lnD, pcr.ln(minimum_depth)), 0.50))
        
        # smoothing per quarter arc degree
        lnD = pcr.windowaverage(lnD, 0.25)

        # thickness in meter
        self.D = pcr.exp(lnD)
 
        #~ # smoothing  bottom elevation 
        #~ dem_bottom = pcr.windowaverage(self.dem_average - self.D, 0.50)
        #~ # thickness in meter
        #~ self.D = pcr.max(0.0, self.dem_average - dem_bottom)

       #~ # smoothing 
        #~ self.D = pcr.windowaverage(self.D, 1.50*vos.getMapAttributes(self.clone_map_file,"cellsize"))
        
        # accuracy until cm only
        self.D = pcr.rounddown(self.D*100.)/100.
        
        self.report(self.D, "damc")

    def postmcloop(self):
    
        logger.info("Step 5: Reporting the results.")

        names= ["damc"]
        mcaveragevariance(names, self.sampleNumbers(), self.timeSteps())
        
        if self.include_percentile:
            percentiles = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
            mcpercentiles(names, percentiles, self.sampleNumbers(), self.timeSteps())

        #~ # debugging average values
        #~ os.system("aguila damc-ave.001")
        
        self.average          = pcr.readmap("damc-ave.001")
        self.average_variance = pcr.readmap("damc-var.001")
        
        self.total_variance = self.average_variance * self.number_of_samples
        self.standard_deviation = (self.total_variance / (self.number_of_samples - 1.0) ) ** (0.5)
        
        if self.include_percentile:
            self.percentiles = {}                                                       
            self.percentileList = percentiles                                           
            for percentile in percentiles:
                filename = "damc_1_%1.1f.map" %(percentile) ; print filename
                self.percentiles[percentile] = pcr.readmap(filename)
