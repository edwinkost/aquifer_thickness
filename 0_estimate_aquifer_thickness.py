#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Edwin Husni Sutanudjaja (EHS, 29 Sep 2014): This is the script for estimating aquifer thickness based on de Graaf et al. (2014). 
#                                             This estimate is then corrected based on the table of Margat (2006). 

import os
import sys
import numpy as np

from pcraster.framework import *
import pcraster as pcr

import outputNetCDF
import virtualOS as vos
from monte_carlo_thickness import MonteCarloAquiferThickness
import margat_correction 

import logging
logger = logging.getLogger("main_script") # get name for the logger
from logger import Logger

# output directory: 
#~ output_directory = "/scratch/edwin/aquifer_thickness_5arcmin_world_final/" 
output_directory = "/p/1209496-pcrglobwb/scratch_edwin/test_aquifer_thickness" 
cleanOutputDir   = True

# clone map 
#~ clone_map_file = "/data/hydroworld/others/RhineMeuse/RhineMeuse05min.clone.map"
#~ clone_map_file = "/data/hydroworld/PCRGLOBWB20/input5min/routing/lddsound_05min.map"
#~ clone_map_file = "/scratch/edwin/tmp/Australia05min.clone.map"
#~ clone_map_file = "/scratch/edwin/processing_whymap/version_19september2014/water_polygon/water-polygons-split-4326/landmask_05min.map"
clone_map_file = "/p/1209496-pcrglobwb/pcrglobwb_data/dfguu/data/hydroworld/others/Java/Java05min.clone.map"

# number_of_samples and option to include_percentile_report
number_of_samples         = 5 # 1000
number_of_cores           = 1 # 10
include_percentile_report = False

# sedimentary basin output file:
sedimentary_basin_netcdf = {}
sedimentary_basin_netcdf['file_name']                 = "sedimentary_basin_05_arcmin.nc"
sedimentary_basin_netcdf['attribute']                 = {} 
sedimentary_basin_netcdf['attribute']['institution']  = "Utrecht University, Dept. of Physical Geography"
sedimentary_basin_netcdf['attribute']['title'      ]  = "Aquifer/sedimentary basin thickness at 5 arcmin resolution"
sedimentary_basin_netcdf['attribute']['source'     ]  = "None"
sedimentary_basin_netcdf['attribute']['history'    ]  = "None"
#
sedimentary_basin_netcdf['attribute']['references' ]  = "Margat, 2006. "
sedimentary_basin_netcdf['attribute']['references' ] += "Sutanudjaja et al, 2011. "
sedimentary_basin_netcdf['attribute']['references' ] += "Sutanudjaja et al, 2014. "
sedimentary_basin_netcdf['attribute']['references' ] += "de Graaf et al., 2014. "
#
sedimentary_basin_netcdf['attribute']['description']  = "Aquifer/sedimentary basin thickness derived based on the method of de Graaf et al. (2014) "
sedimentary_basin_netcdf['attribute']['description'] += "and then corrected based on the reported values from Margat (2003)."
#
sedimentary_basin_netcdf['attribute']['comment'    ]  = "Processed and calculated by Edwin H. Sutanudjaja on 29 September 2014."

# input files:
#
dem_average_netcdf = {}
#~ dem_average_netcdf['file_name']     = "/scratch/edwin/floodplain_5arcmin_world_final/topography_parameters_05_arcmin.nc"
dem_average_netcdf['file_name']        = "/p/1209496-pcrglobwb/pcrglobwb_data/dfguu/users/edwin/data/floodplain_05arcmin_world_final/based_on_daily_runoff/topography_parameters_5_arcmin_october_2015.nc"
dem_average_netcdf['variable_name']    = "dem_average"
#
dem_floodplain_netcdf = {}
#~ dem_floodplain_netcdf['file_name']  = "/scratch/edwin/floodplain_5arcmin_world_final/channel_parameters_05_arcmin.nc"
dem_floodplain_netcdf['file_name']     = "/p/1209496-pcrglobwb/pcrglobwb_data/dfguu/users/edwin/data/floodplain_05arcmin_world_final/based_on_daily_runoff/channel_parameters_5_arcmin_october_2015.nc"
dem_floodplain_netcdf['variable_name'] = "dem_floodplain"
#
ldd_netcdf = {}
ldd_netcdf['file_name']                = "/p/1209496-pcrglobwb/pcrglobwb_data/dfguu/users/edwin/data/floodplain_05arcmin_world_final/based_on_daily_runoff/channel_parameters_5_arcmin_october_2015.nc"
ldd_netcdf['variable_name']            = "lddMap"

# table/figure of Margat
margat_aquifers = {}
#~ margat_aquifers['shapefile'] = "/scratch/edwin/processing_whymap/version_19september2014/whymap_wgs1984.shp"
#~ margat_aquifers['txt_table'] = "/scratch/edwin/processing_whymap/version_19september2014/table/margat_table.txt"
margat_aquifers['shapefile']    = "/p/1209496-pcrglobwb/pcrglobwb_data/dfguu/users/edwin/data/processing_whymap/version_19september2014/whymap_wgs1984.shp"
margat_aquifers['txt_table']    = "/p/1209496-pcrglobwb/pcrglobwb_data/dfguu/users/edwin/data/processing_whymap/version_19september2014/table/margat_table.txt"

# TODO: include the parameterization of kSat and Sy

def main():

    # make output directory
    try:
        os.makedirs(output_directory) 
    except:
        if cleanOutputDir == True: os.system('rm -r '+output_directory+"/*")

    # path for Inge's table
    table_path = str(os.getcwd()+"/")
    print(table_path)
    
    table_thickness = table_path + "table_from_inge/lookupDepth.txt"
    table_zscore    = table_path + "table_from_inge/zscore.txt"
    
    # change the current directory/path to output directory
    os.chdir(output_directory)
    
    # make temporary directory
    tmp_directory = output_directory+"/tmp"
    os.makedirs(tmp_directory)     
    vos.clean_tmp_dir(tmp_directory)
    
    # format and initialize logger
    logger_initialize = Logger(output_directory)
    # Monte Carlo simulation
    #
    logger.info('Performing Monte Carlo simulation to estimate aquifer properties !!!')
    #
    myModel = MonteCarloAquiferThickness(clone_map_file, \
                                         dem_average_netcdf, dem_floodplain_netcdf, ldd_netcdf, \
                                         table_thickness, table_zscore, \
                                         number_of_samples, include_percentile_report)

    dynamic_framework = DynamicFramework(myModel,1)
    mcModel = MonteCarloFramework(dynamic_framework, nrSamples=number_of_samples)
    mcModel.setForkSamples(True, nrCPUs=number_of_cores)
    mcModel.run()
    
    # report average, average variance, standard deviation and percentiles to netcdf files
    #
    sedimentary_basin_netcdf['file_name'] = vos.getFullPath(sedimentary_basin_netcdf['file_name'], output_directory)
    logger.info('Reporting topography parameters to a netcdf file: '+sedimentary_basin_netcdf['file_name'])
    #
    sed_bas_netcdf = outputNetCDF.OutputNetCDF(clone_map_file)
    #
    variable_names = ["average","average_variance","standard_deviation"]
    units = ["m","m2","m"]
    pcr.setclone(clone_map_file)
    variable_fields = [pcr.pcr2numpy(myModel.average           , vos.MV),
                       pcr.pcr2numpy(myModel.average_variance  , vos.MV),
                       pcr.pcr2numpy(myModel.standard_deviation, vos.MV)
                       ]
    sed_bas_netcdf.createNetCDF(   sedimentary_basin_netcdf['file_name'],variable_names,units)
    sed_bas_netcdf.changeAtrribute(sedimentary_basin_netcdf['file_name'],sedimentary_basin_netcdf['attribute'])
    sed_bas_netcdf.data2NetCDF(    sedimentary_basin_netcdf['file_name'],variable_names,variable_fields)
    #
    # reporting percentile values
    if include_percentile_report:
        pcr.setclone(clone_map_file)
        for percentile in myModel.percentileList:
            variable_name = "percentile%04d" %(int(percentile*100)) ; print(variable_name)
            variable_field = pcr.pcr2numpy(myModel.percentiles[percentile], vos.MV)
            variable_unit = "m"
            sed_bas_netcdf.addNewVariable(sedimentary_basin_netcdf['file_name'],variable_name,variable_unit)
            sed_bas_netcdf.data2NetCDF(   sedimentary_basin_netcdf['file_name'],variable_name,variable_field)

    # Correcting or rescaling aquifer thickness map based on Margat's table
    logger.info("Correcting/rescaling based on the table of Margat and van der Gun")
    #
    MargatCorrection = margat_correction.MargatCorrection(clone_map_file,\
                                                          sedimentary_basin_netcdf['file_name'],\
                                                          "average",\
                                                          margat_aquifers,
                                                          tmp_directory)
    average_corrected = pcr.pcr2numpy(\
                        MargatCorrection.aquifer_thickness, vos.MV)
    #
    # saving corrected aquifer thickness value to the netcdf file
    sed_bas_netcdf.addNewVariable(sedimentary_basin_netcdf['file_name'],"average_corrected","m")
    sed_bas_netcdf.data2NetCDF( sedimentary_basin_netcdf['file_name'],"average_corrected",average_corrected)

if __name__ == '__main__':
    sys.exit(main())
