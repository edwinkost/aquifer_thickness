#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Edwin Husni Sutanudjaja (EHS, 29 Sep 2014): This script for reporting aquifer properties 
#                                             at 5 arc-min and 30 arc-min resolution. 

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
#~ output_directory      = "/scratch/edwin/aquifer_properties/" 
#~ output_05min_filename = "/scratch/edwin/aquifer_properties/groundwater_properties_05min.nc"
#~ output_30min_filename = "/scratch/edwin/aquifer_properties/groundwater_properties_30min.nc"
output_directory         = "/p/1209496-pcrglobwb/scratch_edwin/test_aquifer_properties/" 
output_05min_filename    = output_directory + "/groundwater_properties_05min.nc"
#~ output_30min_filename = output_directory + "/groundwater_properties_30min.nc"
cleanOutputDir   = True

# netcdf attributes:
netcdf_attributes = {}
netcdf_attributes['institution']  = "Utrecht University, Dept. of Physical Geography"
netcdf_attributes['title'      ]  = "Groundwater properties" 
netcdf_attributes['source'     ]  = "None" 
netcdf_attributes['history'    ]  = "version 29 September 2014" 
netcdf_attributes['references' ]  = "None" 
netcdf_attributes['description']  = "None" 
netcdf_attributes['comment'    ]  = "Processed and calculated by Edwin H. Sutanudjaja (e-mail: E.H.Sutanudjaja@uu.nl" 

# clone/landmask maps at 05 arc min and 30 arc min resolution
# - for the entire world
#~ clone_map_05min_file = "/scratch/edwin/processing_whymap/version_19september2014/water_polygon/water-polygons-split-4326/landmask_05min.map"
#~ clone_map_30min_file = "/data/hydroworld/PCRGLOBWB20/input30min/routing/lddsound_30min.map"
# - Rhine Meuse example:
clone_map_05min_file    = "/p/1209496-pcrglobwb/pcrglobwb_data/dfguu/data/hydroworld/others/RhineMeuse/RhineMeuse05min.landmask.map"
#~ clone_map_30min_file    = "/p/1209496-pcrglobwb/pcrglobwb_data/dfguu/data/hydroworld/others/RhineMeuse/RhineMeuse30min.landmask.map"

# input file: thickness properties at 30arc min resolution
thickness_05min_netcdf = {}
# - output from the script to estimate "aquifer thickness"
thickness_05min_netcdf['filename']    = "/p/1209496-pcrglobwb/scratch_edwin/test_aquifer_thickness/sedimentary_basin_05_arcmin.nc" 
# - using the corrected values based on Margat's table
thickness_05min_netcdf['var_name']    = "average_corrected"
#~ # - using the non-corrected values
#~ thickness_05min_netcdf['var_name'] = "average"

# input file: aquifer properties at 5arc min resolution  - from Gleeson et al.
aquifer_properties_05min_netcdf = {}
aquifer_properties_05min_netcdf['filename'] = "/p/1209496-pcrglobwb/pcrglobwb_data/dfguu/data/hydroworld/PCRGLOBWB20/input5min/groundwater/groundwaterProperties5ArcMin.nc"

#~ # input file: aquifer properties at 30arc min resolution - from Gleeson et al.
#~ aquifer_properties_30min_netcdf = {}
#~ aquifer_properties_30min_netcdf['filename'] = "/p/1209496-pcrglobwb/pcrglobwb_data/dfguu/data/hydroworld/PCRGLOBWB20/input30min/groundwater/groundwaterProperties.nc"


def main():

    # make output directory
    try:
        os.makedirs(output_directory) 
    except:
        if cleanOutputDir == True: os.system('rm -r '+output_directory+"/*")

    # change the current directory/path to output directory
    os.chdir(output_directory)
    
    # make temporary directory
    tmp_directory = output_directory+"/tmp"
    os.makedirs(tmp_directory)     
    vos.clean_tmp_dir(tmp_directory)
    
    # format and initialize logger
    logger_initialize = Logger(output_directory)

    logger.info('Start processing for 5 arc-min resolution!')
    
    # clone and landmask for 5 arc-min resolution
    pcr.setclone(clone_map_05min_file)
    landmask = pcr.defined(clone_map_05min_file)
    
    # read thickness value (at 5 arc min resolution)
    logger.info('Reading the thickness at 5 arc-min resolution!')
    thickness = pcr.ifthen(landmask,\
                vos.netcdf2PCRobjCloneWithoutTime(thickness_05min_netcdf['filename'], "average_corrected", clone_map_05min_file))
    #
    # update landmask
    landmask = pcr.defined(thickness)            
    
    # read aquifer properties at 5 arc min resolution
    logger.info('Reading saturated conductivity and specific yield at 5 arc-min resolution!')
    saturated_conductivity = pcr.ifthen(landmask,\
                             vos.netcdf2PCRobjSCloneWithoutTime(\
                             aquifer_properties_S05min_netcdf['filename'],\
                             "kSatAquifer"  , clone_map_05min_file))
    specific_yield         = pcr.ifthen(landmask,\
                             vos.netcdf2PCRobjCloneWithoutTime(\
                             aquifer_properties_05min_netcdf['filename'],\
                             "specificYield", clone_map_05min_file))

    # saving 5 min parameters to a netcdf file
    logger.info('Saving groundwater parameter parameters to a netcdf file: '+output_05min_filename)
    #
    output_05min_netcdf = outputNetCDF.OutputNetCDF(clone_map_05min_file)
    #
    variable_names = ["saturated_conductivity","specific_yield","thickness"]
    units = ["m/day","1","m"]
    variable_fields = [
                       pcr.pcr2numpy(saturated_conductivity, vos.MV),
                       pcr.pcr2numpy(specific_yield        , vos.MV),
                       pcr.pcr2numpy(thickness             , vos.MV),
                       ]
    pcr.report(saturated_conductivity, "saturated_conductivity_05min.map")
    pcr.report(specific_yield        , "specific_yield_05min.map")
    pcr.report(thickness             , "thickness_05min.map")
    output_05min_netcdf.createNetCDF(   output_05min_filename,variable_names,units)
    output_05min_netcdf.changeAtrribute(output_05min_filename,netcdf_attributes)
    output_05min_netcdf.data2NetCDF(    output_05min_filename,variable_names,variable_fields)

    #~ logger.info('Start processing for 30 arc-min resolution!')

    #~ # upscaling thickness to 30 arc min resolution
    #~ logger.info('Upscaling thickness from 5 arc-min resolution to 30 arc-min!')
    #~ thickness_05min_array = pcr.pcr2numpy(thickness, vos.MV)
    #~ thickness_30min_array = vos.regridToCoarse(thickness_05min_array, 30./5.,"average")
    #~ #
    #~ # set clone for 30 arc min resolution
    #~ pcr.setclone(clone_map_30min_file)
    #~ #
    #~ landmask  = pcr.defined(clone_map_30min_file)
    #~ thickness = pcr.ifthen(landmask, 
                #~ pcr.numpy2pcr(pcr.Scalar, thickness_30min_array, vos.MV))
    #~ #
    #~ # update landmask
    #~ landmask  = pcr.defined(thickness)            

    #~ # read aquifer properties at 30 arc min resolution
    #~ logger.info('Reading saturated conductivity and specific yield at 30 arc-min resolution!')
    #~ saturated_conductivity = pcr.ifthen(landmask,\
                             #~ vos.netcdf2PCRobjCloneWithoutTime(\
                             #~ aquifer_properties_30min_netcdf['filename'],\
                             #~ "kSatAquifer"  , clone_map_30min_file))
    #~ specific_yield         = pcr.ifthen(landmask,\
                             #~ vos.netcdf2PCRobjCloneWithoutTime(\
                             #~ aquifer_properties_30min_netcdf['filename'],\
                             #~ "specificYield", clone_map_30min_file))

    #~ # saving 30 min parameters to a netcdf file
    #~ logger.info('Saving groundwater parameter parameters to a netcdf file: '+output_30min_filename)
    #~ #
    #~ output_30min_netcdf = outputNetCDF.OutputNetCDF(clone_map_30min_file)
    #~ #
    #~ variable_names = ["saturated_conductivity","specific_yield","thickness"]
    #~ units = ["m/day","1","m"]
    #~ variable_fields = [
                       #~ pcr.pcr2numpy(saturated_conductivity, vos.MV),
                       #~ pcr.pcr2numpy(specific_yield        , vos.MV),
                       #~ pcr.pcr2numpy(thickness             , vos.MV),
                       #~ ]
    #~ pcr.report(saturated_conductivity, "saturated_conductivity_30min.map")
    #~ pcr.report(specific_yield        , "specific_yield_30min.map")
    #~ pcr.report(thickness             , "thickness_30min.map")
    #~ output_30min_netcdf.createNetCDF(   output_30min_filename,variable_names,units)
    #~ output_30min_netcdf.changeAtrribute(output_30min_filename,netcdf_attributes)
    #~ output_30min_netcdf.data2NetCDF(    output_30min_filename,variable_names,variable_fields)


if __name__ == '__main__':
    sys.exit(main())
