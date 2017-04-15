#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np

from pcraster.framework import *
import pcraster as pcr

import logging
logger = logging.getLogger(__name__)

import virtualOS as vos

class MargatCorrection(object):

    def __init__(self, clone_map_file,\
                       input_thickness_netcdf_file,\
                       input_thickness_var_name   ,\
                       margat_aquifers,\
                       tmp_directory,
                       landmask = None,
                       arcdegree = True):

        object.__init__(self)

        # aquifer table from Margat and van der Gun 
        self.margat_aquifers = margat_aquifers

        # clone map
        self.clone_map_file = clone_map_file
        self.clone_map_attr = vos.getMapAttributesALL(self.clone_map_file)
        if arcdegree == True:
            self.clone_map_attr['cellsize'] = round(self.clone_map_attr['cellsize'] * 360000.)/360000.
        xmin = self.clone_map_attr['xUL']
        xmax = xmin + self.clone_map_attr['cols'] * self.clone_map_attr['cellsize']
        ymax = self.clone_map_attr['yUL']
        ymin = ymax - self.clone_map_attr['rows'] * self.clone_map_attr['cellsize']
        pcr.setclone(self.clone_map_file)

        # temporary directory 
        self.tmp_directory = tmp_directory

        # thickness approximation (unit: m, file in netcdf with variable name = average 
        self.approx_thick = vos.netcdf2PCRobjCloneWithoutTime(input_thickness_netcdf_file,\
                                                              input_thickness_var_name,\
                                                              self.clone_map_file)
        # set minimum value to 0.1 mm
        self.approx_thick = pcr.max(0.0001, self.approx_thick)

        # rasterize the shape file 
        #               -        
        # save current directory and move to temporary directory
        current_dir = str(os.getcwd()+"/")
        os.chdir(str(self.tmp_directory))
        #
        cmd_line  = 'gdal_rasterize -a MARGAT '                                     # layer name = MARGAT
        cmd_line += '-te '+str(xmin)+' '+str(ymin)+' '+str(xmax)+' '+str(ymax)+ ' '       
        cmd_line += '-tr '+str(self.clone_map_attr['cellsize'])+' '+str(self.clone_map_attr['cellsize'])+' '
        cmd_line += str(margat_aquifers['shapefile'])+' '
        cmd_line += 'tmp.tif'
        print(cmd_line); os.system(cmd_line)
        #
        # make it nomial
        cmd_line = 'pcrcalc tmp.map = "nominal(tmp.tif)"' 
        print(cmd_line); os.system(cmd_line)
        #
        # make sure that the clone map is correct
        cmd_line = 'mapattr -c '+str(self.clone_map_file)+' tmp.map'
        print(cmd_line); os.system(cmd_line)
        #
        # read the map
        self.margat_aquifer_map = pcr.nominal(pcr.readmap("tmp.map"))
        #
        # clean temporary directory and return to the original directory
        vos.clean_tmp_dir(self.tmp_directory)
        os.chdir(current_dir)
        
        # extend the extent of each aquifer
        self.margat_aquifer_map = pcr.cover(self.margat_aquifer_map, 
                                  pcr.windowmajority(self.margat_aquifer_map, 1.25))

        # assign aquifer thickness, unit: m (lookuptable operation) 
        self.margat_aquifer_thickness = pcr.lookupscalar(margat_aquifers['txt_table'], self.margat_aquifer_map)
        self.margat_aquifer_thickness = pcr.ifthen(self.margat_aquifer_thickness > 0., \
                                                   self.margat_aquifer_thickness)
        #~ pcr.report(self.margat_aquifer_thickness,"thick.map"); os.system("aguila thick.map")

        # aquifer map
        self.margat_aquifer_map       = pcr.ifthen(self.margat_aquifer_thickness > 0., self.margat_aquifer_map)        
        
        # looping per aquifer: cirrecting or rescaling 
        aquifer_ids = np.unique(pcr.pcr2numpy(pcr.scalar(self.margat_aquifer_map), vos.MV))
        aquifer_ids = aquifer_ids[aquifer_ids > 0]
        aquifer_ids = aquifer_ids[aquifer_ids < 10000]
        self.rescaled_thickness = None
        for id in aquifer_ids:
            rescaled_thickness = self.correction_per_aquifer(id)
            try:
                self.rescaled_thickness = pcr.cover(self.rescaled_thickness, rescaled_thickness)
            except:
                self.rescaled_thickness = rescaled_thickness
        
        # integrating
        ln_aquifer_thickness  = self.mapFilling( pcr.ln(self.rescaled_thickness), pcr.ln(self.approx_thick) )
        self.aquifer_thickness = pcr.exp(ln_aquifer_thickness)
        #~ pcr.report(self.aquifer_thickness,"thick.map"); os.system("aguila thick.map")

        # cropping only in the landmask region
        if landmask == None: landmask = self.clone_map_file
        self.landmask = pcr.defined(vos.readPCRmapClone(landmask,self.clone_map_file,self.tmp_directory))
        #~ pcr.report(self.landmask,"test.map"); os.system("aguila test.map")

        self.aquifer_thickness = pcr.ifthen(self.landmask, self.aquifer_thickness)
        #~ pcr.report(self.aquifer_thickness,"thick.map"); os.system("aguila thick.map")

    def correction_per_aquifer(self, id):
        
        id = float(id); print id
        
        # identify aquifer mask  
        aquifer_landmask = pcr.ifthen(self.margat_aquifer_map == pcr.nominal(id), pcr.boolean(1))
        
        # obtain the logarithmic value of Margat value
        exp_margat_thick = pcr.cellvalue(\
                           pcr.mapmaximum(\
                           pcr.ifthen(aquifer_landmask, pcr.ln(self.margat_aquifer_thickness))), 1)[0]

        # obtain the logarithmic values of 'estimated thickness'
        exp_approx_thick = pcr.ifthen(aquifer_landmask, pcr.ln(self.approx_thick)) 
                       
        exp_approx_thick_array = pcr.pcr2numpy(exp_approx_thick, vos.MV)
        exp_approx_thick_array = exp_approx_thick_array[exp_approx_thick_array <> vos.MV]
        exp_approx_thick_array = exp_approx_thick_array[exp_approx_thick_array < 1000000.]
        
        # identify percentile
        exp_approx_minim = np.percentile(exp_approx_thick_array,  2.5);
        exp_approx_maxim = np.percentile(exp_approx_thick_array, 97.5); 

        # correcting
        exp_approx_thick_correct  = ( exp_approx_thick - exp_approx_minim ) / \
                                    ( exp_approx_maxim - exp_approx_minim )   
        exp_approx_thick_correct  = pcr.max(0.0, exp_approx_thick_correct )
        exp_approx_thick_correct *= pcr.max(0.0,\
                                    ( exp_margat_thick - exp_approx_minim ) )
        exp_approx_thick_correct += pcr.min(exp_approx_minim, exp_approx_thick)

        # maximum thickness
        exp_approx_thick_correct  = pcr.min(exp_margat_thick, exp_approx_thick_correct)

        # corrected thickness
        correct_thickness = pcr.exp(exp_approx_thick_correct)
        
        return correct_thickness
      
    def mapFilling(self, map_with_MV, map_without_MV, method = "window_average"):
    
        # ----- method 1: inverse distance method (but too slow)
        if method == "inverse_distance":
            logger.info('Extrapolation using "inverse distance" in progress!')
            #
            # - interpolation mask for cells without values
            interpolatedMask = pcr.ifthenelse(\
                               pcr.defined(map_with_MV),\
                               pcr.boolean(0),\
                               pcr.boolean(1),)
            map_with_MV_intrpl = pcr.inversedistance(interpolatedMask, \
                                                   map_with_MV, 2, 1.50, 25)
        #
        else: # method 2: using window average
            logger.info('Extrapolation using "modified window average" in progress!')
            #
            map_with_MV_intrpl = 0.70 * pcr.windowaverage(map_with_MV, 1.50) + \
                                 0.25 * pcr.windowaverage(map_with_MV, 2.00) + \
                                 0.05 * pcr.windowaverage(map_with_MV, 2.50) + \
                                 pcr.scalar(0.0)
        #
        # - interpolated values are only introduced in cells with MV 
        map_with_MV_intrpl = pcr.cover(map_with_MV, map_with_MV_intrpl)
        #
        # - calculating weight factor:
        weight_factor = pcr.scalar(pcr.defined(map_with_MV))
        weight_factor = pcr.windowaverage(0.70*weight_factor, 1.50) +\
                        pcr.windowaverage(0.25*weight_factor, 2.00) +\
                        pcr.windowaverage(0.05*weight_factor, 2.50)
        weight_factor = pcr.min(1.0, weight_factor)
        weight_factor = pcr.max(0.0, weight_factor)
        weight_factor = pcr.cover(weight_factor, 0.0)
        #
        # merge with weight factor
        merged_map = weight_factor  * map_with_MV_intrpl + \
              (1.0 - weight_factor) * map_without_MV
        #
        # retain the original values and make sure that all values are covered
        filled_map = pcr.cover(map_with_MV, merged_map)
        filled_map = pcr.cover(filled_map, map_without_MV)
    
        logger.info('Extrapolation is done!')
        return filled_map

