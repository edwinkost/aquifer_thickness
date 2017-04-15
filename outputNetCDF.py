#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import datetime
import time
import re
import subprocess
import netCDF4 as nc
import numpy as np
import pcraster as pcr
import virtualOS as vos

class OutputNetCDF():
    
    def __init__(self, cloneMapFileName_or_latlonDict, attributeDictionary = None):
        		
        # cloneMap
        if isinstance( cloneMapFileName_or_latlonDict, str):
            # define latitudes and longitudes based on cloneMap
            cloneMapFileName = cloneMapFileName_or_latlonDict
            pcr.setclone(cloneMapFileName)
            cloneMap = pcr.boolean(1.0)
            # latitudes and longitudes
            self.latitudes  = np.unique(pcr.pcr2numpy(pcr.ycoordinate(cloneMap), vos.MV))[::-1]
            self.longitudes = np.unique(pcr.pcr2numpy(pcr.xcoordinate(cloneMap), vos.MV))
        else:
            # define latitudes and longitudes based on latlonDict       # NOT TESTED YET
            latlonDict = cloneMapFileName_or_latlonDict
            self.latitudes  = latlonDict['lat']
            self.longitudes = latlonDict['lon']

        # make sure that latitudes are from high to low
        if self.latitudes[-1]  >  self.latitudes[0]: self.latitudes = self.latitudes[::-1]
        if self.longitudes[-1] < self.longitudes[0]: self.longitudes = self.longitudes[::-1]
        
        # netcdf format:
        self.format = 'NETCDF3_CLASSIC'
        
        self.attributeDictionary = {}
        if attributeDictionary == None:
            self.attributeDictionary['institution'] = "None"
            self.attributeDictionary['title'      ] = "None"
            self.attributeDictionary['source'     ] = "None"
            self.attributeDictionary['history'    ] = "None"
            self.attributeDictionary['references' ] = "None"
            self.attributeDictionary['description'] = "None"
            self.attributeDictionary['comment'    ] = "None"
        else:
            self.attributeDictionary = attributeDictionary
        
    def createNetCDF(self,ncFileName,varName,varUnit=None,varLongName=None,timeAttribute=None):

        rootgrp= nc.Dataset(ncFileName,'w',format= self.format)

        #-create dimensions - time is unlimited, others are fixed
        rootgrp.createDimension('lat',len(self.latitudes))
        rootgrp.createDimension('lon',len(self.longitudes))

        lat= rootgrp.createVariable('lat','f4',('lat',))
        lat.long_name = 'latitude'
        lat.units = 'degrees_north'
        lat.standard_name = 'latitude'

        lon= rootgrp.createVariable('lon','f4',('lon',))
        lon.standard_name = 'longitude'
        lon.long_name = 'longitude'
        lon.units = 'degrees_east'

        lat[:] = self.latitudes
        lon[:] = self.longitudes

        if timeAttribute != None:                                                                                                                                                      
            rootgrp.createDimension('time',None)                                                                                                                                       
            date_time= rootgrp.createVariable('time','f4',('time',))                                                                                                                   
            date_time.standard_name = 'time'                                                                                                                                            
            date_time.long_name = 'Days since 1901-01-01'                                                                                                                               
            date_time.units = 'Days since 1901-01-01'                                                                                                                                   
            date_time.calendar = 'standard'                                                                                                                                             

        # variable short and long names
        if isinstance(varName,list) == False: varName = [varName] 
        if varLongName == None: varLongName = varName 
        if isinstance(varLongName,list) == False: varLongName = [varLongName]
        
        # variable unit
        if varUnit == None: varUnit = [None] * length(varName)
        if isinstance(varUnit, list) == False: varUnit = [varUnit] 
        
        for i in range(0, len(varName)):                                                                                                                                         
            shortVarName = varName[i]                                                                                                                                            
            longVarName  = varLongName[i]                                                                                                                                          
            unitVar      = varUnit[i]                                                                                                                                                  
            if unitVar == None: unitVar = 'undefined'                                                                                                                                                  
            if timeAttribute != None:                                                                                                                                                  
                var= rootgrp.createVariable(shortVarName,'f4',('time','lat','lon',) ,fill_value=vos.MV,zlib=False)                                                                      
            else:                                                                                                                                                                      
                var= rootgrp.createVariable(shortVarName,'f4',('lat','lon',) ,fill_value=vos.MV,zlib=False)                                                                             
            var.standard_name = shortVarName                                                                                                                                           
            var.long_name = longVarName                                                                                                                                                
            var.units = unitVar

        for k, v in self.attributeDictionary.items(): setattr(rootgrp,k,v)

        rootgrp.sync()
        rootgrp.close()

    def addNewVariable(self,ncFileName,varName,varUnit=None,varLongName=None,timeAttribute=None):

        rootgrp= nc.Dataset(ncFileName,'a',format= self.format)

        # variable short and long names
        if isinstance(varName,list) == False: varName = [varName] 
        if varLongName == None: varLongName = varName 
        if isinstance(varLongName,list) == False: varLongName = [varLongName]
        
        # variable unit
        if varUnit == None: varUnit = [None] * length(varName)
        if isinstance(varUnit, list) == False: varUnit = [varUnit] 
        
        for i in range(0, len(varName)):                                                                                                                                         
            shortVarName = varName[i]                                                                                                                                            
            longVarName  = varLongName[i]                                                                                                                                          
            unitVar      = varUnit[i]                                                                                                                                                  
            if unitVar == None: unitVar = 'undefined'                                                                                                                                                  
            if timeAttribute != None:                                                                                                                                                  
                var = rootgrp.createVariable(shortVarName,'f4',('time','lat','lon',) ,fill_value=vos.MV,zlib=False)                                                                      
            else:                                                                                                                                                                      
                var = rootgrp.createVariable(shortVarName,'f4',('lat','lon',) ,fill_value=vos.MV,zlib=False)                                                                             
            var.standard_name = shortVarName                                                                                                                                           
            var.long_name = longVarName                                                                                                                                                
            var.units = unitVar

        rootgrp.sync()
        rootgrp.close()

    def changeAtrribute(self,ncFileName,attributeDictionary):

        rootgrp= nc.Dataset(ncFileName,'a',format= self.format)

        for k, v in attributeDictionary.items():
          setattr(rootgrp,k,v)

        rootgrp.sync()
        rootgrp.close()

    def data2NetCDF(self,ncFile,varName,varField,timeStamp=None,posCnt=None):

        #-write data to netCDF
        rootgrp= nc.Dataset(ncFile,'a')    

        if isinstance(varName,list) == False: varName = [varName] 
        if isinstance(varField,list) == False: varField = [varField]

        # index for time
        if timeStamp != None:
            date_time = rootgrp.variables['time']
            if posCnt == None: posCnt = len(date_time)
            date_time[posCnt] = nc.date2num(timeStamp,date_time.units,date_time.calendar)	                                                                                                                                                  

        for i in range(0, len(varName)):
            shortVarName = varName[i]
            if timeStamp != None:                                                                                                                                                  
                rootgrp.variables[shortVarName][posCnt,:,:] = varField[i]                                                                      
            else:                                                                                                                                                                      
                rootgrp.variables[shortVarName][:,:]        = varField[i]

        rootgrp.sync()
        rootgrp.close()
