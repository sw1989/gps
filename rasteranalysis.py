#-------------------------------------------------------------------------------
# Name:         rasteranalysis.py
# Purpose:      This scipt enriches GPS data with distances to landuse areas
#
# Author:      simon
#
# Created:     05/06/2018
# Copyright:   (c) simon 2018
# Licence:
#-------------------------------------------------------------------------------

import numpy as np
import time
import os
import csv

#These geo modules require the installation of GDAL and other libraries. for windows binaries, see https://www.lfd.uci.edu/~gohlke/pythonlibs/
import georasters as gr
#see https://github.com/ozak/georasters
import pyproj
#see https://github.com/jswhit/pyproj
import shapely
#See https://github.com/Toblerity/Shapely
from shapely.ops import transform

#import arcpy
#This module is not open source (ArcGIS), but can be left away if rasters are not generated




"""generates distance raster files for diverse landuse categories on the BBG dataset"""
def generatedistraster(listofbbgvalues=[40,42], name="parks"):
    #0) Setting computing environment
    arcpy.env.overwriteOutput = True #Such that files can be overwritten
    arcpy.env.workspace = r"C:\\Users\\simon\\Documents\\GitHub\\gps\\distrast" #Setting the workspace
    if arcpy.CheckExtension("Spatial") == "Available": #Check out spatial analyst extension
        arcpy.CheckOutExtension("Spatial")
    listofbbgvalues = (str(t) for t in listofbbgvalues)
    # Local variables:

    BBG2012_Publicatiebestand = 'C:/Temp/BestandBodemgebruik2012/BestandBodemgebruik2012.gdb/Publicatie/BBG2012_Publicatiebestand'
    arcpy.env.extent = BBG2012_Publicatiebestand

    # Process: Select landuse Layer By Attribute
    landuseregions = os.path.join(arcpy.env.workspace, 'landuse.shp')
    arcpy.MakeFeatureLayer_management(BBG2012_Publicatiebestand, 'BBGlayer')
    arcpy.SelectLayerByAttribute_management('BBGlayer', "NEW_SELECTION", "BG2012 IN ("+','.join(listofbbgvalues)+")")
    arcpy.CopyFeatures_management('BBGlayer', landuseregions)

    # Process: Euclidean Distance
    distrast = os.path.join(arcpy.env.workspace,"dist2"+name)
    arcpy.gp.EucDistance_sa(landuseregions, distrast, '', "100", '')
    arcpy.Delete_management(landuseregions)
    try:
        arcpy.RasterToOtherFormat_conversion(distrast,arcpy.env.workspace,"TIFF")
        arcpy.Delete_management(distrast)
    except:
        print "Raster To Other Format exsample failed."
        print arcpy.GetMessages()
    return distrast+'.tif'



"""A function for looking up Raster cell row/colum projected in RD_new, based on a WGS84 coordinate pair as input as well as a Geo raster tile"""
project = lambda x, y: pyproj.transform(pyproj.Proj(init='EPSG:4326'), pyproj.Proj(init='EPSG:28992'), x, y)
def lookupWGS84(x,y,GeoT):
    p = shapely.geometry.point.Point(x,y)
    print(str(p))
    pn = transform(project, p)
    print(str(pn))
    try:
        # Find location of point (x,y) on raster, e.g. to extract info at that location
        col, row = gr.map_pixel(pn.x,pn.y,GeoT[1],GeoT[-1], GeoT[0],GeoT[3])
        #value = gr.map_pixel(pn.x,pn.y, )
    except:
        print "Coordinates out of bounds!"
        return 'NN','NN'
    return row, col

"""Enriches a track csv file with distances to landuse categories, both given as input. Writes out a novel csv file. If coordinates are out of bounds, returns NN"""
def enrichtrack(csvf, categories):
    clist = categories.keys()
    #testpoint = 123907.175,  496407.706
    #print("\n  Mean: {}, Std: {}, max: {}, min: {}".format(str(data.mean()), str(data.std()),str(data.max()),str(data.min())))
    NDV, xsize, ysize, GeoT, Projection, DataType = gr.get_geo_info("distrast\\dist2sports.tif")
    data = gr.load_tiff("distrast\\dist2sports.tif")
    s = data.shape
    rows = s[0]
    cols = s[1]
    #load raster data once
    rasters = []
    for c in clist:
        rasters.append(gr.load_tiff("distrast\\dist2"+c+".tif"))
    with open(csvf, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        with open(csvf.split('.')[0]+'en.csv', 'wb') as out:
            writer = csv.writer(out, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['Y', 'X']+clist)
            for i,line in enumerate(reader):
                x = float(line[1])
                y = float(line[0])
                valuelist =[]
                row, col = lookupWGS84(x, y, GeoT)
                if row != 'NN' and (row < rows and col < cols):
                    for data in rasters:
                        valuelist.append(str(data[row,col]))
                    writer.writerow([line[0], line[1]]+valuelist)
                else:
                    print("Coordinates out of bounds!")
                    writer.writerow([line[0], line[1]]+['NN' for i in clist])
                #if i ==20:
                #    break
    csvfile.close()
    out.close()







def main():
    bbgcategories = {'parks':[40,42], 'agric':[51], 'sports':[41], 'recr': [43,44], 'frst':[60], 'ntr':[61,62], 'wtr':[70,71,72,73,74,75,76,77,78,80,81,82,83], 'trfc':[11]}
    """meaning of BBG 2012 landuse categories: agric (agriculture areas), sports (sports areas) recr (recreation areas),frst (forest areas), ntr (nature areas), wtr (water areas), trfc (traffic areas)"""
    start = time.time()
##    for name,v in bbgcategories.items():
##        generatedistraster(v,name)
    #store distances to these landuse areas for points in some track
    enrichtrack('GPS.csv',bbgcategories)
    end = time.time()
    print("Duration = "+str(end - start))






if __name__ == '__main__':
    main()
