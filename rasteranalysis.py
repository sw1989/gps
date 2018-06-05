#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      simon
#
# Created:     05/06/2018
# Copyright:   (c) simon 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------


import georasters as gr
import numpy as np
import arcpy
import time
import pyproj
import shapely
from shapely.ops import transform
from functools import partial
import matplotlib.pyplot as plt




#prepare distance raster

def generatedistraster(listofbbgvalues=[40,42]):
    #0) Setting computing environment
    arcpy.env.overwriteOutput = True #Such that files can be overwritten
    arcpy.env.workspace = r"C:\\Users\\simon\\Documents\\GitHub\\gps" #Setting the workspace
    if arcpy.CheckExtension("Spatial") == "Available": #Check out spatial analyst extension
        arcpy.CheckOutExtension("Spatial")
    listofbbgvalues = (str(t) for t in listofbbgvalues)
    # Local variables:
    distrast = "distrast"
    BBG2012_Publicatiebestand = 'C:/Temp/BestandBodemgebruik2012/BestandBodemgebruik2012.gdb/Publicatie/BBG2012_Publicatiebestand'
    arcpy.env.extent = BBG2012_Publicatiebestand



    # Process: Select Layer By Attribute
    arcpy.MakeFeatureLayer_management(BBG2012_Publicatiebestand, 'BBGlayer')
    arcpy.SelectLayerByAttribute_management('BBGlayer', "NEW_SELECTION", "BG2012 IN ("+','.join(listofbbgvalues)+")")

    # Process: Euclidean Distance
    arcpy.gp.EucDistance_sa(BBG2012_Publicatiebestand, distrast, "", "100")
    try:
        arcpy.RasterToOtherFormat_conversion(distrast,arcpy.env.workspace,"TIFF")
    except:
        print "Raster To Other Format exsample failed."
        print arcpy.GetMessages()

project = lambda x, y: pyproj.transform(pyproj.Proj(init='EPSG:4326'), pyproj.Proj(init='EPSG:28992'), x, y)

def lookupWGS84(x,y,raster="distrast.tif"):
    p = shapely.geometry.point.Point(x,y)
    print(str(p))
    pn = transform(project, p)
    print(str(pn))
    # Load data
    data = gr.from_file(raster)
    value = data.map_pixel(pn.x,pn.y)
    print(str(value))
    plt.figure(figsize=(12, 8))
    data.plot()




def main():
##    start = time.time()
##    generatedistraster()
##    end = time.time()
##    print("Duration = "+str(end - start))

    lookupWGS84(6.29, 52.11)


if __name__ == '__main__':
    main()
