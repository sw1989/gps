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
import pandas

#These geo modules require the installation of GDAL and other libraries. for windows binaries, see https://www.lfd.uci.edu/~gohlke/pythonlibs/
import georasters as gr
import geopandas as gp
import rasterstats
#see https://github.com/ozak/georasters
import pyproj
#see https://github.com/jswhit/pyproj
import shapely
#See https://github.com/Toblerity/Shapely
from shapely.ops import transform
from shapely.geometry import shape
import rtree
import fiona


#import arcpy
#from arcpy.sa import *

#This module is not open source (ArcGIS), but can be left away if rasters are not generated




"""generates distance raster files for diverse landuse categories on the BBG dataset"""
def generatedistraster(listofbbgvalues=[40,42], name="parks"):
    #0) Setting computing environment
    arcpy.env.overwriteOutput = True #Such that files can be overwritten
    ##arcpy.env.workspace = r"C:\\Users\\schei008\\Documents\\GitHub\\gps\\distrast" #Setting the workspace
    arcpy.env.workspace = r"C:\\Users\\simon\\Documents\\GitHub\\gps\\distrast" #Setting the workspace

    if arcpy.CheckExtension("Spatial") == "Available": #Check out spatial analyst extension
        arcpy.CheckOutExtension("Spatial")
    listofbbgvalues = (str(t) for t in listofbbgvalues)
    # Local variables:

    #BBG2012_Publicatiebestand = 'C:/Users/schei008/surfdrive/Temp/BestandBodemgebruik2012/BestandBodemgebruik2012.gdb/Publicatie/BBG2012_Publicatiebestand'
    BBG2012_Publicatiebestand = 'C:/Users/simon/surfdrive/Temp/BestandBodemgebruik2012/BestandBodemgebruik2012.gdb/Publicatie/BBG2012_Publicatiebestand'
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

"""generates raster that captures the coverage of a landuse class in a neighborhood around each cell"""
def generateCoverageRaster(name="parks"):
    arcpy.env.overwriteOutput = True #Such that files can be overwritten
    #arcpy.env.workspace = r"C:\\Users\\schei008\\Documents\\Github\\gps\\distrast" #Setting the workspace
    arcpy.env.workspace = r"C:\\Users\\simon\\Documents\\Github\\gps\\distrast" #Setting the workspace
    if arcpy.CheckExtension("Spatial") == "Available": #Check out spatial analyst extension
        arcpy.CheckOutExtension("Spatial")
    #distrast = os.path.join(r"C:\\Users\\schei008\\Documents\\Github\\gps\\distrast","dist2"+name)
    distrast = os.path.join(r"C:\\Users\\simon\\Documents\\Github\\gps\\distrast","dist2"+name)
    #reclrast  = os.path.join(r"C:\\Users\\schei008\\Documents\\Github\\gps\\distrast","in2"+name)
    reclrast  = os.path.join(r"C:\\Users\\simon\\Documents\\Github\\gps\\distrast","in2"+name)
    covrast = os.path.join(arcpy.env.workspace,"covOf"+name)

    #Turns cells within landuse region into 1, and all other into 0
    recl = arcpy.sa.Reclassify(distrast+'.tif', "VALUE", arcpy.sa.RemapValue([[0,1]]), "NODATA")
    recl2 = Con(IsNull(recl), 0, recl)

    recl2.save(reclrast)
    arcpy.Delete_management(recl)

    #Generates the number (= percentage) of landuse rasters within a focal window (of 1 km * 1km) around each cell
    outfocal = arcpy.sa.FocalStatistics(reclrast,arcpy.sa.NbrRectangle(10,10,"CELL"), "SUM")
    outfocal.save(covrast)
    try:
        arcpy.RasterToOtherFormat_conversion(covrast,arcpy.env.workspace,"TIFF")
        arcpy.Delete_management(covrast)
        arcpy.Delete_management(reclrast)
        arcpy.Delete_management(outfocal)
    except:
        print "Raster To Other Format exsample failed."
        print arcpy.GetMessages()
    return covrast+'.tif'






"""A function for looking up Raster cell row/colum projected in RD_new, based on a WGS84 coordinate pair as input as well as a Geo raster tile"""
project = lambda x, y: pyproj.transform(pyproj.Proj(init='EPSG:4326'), pyproj.Proj(init='EPSG:28992'), x, y)
def lookupWGS84(x,y,raster):
    p = shapely.geometry.point.Point(x,y)
    print(str(p))
    x,y = project(x,y)
    print(x,y)
    try:
        # Find location of point (x,y) on raster, e.g. to extract info at that location
        #col, row = gr.map_pixel(pn.x,pn.y,GeoT[1],GeoT[-1], GeoT[0],GeoT[3])
        value = raster.map_pixel(x,y)
    except:
        print "Coordinates out of bounds!"
        return 'NN'#,'NN'
    return value

"""Enriches a track csv file with distances to landuse categories, both given as input. Writes out a novel csv file. If coordinates are out of bounds, returns NN"""
def enrichtrack(csvf, categories, rastertype='dist2'):
    clist = categories.keys()
    #testpoint = 123907.175,  496407.706
    #print("\n  Mean: {}, Std: {}, max: {}, min: {}".format(str(data.mean()), str(data.std()),str(data.max()),str(data.min())))
    #NDV, xsize, ysize, GeoT, Projection, DataType = gr.get_geo_info("distrast\\dist2rsda.tif")
    #data = gr.load_tiff("distrast\\dist2rsda.tif")
    #s = data.shape
    #rows = s[0]
    #cols = s[1]
    #load raster data once
    rasters = []
    for c in clist:
        rasters.append(gr.from_file("distrast\\"+rastertype+c+".tif"))
    with open(csvf, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        with open(csvf.split('.')[0]+'en.csv', 'wb') as out:
            writer = csv.writer(out, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['Y', 'X']+clist)
            for i,line in enumerate(reader):
                x = float(line[1])
                y = float(line[0])
                valuelist =[]
                #row, col = lookupWGS84(x, y, GeoT)
                for data in rasters:
                    value = lookupWGS84(x, y, data)
                    if value != 'NN': #and (row < rows and col < cols):
                        valuelist.append(str(value))
                    else:
                        print("Coordinates out of bounds!")
                        valuelist.append(str('NN'))
                    writer.writerow([line[0], line[1]]+valuelist)
                #if i ==20:
                #    break
    csvfile.close()
    out.close()


#see https://blog.maptiks.com/spatial-queries-in-python/
def generate_index(records, index_path=None):
    prop = rtree.index.Property()
    if index_path is not None:
        prop.storage = rtree.index.RT_Disk
        prop.overwrite = index_path

    sp_index = rtree.index.Index(index_path, properties=prop)
    for n, ft in enumerate(records):
        if ft['geometry'] is not None:
            sp_index.insert(int(n), shape(ft['geometry']).bounds)
    return sp_index

#Selects a polygon in a projected layer based on x y coordinates in WGS84
def selectPolygonData(x,y, index, layer):
    p = shapely.geometry.point.Point(x,y)
    print(str(p))
    pn = transform(project, p)
    print(str(pn))
    for fid in index.intersection(pn.bounds):
        if fid != None and fid in layer.keys():
            feat = layer[int(fid)]
            #print feat
            shp2 = shape(feat['geometry'])
            if shp2.intersects(pn):
                    return feat
                    break
    return None


def enrichCBS(csvf):
    attributes = ['BEV_DICHTH','P_65_EO_JR','P_EENP_HH', 'P_HH_Z_K','P_N_W_AL']
    """Population density, percentage people over 65 years, percentage one person households, percentage households without children, percentage non-western immigrants"""
    layer2 = 'stats\\buurt_2014.shp'
    with fiona.open(layer2) as datasrc2:
        fname, _ = os.path.splitext(layer2)
        if os.path.exists(fname + '.dat') and os.path.exists(fname + '.idx'):
            layer2_index = rtree.index.Index(fname)
        else:
            layer2_index = generate_index(datasrc2, fname)

        with open(csvf, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            with open(csvf.split('.')[0]+'en.csv', 'wb') as out:
                writer = csv.writer(out, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                writer.writerow(['Y', 'X']+attributes)
                for i,line in enumerate(reader):
                    x = float(line[1])
                    y = float(line[0])
                    feat = selectPolygonData(x, y, layer2_index, datasrc2)
                    if feat != None:
                        attr = feat['properties']
                        #print("row: "+line[0]+" "+ line[1]+" "+str(attr['BEV_DICHTH']) )
                        writer.writerow([line[0], line[1]]+[str(attr[a]) for a in attributes])
                    else:
                        print("Coordinates out of bounds!")
                        writer.writerow([line[0], line[1]]+['NN' for i in attributes])
                    #if i ==5:
                    #    break
        csvfile.close()
    datasrc2.close()


def generateCellStats(categories, rastertype='dist2', sourcefolder= 'distrast'):

    clist = categories.keys()
    #testpoint = 123907.175,  496407.706
    #print("\n  Mean: {}, Std: {}, max: {}, min: {}".format(str(data.mean()), str(data.std()),str(data.max()),str(data.min())))
    #load raster data once
    rasters = []
    NDV, xsize, ysize, GeoT, Projection, DataType = gr.get_geo_info(sourcefolder+"\\"+rastertype+clist[0]+".tif")
    netherlands = gp.read_file('Netherlands_withoutwater_s.shp')
    netherlands_w = gp.read_file('Netherlands.shp')
    stats = []
    import matplotlib.pyplot as plt
    for c in clist:
        print sourcefolder+"\\"+rastertype+c+".tif"
        r =gr.from_file(sourcefolder+"\\"+rastertype+c+".tif")

        #print str(r.mean())
        #print r.projection
        #r.plot()
        #plt.show()

        r = r.clip(netherlands)[0]
        rf = r.flatten(order='A')
        print c +'min:{0} max:{1} avg:{2} stdd:{3}'.format(str(r.min()), str(r.max()),str(r.mean()), str(r.std()))
        #np.hist(rf)
        #import seaborn as sns


        #sns.set_style('darkgrid')
        #sns_plot = sns.distplot(rf)
        #sns_plot.savefig('cellstats/'+rastertype+c+'.png')
        #np.
        plt.clf()
        df = r.to_pandas().value
        df.plot.hist(bins=20, color='#0504aa', alpha=0.7, rwidth=0.85)
        #plt.ticklabel_format(useOffset=False, style='plain')
        #plt.hist(x=rf, bins= 40, color='#0504aa', alpha=0.7, rwidth=0.85)
        plt.grid(axis='y', alpha=0.75)
        plt.xlabel(rastertype+c)
        plt.ylabel('Frequency')
        plt.title(rastertype+c)
        #plt.text(23, 45, r'$\mu=15, b=3$')
        #maxfreq = rf.max()
        # Set a clean upper y-axis limit.
        #plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
        #plt.show()
        plt.savefig('cellstats/'+rastertype+c+'.png')
        r.raster[r.raster == 0] =np.nan #This is in order to computre quantiles for non-zero values only
        df = r.stats(netherlands_w, stats='percentile_10 percentile_20 percentile_30 percentile_40 percentile_50 percentile_60 percentile_70 percentile_80 percentile_90 percentile_100')
        df.insert(0, 'cat', rastertype+c)
        stats.append(df)
        #print rasterstats.zonal_stats(netherlands,  r, stats="percentile_10 percentile_50 percentile_80")

       # print " median = "+str(df.quantile(0.5))
    res = pandas.concat(stats).set_index('cat').drop(['GeoRaster','Id'],1)
    print res
    res.to_csv('cellstats/'+rastertype+'percentiles.csv')


def generateCBSStats():
    inFeatures = 'stats\\buurt_2014.shp'
    attributes = ['BEV_DICHTH','P_65_EO_JR','P_EENP_HH', 'P_HH_Z_K','P_N_W_AL']
    arcpy.env.overwriteOutput = True #Such that files can be overwritten
    arcpy.env.outputCoordinateSystem = inFeatures
    print attributes
    #arcpy.env.workspace = r"C:\Users\schei008\Documents\github\gps\cellstats" #Setting the workspace
    arcpy.env.extent =inFeatures
    # Set local variables
    for a in attributes:
        valField = a
        rastname = (("r"+a)[:13])+".tif"
        outRaster = os.path.join(r"C:\Users\schei008\Documents\github\gps\cellstats",rastname)
        assignmentType = "CELL_CENTER"
        cellSize = 100
        print valField
        # Execute PolygonToRaster
        arcpy.PolygonToRaster_conversion(inFeatures, valField, outRaster, assignmentType, "", cellSize)
        try:
            arcpy.RasterToOtherFormat_conversion(outRaster,arcpy.env.workspace,"TIFF")
            #arcpy.Delete_management(outRaster)
        except:
            print "Raster To Other Format exsample failed."
            print arcpy.GetMessages()
    #return outRaster+'.tif'




def main():
    bbgcategories = {'rsda': [20],'parks':[40,42], 'agric':[51], 'sports':[41], 'recr': [43,44], 'frst':[60], 'ntr':[61,62], 'wtr':[70,71,72,73,74,75,76,77,78,80,81,82,83], 'trfc':[11]}
    """meaning of BBG 2012 landuse categories: rsda (residential areas), agric (agriculture areas), sports (sports areas) recr (recreation areas),frst (forest areas), ntr (nature areas), wtr (water areas), trfc (traffic areas)"""
    attributes = {'BEV_DICHTH':None,'P_65_EO_JR':None,'P_EENP_HH':None, 'P_HH_Z_K':None,'P_N_W_AL':None}
    bbgcategories = {'rsda': [20]} #residential areas
    start = time.time()
    #for name,v in bbgcategories.items():
    #    generatedistraster(v,name)
    #    generateCoverageRaster(name)
    #store distances to these landuse areas for points in some track
    enrichtrack('GPS.csv',bbgcategories,'dist2')
    #enrichCBS('GPS.csv')
    #generateCBSStats()
    #generateCellStats(attributes, rastertype='r', sourcefolder = 'cellstats')
    #generateCellStats(bbgcategories, rastertype='covOf')
    #generateCellStats(bbgcategories, rastertype='dist2')

    end = time.time()
    print("Duration = "+str(end - start))







if __name__ == '__main__':
    main()
