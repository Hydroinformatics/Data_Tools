# -*- coding: utf-8 -*-

import shapely, os
from shapely.geometry import shape
from osgeo import ogr
import numpy as np


def centeroidnp(arr):
    
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    
    return sum_x/length, sum_y/length

#%%

# Save extent to a new Shapefile
#outFileName = '\wr_v_pou_public_proj_REGION_IRRIGATION_ConvexHull.shp'
outFileName = '\states_convexhull.shp'
outShapefile = r'Z:\Projects\INFEWS\Modeling\GIS_Data\Data\wr_state' + outFileName
outDriver = ogr.GetDriverByName("ESRI Shapefile")

# Remove output shapefile if it already exists
if os.path.exists(outShapefile):
    outDriver.DeleteDataSource(outShapefile)

# Create the output shapefile
outDataSource = outDriver.CreateDataSource(outShapefile)
#outLayer = outDataSource.CreateLayer("wr_v_pou_public_proj_REGION_IRRIGATION_ConvexHull", geom_type=ogr.wkbPolygon)
outLayer = outDataSource.CreateLayer("states_convexhull", geom_type=ogr.wkbPolygon)

# Add an ID field
idField = ogr.FieldDefn("id", ogr.OFTInteger)
outLayer.CreateField(idField)

idField = ogr.FieldDefn("obj_wrid", ogr.OFTInteger)
outLayer.CreateField(idField)

idField = ogr.FieldDefn("snp_id", ogr.OFTInteger)
outLayer.CreateField(idField)

#%%

path_shp = r'Z:\Projects\INFEWS\Modeling\GIS_Data\Data\wr_state\wr_v_pou_public_proj_REGION_IRRIGATION.shp'

driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(path_shp, 0)
layer = dataSource.GetLayer()

ft_idc = 0
poly_area = dict()
pts_len =[]

for feature in layer:
    geom = feature.geometry()
    points = feature.GetGeometryRef()
    ft_id = feature.GetField("OBJECTID")
    snp_id = feature.GetField("snp_id")
    
    # Create the feature and set values
    featureDefn = outLayer.GetLayerDefn()
    featureo = ogr.Feature(featureDefn)
    
    temp_area = []
    diam_bool = 0
    for pt in points:
#        if pt.GetArea() > 220.0 and pt.GetArea() < 225.0:
#            diam_bool = 1
        temp_area.append(pt.GetArea())
        
    if np.mean(temp_area) > 220.0 and np.mean(temp_area) < 225.0:
        diam_bool = 1
        
    pts = []
    area = []
    cpoly = 1
    if diam_bool == 1:
        for pt in points:
    
            featureo.SetField("id", ft_idc)
            featureo.SetField("obj_wrid", ft_id)
            featureo.SetField("snp_id", snp_id)
        
            #diam_bool = 1
            bufferDistance = 500
            
            if not pt.GetPoints():
                ptsc = pt.GetGeometryRef(0)
                ptsc = np.asarray(ptsc.GetPoints())
            else:
                ptsc = np.asarray(pt.GetPoints())
            
            
            xc, yc = centeroidnp(ptsc)
            wkt = 'POINT (' + str(xc) + ' ' + str(yc) + ')'
            ptsc = ogr.CreateGeometryFromWkt(wkt)
            poly = ptsc.Buffer(bufferDistance)
            ##print poly
            ##poly = poly.GetGeometryRef(0)
            
            #featureo.SetGeometry(poly)
            
            ##pts.append(poly.GetPoints())
            ##area.append(pt.GetArea())
            #outLayer.CreateFeature(featureo)
            #ft_idc = ft_idc + 1
            
            
            if cpoly == 1:
                poly2 = poly
                #poly2 = ogr.Geometry(ogr.wkbMultiPolygon)
                cpoly = 0
            else:
                poly2 = poly2.Union(poly)
            
            #poly2.AddGeometry(poly)
                
        
        featureo.SetGeometry(poly2)
        
        #convexhull = poly2.ConvexHull()
        #featureo.SetGeometry(convexhull)
        
        outLayer.CreateFeature(featureo)
        ft_idc = ft_idc + 1
            
    else:
        featureo.SetField("id", ft_idc)
        featureo.SetField("obj_wrid", ft_id)
        featureo.SetField("snp_id", snp_id)
        
        pts.append(pt.GetPoints())
        #area.append(pt.GetArea())
        featureo.SetGeometry(points)
        ft_idc = ft_idc + 1
            
            
    outLayer.CreateFeature(featureo)
        
        
    pts_len.append((ft_id,len(pts)))
    
    
#    if diam_bool == -1:
#        convexhull = points.ConvexHull()
#        featureo.SetGeometry(convexhull)
#    else:
#        featureo.SetGeometry(points)
        
    
    #poly_area[ft_id] = area
    poly_area[ft_id] = diam_bool
    
#    print len(pts)
#    #print feature.GetField("use_code_d")

#for i in range(layer.GetFeatureCount()):
#    #feature = layer.GetFeature(i)
#    #name = feature.GetField("NAME")
#    #geom = feature.GetGeometryRef()
#    #ring = geom.GetGeometryRef(0)
#    geom = feature.geometry()
#    
#    #pointx = geom.GetY()
#    points = geom.GetPoints()
#    print i, points
#    
layer.ResetReading()

#feature = shape.GetFeature(0)
#first = feature.ExportToJson()
#
#shp_geom = shape(first['geometry']) # or shp_geom = shape(first) with PyShp)
#print shp_geom


## Calculate convex hull
#convexhull = geomcol.ConvexHull()

feature = None

# Save and close DataSource
inDataSource = None
outDataSource = None
