# -*- coding: utf-8 -*-

from osgeo import ogr
from shapely.geometry import shape

import gdal, os, shutil, zipfile, json, shapely
import matplotlib.pyplot as plt
import numpy as np

import rasterio
from rasterio.mask import mask

def Get_No_Crop_Ids():
    # Land use IDS of non-crops
    # Use to convert these land uses to background
    no_crop_ids = [0,60,63,64,65,81,82,83,87,88,92,111,112,121,
                   122,123,124,131,141,142,143,152,176,190,195]
               
    # Only consider major Crops
    no_crop_ids =[no_crop_ids,range(2,7),10,11,13,22,range(25,31),range(32,36),38,39,
                  41,42,45,46,47,48,50,51,52,54,55,56,58,60,range(67,71),72,74,
                  75,76,77,204,range(206,215),range(216,228),range(229,251),254]
    
    temp = []
    for i in no_crop_ids:
        if np.size(i) > 1:
            for ii in i:
                temp.append(ii)
        else:
            temp.append(i)
    
    no_crop_ids = temp
    
    no_crop_ids.remove(5)
    no_crop_ids.remove(27)
    
    return no_crop_ids

#%%
def records(file):  
    # generator 
    reader = ogr.Open(file)
    layer = reader.GetLayer(0)
    for i in range(layer.GetFeatureCount()):
        feature = layer.GetFeature(i)
        yield json.loads(feature.ExportToJson())

#%%
def Raster_row_boundaries(raster):
    
    row_sum_index = np.where(np.sum(raster,axis=1) != 0)
    top_row = row_sum_index[0][0]
    last_row = row_sum_index[0][-1]
    
    if top_row < 0:
        top_row = 0
    if last_row < 0:
        last_row = 0
        
    return top_row, last_row

def Raster_col_boundaries(raster):
    
    col_sum_index = np.where(np.sum(raster,axis=0) != 0)
    left_col = col_sum_index[0][0]
    right_col = col_sum_index[0][-1]
    
    if left_col < 0:
        left_col = 0
    if right_col < 0:
        right_col = 0
        
    return left_col, right_col

#%%
def Read_Raster(raster_file, raster_NoData = -999):
    
    raster_ds = gdal.Open(raster_file, gdal.GA_ReadOnly)
    raster_NoData = raster_ds.GetRasterBand(1).GetNoDataValue()
    raster = raster_ds.GetRasterBand(1).ReadAsArray()     
    
    return raster, raster_NoData, raster_ds

def GetPixelCentroids(raster):
    upx, xres, xskew, upy, yskew, yres = raster.GetGeoTransform()
    cols = raster.RasterXSize
    rows = raster.RasterYSize
    
    centroids_X = np.ones((1,cols))*-999
    centroids_Y = np.ones((rows,1))*-999
    
    for j in range(rows):
        if j == 0:
            centroids_Y[j,0] = upy - np.abs(yres)/2.0
        else:
            centroids_Y[j,0] = centroids_Y[j-1,0] - np.abs(yres) 
            
    for i in range(cols):
        if i == 0:
            centroids_X[0,i] = upx + np.abs(xres)/2.0 
        else:
            centroids_X[0,i] = centroids_X[0,i-1] + np.abs(xres)
    
    return centroids_X, centroids_Y

#%%

cdl_path = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\Umatilla_Input_data\CDL_Data\CDL_Data_Filter'

years = range(2008,2018)

fnames = os.listdir(cdl_path)
findex = []
c = 0

for fname in fnames:
    if '.tif.' not in fname and '.tfw' not in fname and '.aux.xml' not in fname and '_mfv2.tif' in fname:
        findex.append((c,int(fname[4:8])))
    c += 1    

findex = np.asarray(findex)
#for year in years: 
#            
#    print ('Zonal_Stats of ' + fnames[findex[findex[:,1]==year,0][0]])
#            cdl_raster = files['cdl_path'] + '\\' + fnames[findex[findex[:,1]==year,0][0]]            

#%%

no_crop_ids = Get_No_Crop_Ids()

cdl_rasters = dict()
for fname in fnames:
    raster_file = cdl_path + '\\' + fname
    region, region_NoData, region_obj = Read_Raster(raster_file)
    cdl_rasters[int(fname[4:8])] = region
    
    
raster_file = cdl_path + '\\' + fnames[0]
region, region_NoData, region_obj = Read_Raster(raster_file)

centroids_X, centroids_Y = GetPixelCentroids(region_obj)

top_row, last_row = Raster_row_boundaries(region)
left_col, right_col = Raster_col_boundaries(region)
#region = region[top_row:last_row+1,left_col:right_col+1]

#plt.matshow(region)
#
#multipoint = ogr.Geometry(ogr.wkbMultiPoint)
##multipoint = ogr.Geometry(ogr.wkbPoint)
#
##CDL_centroids = np.zeros((len(region.flatten()),3))
#
#n=0
#for i in range(top_row,last_row+1):
#    for j in range(left_col,right_col+1):
##        CDL_centroids[n,0] = n 
##        CDL_centroids[n,1] = centroids_X[0,j]
##        CDL_centroids[n,2] = centroids_Y[i,0] 
#        
#        point1 = ogr.Geometry(ogr.wkbPoint)
#        point1.AddPoint(centroids_X[0,j], centroids_Y[i,0])
#        multipoint.AddGeometry(point1)
#        
#        #multipoint.AddPoint(centroids_X[0,j], centroids_Y[i,0])
#        
#        n = n + 1

#np.savetxt(cdl_path + '/' + 'CDL_Raster_Centroids.csv', CDL_centroids, delimiter=",") 
      
#points  = [pt for pt in records(CDL_centroids)]
##%%      


#path_shp = r'Z:\Projects\INFEWS\Modeling\GIS_Data\Data\wr_state\wr_v_pou_public_proj_REGION_IRRIGATION_ConvexHullV2.shp'
#wr_shp = records(path_shp)

csv_path = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\Umatilla_Input_data\CDL_Data'
csv_file = csv_path + '\OWRD_CDL_Match_CropIDs.csv'
csv_file2 = csv_path + '\OWRD_CDL_Match_CropCounts.csv'

filein = open(csv_file,'w')
filein2 = open(csv_file2,'w')

path_shp = r'Z:\Projects\INFEWS\Modeling\GIS_Data\Data\wr_state\wr_v_pou_public_proj_REGION_IRRIGATION_ConvexHullV2.shp'

driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(path_shp, 0)
layer = dataSource.GetLayer()

poly_area = dict()
pts_len =[]

CDL_dict = dict()
for feature in layer:
    geom = feature.geometry()
    points = feature.GetGeometryRef()
    
    ftid = feature.GetField("id")
    wrid = feature.GetField("wrid")
    print wrid
    
    if wrid not in CDL_dict.keys():
        CDL_dict[wrid] = dict()
        
        if ftid not in CDL_dict[wrid].keys():
            CDL_dict[wrid][ftid] = dict()
    
            points_json = [json.loads(feature.GetGeometryRef().ExportToJson())]
            
            if str(points_json[0]['type']) == 'MultiPolygon':
                pts = []
                for pp in range(len(points_json[0]['coordinates'])):
                    for ppp in range(len(points_json[0]['coordinates'][pp])):
                        for xpt, ypt in points_json[0]['coordinates'][pp][ppp]:
                            pts.append((xpt,ypt))
                pts = np.asarray(pts)   
                
            else:
                pts = np.asarray(points_json[0]['coordinates'][0])
        #for multi in wr_shp:
            #wrid = multi['properties']['wrid']
            
            xmin = np.min(pts[:,0])
            xmax = np.max(pts[:,0])
            
            ymin = np.min(pts[:,1])
            ymax = np.max(pts[:,1])
            
            
            xrows_ids = np.where((centroids_X[0,:] >= xmin) & (centroids_X[0,:] <= xmax))[0]
            ycols_ids = np.where((centroids_Y[:,0] >= ymin) & (centroids_Y[:,0] <= ymax))[0]
            
            
            #if wrid not in CDL_dict.keys():
            #    CDL_dict[wrid] = dict()
                
            ##ftid = multi['properties']['id']
            #if ftid not in CDL_dict[wrid].keys():
            #     CDL_dict[wrid][ftid] = dict()
            
            temp_points = []
            
        #    for i in range(top_row,last_row+1):
        #        for j in range(left_col,right_col+1):
            for i in ycols_ids:
                for j in xrows_ids:
                    point = ogr.Geometry(ogr.wkbPoint)
                    point.AddPoint(centroids_X[0,j], centroids_Y[i,0])
                    
                    if point.Within(points):
                        #print i, shape(points[i]['geometry'])
                        temp_points.append((i,j))
                        #temp_data.append(region[i,j])
                        
                        temp_data = []
            temp_points = np.asarray(temp_points)
            
            
            for yr in years:
                temp_data = []
                #CDL_dict[wrid][ftid][yr] = dict()
                for i,j in temp_points:
                    #if cdl_rasters[yr][i,j] not in no_crop_ids: # Eliminate No Crops
                    temp_data.append(cdl_rasters[yr][i,j])
                
                atxt = str(wrid) + ',' + str(ftid) + ',' + str(yr) + ','
                atxt2 = atxt
                if len(temp_data) > 0:
                    temp_data, temp_count = np.unique(np.asarray(temp_data),return_counts=True)
                    for i in range(len(temp_data)):
                        atxt = atxt + str(temp_data[i]) + ','
                        atxt2 = atxt2 + str(temp_count[i]) + ','
                
                filein.write(atxt + '\n')
                filein2.write(atxt2 + '\n')
                #CDL_dict[wrid][ftid][yr]['CDL'] = temp_data
                #CDL_dict[wrid][ftid][yr]['count'] = temp_count
        
filein.close()
filein2.close()

#path_shp = r'Z:\Projects\INFEWS\Modeling\GIS_Data\Data\wr_state\wr_v_pou_public_proj_REGION_IRRIGATION_ConvexHullV2.shp'
#
#driver = ogr.GetDriverByName("ESRI Shapefile")
#dataSource = driver.Open(path_shp, 0)
#layer = dataSource.GetLayer()
#
#ft_idc = 0
#poly_area = dict()
#pts_len =[]
#
#for feature in layer:
#    geom = feature.geometry()
#    points = feature.GetGeometryRef()
#    
#    ft_id = feature.GetField("id")
#    wr_id = feature.GetField("wrid")
#    
#    
#    
##    pts = []
##    for pt in points:
##        if not pt.GetPoints():
##            ptsc = pt.GetGeometryRef(0)
##            ptsc = np.asarray(ptsc.GetPoints())
##        else:
##            ptsc = np.asarray(pt.GetPoints())
##        
##        pts = []
##        for x,y in ptsc:
##            pts.append((x,y))
#    
#    #points = ["{geometry':" + feature.GetGeometryRef().ExportToJson() +"}" ]
#    
#    
#    points_json = [json.loads(feature.GetGeometryRef().ExportToJson())]
##    points_json[0]['coordinates'] = pts
#    
#    #points = ogr.CreateGeometryFromJson(points)
#    raster_file = cdl_path + '\\' + fnames[0]
#    #%%
#    with rasterio.open(raster_file) as src:
#        out_image, out_transform = mask(src, points_json, crop=True)
#        
#        plt.matshow(out_image[0,:,:])
#        x_min = out_transform[0]
#        y_max = out_transform[3]
#        #x_max = x_min + geo_transform[1] * Image.RasterXSize
#        #y_min = y_max + geo_transform[5] * Image.RasterYSize
#        x_max = out_transform[1]
#        y_min = out_transform[2]
#
#        x_res = int(np.round((x_max - x_min)/pixel_width))
#        y_res = int(np.round((y_max - y_min)/pixel_width))
#        
##%%
#
##shapefile = [{'geometry': {'coordinates': [[(-4.663611, 51.158333),(-4.669168, 51.159439),(-4.673334, 51.161385),(-4.674445, 51.165276),(-4.67139, 51.185272),(-4.669445, 51.193054),(-4.665556, 51.195),(-4.65889, 51.195),(-4.656389, 51.192215),(-4.646389, 51.164444),(-4.646945, 51.160828),(-4.651668, 51.159439),(-4.663611, 51.158333)]], 'type': 'Polygon'},
## 'id': '1','properties': '','type': 'Feature'}]
##    
##shapes = [feature["geometry"] for feature in shapefile]