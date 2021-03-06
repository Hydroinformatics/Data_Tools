# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 13:18:45 2019

@author: sammy
"""

from osgeo import ogr
import os, json
from shapely.geometry import shape

import gdal, os, shutil, zipfile
import matplotlib.pyplot as plt
import numpy as np


#%%
def records(file):  
    # generator 
    reader = ogr.Open(file)
    layer = reader.GetLayer(0)
    for i in range(layer.GetFeatureCount()):
        feature = layer.GetFeature(i)
        yield json.loads(feature.ExportToJson())

#%%
def UnzipData(path,folderpath):

    if os.path.isdir(folderpath):
        shutil.rmtree(folderpath)
    
    #print folderpath
    os.makedirs(folderpath)
    #print path
    with zipfile.ZipFile(path, "r") as z:
        z.extractall(folderpath)
    
  
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

def Read_Raster(raster_file, raster_NoData = -999):
    
    raster_ds = gdal.Open(raster_file, gdal.GA_ReadOnly)
    raster_NoData = raster_ds.GetRasterBand(1).GetNoDataValue()
    raster = raster_ds.GetRasterBand(1).ReadAsArray()     
    
    return raster, raster_NoData, raster_ds


def ReadBilFile(bil):
    
    gdal.GetDriverByName('EHdr').Register()
    img = gdal.Open(bil)
    band = img.GetRasterBand(1)
    data = band.ReadAsArray()
    return data

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
prism_centroids = r"C:\Users\sammy\Downloads\PRISM_RecentYears\PRISM_Stations.shp"
points  = [pt for pt in records(prism_centroids)]

watershed_shp = r"C:\Users\sammy\Downloads\PRISM_RecentYears\Watersheds_proj.shp"
watershed_shp = records(watershed_shp)

PRISM_dict = dict()



for multi in watershed_shp:
    temp_points = []
    for i, pt in enumerate(points):
         point = shape(pt['geometry'])
         if point.within(shape(multi['geometry'])):
              #print i, shape(points[i]['geometry'])
              temp_points.append(i)
              
    PRISM_dict[multi['id']] = temp_points
          

main_path = r'C:\Users\sammy\Downloads\PRISM_RecentYears'
pathunzip = r'C:\Users\sammy\Downloads\PRISM_RecentYears\temp'
begin_year = 1981
end_year = 2019
system_boundary = 'System_Boundary_16km.tif'

#%%
path = main_path + '/' + system_boundary
region, region_NoData, region_obj = Read_Raster(path)

centroids_X, centroids_Y = GetPixelCentroids(region_obj)

top_row, last_row = Raster_row_boundaries(region)
left_col, right_col = Raster_col_boundaries(region)
region = region[top_row:last_row+1,left_col:right_col+1]

PRISM_centroids = np.zeros((len(region.flatten()),3))

n=0
for i in range(top_row,last_row+1):
    for j in range(left_col,right_col+1):
        PRISM_centroids[n,0] = n 
        PRISM_centroids[n,1] = centroids_X[0,j]
        PRISM_centroids[n,2] = centroids_Y[i,0] 
        
        n = n + 1
    
#np.savetxt(main_path + '/' + 'PRISM_Raster_Centroids.csv', PRISM_centroids, delimiter=",")

#%%
nrows = 0
for year in range(begin_year,end_year+1):
    path = main_path + '/' + str(year)
    fnames = os.listdir(path)
    nrows = nrows + len(fnames)
    
ncols = len(region.flatten()) + 3 # Year, Mon, Day

#%%
PRISM_Data = np.ones((nrows,ncols))*-999
nrow = 0
for year in range(begin_year,end_year+1):

    path = main_path + '/' + str(year)
    fnames = os.listdir(path)
    
    for name in fnames:
        dash_index = [i for i in range(len(name)) if name[i] == '_']
        datestr = name[dash_index[-2]+1:dash_index[-1]]
        PRISM_Data[nrow,0] = int(datestr[0:4])
        PRISM_Data[nrow,1] = int(datestr[4:6])
        PRISM_Data[nrow,2] = int(datestr[6:])
        
        pathzip = path + '/' + name
        UnzipData(pathzip,pathunzip)
        
        temp_fnames = os.listdir(pathunzip)
        main_bil_file = [tname for tname in temp_fnames if '.bil' in tname[-4:]]
        print main_bil_file[0]
        data = ReadBilFile(pathunzip + '/' + main_bil_file[0])
        data = data[top_row:last_row+1,left_col:right_col+1]
        data = np.transpose(data.flatten())
        PRISM_Data[nrow,3:3+len(data)] = data
        nrow = nrow + 1

#np.savetxt(main_path + '/' + 'PRISM_PPT_Data.csv', PRISM_Data, delimiter=",")
        
#%%
yearly_total = dict()
for wid in PRISM_dict.keys():
    yearly_total[wid] = dict()
    ny = 1
    for year in range(begin_year, end_year+1):
        temp_total_precip = []
        for stid in PRISM_dict[wid]:
            temp_total_precip.append(np.sum(PRISM_Data[np.where(PRISM_Data[:,0] == year),stid+3]))
            
        yearly_total[wid][ny] = np.mean(temp_total_precip)
        ny = ny + 1
            
        
with open('PRISM_Mean_Yearly_Totals.json', 'w') as fp:
    json.dump(yearly_total, fp)      
