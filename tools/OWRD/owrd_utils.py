# -*- coding: utf-8 -*-

import sys, datetime, gdal, os, shutil, zipfile, json, shapely
from statistics import mode
from scipy import stats, ndimage
qswat_utils_path = r'C:\Users\sammy\Documents\GitHub\InterACTWEL\src'
sys.path.append(qswat_utils_path)

from qswat import QSWAT_utils

from osgeo import ogr
from shapely.geometry import Point, shape

import matplotlib.pyplot as plt
import numpy as np

import rasterio
from rasterio.mask import mask


#%%
def records(file):  
    # generator 
    reader = ogr.Open(file)
    layer = reader.GetLayer(0)
    for i in range(layer.GetFeatureCount()):
        feature = layer.GetFeature(i)
        yield json.loads(feature.ExportToJson())

#%%
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
def centeroidnp(arr):
    
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    
    return sum_x/length, sum_y/length

#%% Get Basins for PODs
        
def GetSWATBasinForPODs(root_path, subbasins_file, pod_file):
    sub_pod_dict = dict()
    points  = [pt for pt in records(root_path + '/' + pod_file)]
    subbasin =[sb for sb in records(root_path + '/' + subbasins_file)]
    for i, pt in enumerate(points):
        point = shape(pt['geometry'])
        if int(pt['properties']['snp_id']) not in sub_pod_dict.keys():
            sub_pod_dict[int(pt['properties']['snp_id'])] = dict()
        sub_pod_dict[int(pt['properties']['snp_id'])][int(pt['properties']['pod_use_id'])] = 'nan'
        
        for ii in range(len(subbasin)):
            if point.within(shape(subbasin[ii]['geometry'])):
                sub_pod_dict[int(pt['properties']['snp_id'])][int(pt['properties']['pod_use_id'])] = int(subbasin[ii]['properties']['Subbasin'])
                #break
        
    return sub_pod_dict


def GetSWATBasinForPOUs(root_path, subbasins_file, pou_file):
    sub_pou_dict = dict()
    points  = [pt for pt in records(root_path + '/' + pou_file)]
    subbasin =[sb for sb in records(root_path + '/' + subbasins_file)]
    
    for i, pt in enumerate(points):
        point = shape(pt['geometry'])
        if int(pt['properties']['snp_id']) not in sub_pou_dict.keys():
            sub_pou_dict[int(pt['properties']['snp_id'])] = dict()
        sub_pou_dict[int(pt['properties']['snp_id'])][int(pt['properties']['pou_use_id'])] = []
        
        for ii in range(len(subbasin)):
            if point.within(shape(subbasin[ii]['geometry'])) or point.intersects(shape(subbasin[ii]['geometry'])):
                sub_pou_dict[int(pt['properties']['snp_id'])][int(pt['properties']['pou_use_id'])].append(int(subbasin[ii]['properties']['Subbasin']))
                #break
        
    return sub_pou_dict

#%%
def Eliminate_Duplicate_Records(root_path, input_file_name, output_file_name):    
    
    #Open input shapefile and get spatial reference of it
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(root_path + '/' + input_file_name, 0)
    wr_layer = dataSource.GetLayer()
    
    spatialRef = wr_layer.GetSpatialRef()
    
    #Define output shapefile
    outShapefile = root_path + '/' + output_file_name + '.shp'
    outDriver = ogr.GetDriverByName("ESRI Shapefile")

    # Remove output shapefile if it already exists
    if os.path.exists(outShapefile):
        outDriver.DeleteDataSource(outShapefile)
    
    proj_outFileName = output_file_name + '.prj'
    with open(root_path + '/' + proj_outFileName, 'w') as prj_file:
        prj_file.write(str(spatialRef.ExportToWkt()))
    prj_file.close()

    # Create the output shapefile
    outDataSource = outDriver.CreateDataSource(outShapefile)
    outLayer = outDataSource.CreateLayer("wr_v_pou_public", geom_type = ogr.wkbPolygon)

    ldefn = wr_layer.GetLayerDefn()
    field_names = []
    for n in range(ldefn.GetFieldCount()):
        fdefn = ldefn.GetFieldDefn(n)
        ftype = fdefn.GetFieldTypeName(fdefn.GetType())
        
        # Add an ID field
        if 'real' in ftype.lower():
            idField = ogr.FieldDefn(fdefn.name, ogr.OFTReal)
            
        elif 'integer' in ftype.lower():
            idField = ogr.FieldDefn(fdefn.name, ogr.OFTInteger)
            
        elif 'string' in ftype.lower():
            idField = ogr.FieldDefn(fdefn.name, ogr.OFTString)
            
        elif 'date' in ftype.lower():
            idField = ogr.FieldDefn(fdefn.name, ogr.OFTDate)
        
        outLayer.CreateField(idField)
        field_names.append(fdefn.name)
    
    # Go over all records in input shapefile and identify duplicates using
    # similarity of all fields. In the special case in which only the
    # shape_leng and shape_area are the only different fields, then a union
    # of the polygons associated with the corresponding WR is used.
    
    feature_dict = dict() 
    remove_features = []
    for feature in wr_layer:     
        match_bool = 0
        snp_id = feature.GetField("snp_id")
        obj_id = feature.GetField("OBJECTID")
        
        if snp_id not in feature_dict.keys():
            feature_dict[snp_id] = dict()
            feature_dict[snp_id][obj_id] = dict()
            
        else:
            for objs in feature_dict[snp_id].keys():
                feat_sim_count = 0
                non_equal_feat = []
                for fname in feature_dict[snp_id][objs].keys():
                    field_val = feature.GetField(fname)
                    
                    # Change use code 'SUPPLEMENTAL IRRIGATION' to 'IRRIGATION'
                    # (Often, same polygons and same records, only use code is different)
                    if fname == 'use_code_d' and field_val == 'SUPPLEMENTAL IRRIGATION':
                        field_val = 'IRRIGATION'
                    if fname == 'use_code' and field_val == 'IS':
                        field_val = 'IR'
                        
                    if field_val == feature_dict[snp_id][objs][fname]: 
                        feat_sim_count = feat_sim_count + 1
                    else:
                        non_equal_feat.append(fname)
                
                # Check if all fields have the same value. If yes, remove record.
                if feat_sim_count == len(feature_dict[snp_id][objs].keys()):
                    match_bool = 1
                    remove_features.append([1,snp_id,obj_id,objs])
                    break
        
            if len(non_equal_feat) == 2:
                if non_equal_feat[0].lower() == 'shape_leng' and non_equal_feat[1].lower() == 'shape_area':
                    remove_features.append([2,snp_id,obj_id,objs])
                    match_bool = 1
            
        if match_bool == 0:
            if obj_id not in feature_dict[snp_id].keys():
                feature_dict[snp_id][obj_id] = dict()
                
            for fname in field_names:
                field_val = feature.GetField(fname)   
                 
                if fname == 'use_code_d' and field_val == 'SUPPLEMENTAL IRRIGATION':
                    field_val = 'IRRIGATION'
                if fname == 'use_code' and field_val == 'IS':
                    field_val = 'IR'
                    
                if fname != 'OBJECTID' and fname != 'snp_id':
                    feature_dict[snp_id][obj_id][fname] = field_val
 
################ Re-write layer with only unique records ##################### 
    remove_features_arr = np.array(remove_features)
    union_features = remove_features_arr[np.where(remove_features_arr[:,0]==2)[0],:]
    remove_features_arr = remove_features_arr[np.where(remove_features_arr[:,0]==1)[0],:]
    
    wr_layer.ResetReading()
    match_bool = 0
    obj_id_sum = 0
    for feature in wr_layer:
        points = feature.GetGeometryRef()
        obj_id = feature.GetField("OBJECTID")
        
        if obj_id not in remove_features_arr[:,2]:
            if match_bool == 0:
                featureDefn = outLayer.GetLayerDefn()
                featureo = ogr.Feature(featureDefn)
            
                for fname in field_names:
                    field_val = feature.GetField(fname)
                    
                    if fname == 'use_code_d' and field_val == 'SUPPLEMENTAL IRRIGATION':
                        field_val = 'IRRIGATION'
                    if fname == 'use_code' and field_val == 'IS':
                        field_val = 'IR'
                    
                    featureo.SetField(fname, field_val)
        
            if obj_id in union_features[:,3]:
                obj_id_sum = len(np.where(obj_id == union_features[:,3])[0])+2
                poly = ogr.CreateGeometryFromWkt(points.ExportToWkt())
                match_bool = 1
            
            if match_bool == 1:
                match_bool = match_bool + 1
                
            elif match_bool > 1 and match_bool <= obj_id_sum:
                points = ogr.CreateGeometryFromWkt(points.ExportToWkt())
                poly = poly.Union(points)
                match_bool = match_bool + 1
                points = poly
            
            if match_bool == obj_id_sum:
                match_bool = 0
                obj_id_sum = 0
                
            if match_bool == 0:
                featureo.SetGeometry(points)
                outLayer.CreateFeature(featureo)
            
    outDataSource = None
    
    return feature_dict, remove_features


#%%
def GetWR_POU_Data(file_name, root_path):
    
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(root_path + '/' + file_name, 0)
    wr_layer = dataSource.GetLayer()
    
    wr_pou_dict = dict()
    for feature in wr_layer:
        snap_id = int(feature.GetField("snp_id"))
        if snap_id not in wr_pou_dict.keys():
            wr_pou_dict[snap_id]= dict()
        
        temp_len = int(feature.GetField("pou_use_id"))
        wr_pou_dict[snap_id][temp_len] = dict()
        
        wr_pou_dict[snap_id][temp_len]['Source'] = feature.GetField("wr_type").strip()
        
        linesplit = feature.GetField("use_code_d").strip()
        wr_pou_dict[snap_id][temp_len]['Unmod_Use'] = linesplit
        if 'IRRIGATION LIVESTOCK' in linesplit.strip():
            linesplit = linesplit.replace('N L','N, L')
        
        if ' AND ' in linesplit.strip() or ', ' in linesplit.strip():
            temp_use = linesplit.strip().split(', ')
            and_split =  temp_use[len(temp_use)-1].split(' AND ')
            temp_use[len(temp_use)-1] = and_split[0]
            temp_use.append(and_split[1])
            wr_pou_dict[snap_id][temp_len]['Use'] = temp_use
        else:
            wr_pou_dict[snap_id][temp_len]['Use'] = [linesplit.strip()]
        
        if feature.GetField("priority_d"):
            celltext = feature.GetField("priority_d")
            slashid = [i for i in range(len(celltext)) if celltext[i] is '/' or celltext[i] is ':']
            
            day = int(celltext[slashid[1]+1:slashid[1]+5])
            if len(str(day)) == 1:
                day = '0' + str(day)
            mon = int(celltext[slashid[0]+1:slashid[1]])
            if len(str(mon)) == 1:
                mon = '0' + str(mon)
            year = int(celltext[0:slashid[0]])
            
            #wr_pou_dict[snap_id][temp_len]['priority'] = datetime.datetime(year,mon,day)
            wr_pou_dict[snap_id][temp_len]['priority'] = str(year) + '-' + str(mon) + '-' + str(day)
        else:
            wr_pou_dict[snap_id][temp_len]['priority'] = '-'
         
        if feature.GetField("wris_acres") > 0:
            wr_pou_dict[snap_id][temp_len]['acre_feet'] = float(feature.GetField("wris_acres"))
        else:
            wr_pou_dict[snap_id][temp_len]['acre_feet'] = float('nan')

    return wr_pou_dict


#%%
def GetWR_POD_Data(file_name, root_path, pou_keys, sub_pod_dict):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(root_path + '/' + file_name, 0)
    pod_layer = dataSource.GetLayer()
    wr_pod_dict = dict()
    
    for feature in pod_layer:
        snp_id = feature.GetField("snp_id")        
        if int(snp_id) in pou_keys:
            if int(snp_id) not in wr_pod_dict.keys():
                wr_pod_dict[int(snp_id)] = dict()
            
            pod_id = feature.GetField("pod_use_id")
            wr_pod_dict[int(snp_id)][int(pod_id)] = dict()
            wr_pod_dict[int(snp_id)][int(pod_id)]['pod_loc_id'] = int(feature.GetField("pod_locati"))
            wr_pod_dict[int(snp_id)][int(pod_id)]['pod_use_id'] = int(feature.GetField("pod_use_id"))
            wr_pod_dict[int(snp_id)][int(pod_id)]['wr_type'] = feature.GetField("wr_type").strip()
            wr_pod_dict[int(snp_id)][int(pod_id)]['wr_type_2'] = feature.GetField("wr_type").strip()
            wr_pod_dict[int(snp_id)][int(pod_id)]['source_type'] = feature.GetField("source_typ").strip()
            
            use_code = feature.GetField("use_code_d").strip()
            wr_pod_dict[int(snp_id)][int(pod_id)]['Unmod_Use'] = use_code
            
            if use_code.strip() == 'PRIMARY AND SUPPLEMENTAL IRRIGATION':
                use_code = 'IRRIGATION'
                
            if use_code.strip() == 'IRRIGATION AND DOMESTIC':
                use_code = 'IRRIGATION'
            
            if use_code.strip() == 'DOMESTIC EXPANDED':
                use_code = 'DOMESTIC'
            
            wr_pod_dict[int(snp_id)][int(pod_id)]['use_code'] = use_code
            
            if feature.GetField("priority_d"):
                celltext = feature.GetField("priority_d")
                
                slashid = [i for i in range(len(celltext)) if celltext[i] is '/' or celltext[i] is ':']
                day = int(celltext[slashid[1]+1:slashid[1]+5])
                if len(str(day)) == 1:
                    day = '0' + str(day)
    
                mon = int(celltext[slashid[0]+1:slashid[1]])
                if len(str(mon)) == 1:
                    mon = '0' + str(mon)
                
                year = int(celltext[0:slashid[0]])
                #wr_pod_dict[snp_id][int(pod_id)]['priority'] = datetime.datetime(year,mon,day)
                wr_pod_dict[snp_id][int(pod_id)]['priority'] = str(year) + '-' + str(mon) + '-' + str(day)
            else:
                wr_pod_dict[snp_id][int(pod_id)]['priority'] = '-'
            
            #if float(feature.GetField("duty")) != -999 and float(feature.GetField("duty")) != 0:
            if float(feature.GetField("duty")) != 0:
                wr_pod_dict[int(snp_id)][int(pod_id)]['duty'] = float(feature.GetField("duty"))
            #else:
            #    wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['duty'] = float('nan')
                
            #if float(feature.GetField("rate_cfs")) != -999:
            if float(feature.GetField("rate_cfs")) != 0:
                wr_pod_dict[int(snp_id)][int(pod_id)]['rate'] = float(feature.GetField("rate_cfs"))
            else:
                wr_pod_dict[int(snp_id)][int(pod_id)]['rate'] = float('nan')
                
            #if float(feature.GetField("max_rate_c")) != -999:
            if float(feature.GetField("max_rate_c")) != 0:  
                wr_pod_dict[int(snp_id)][int(pod_id)]['max_rate'] = float(feature.GetField("max_rate_c"))
            else:
                wr_pod_dict[int(snp_id)][int(pod_id)]['max_rate'] = '-'
            
#            if wr_pod_dict[int(snp_id)][int(pod_id)]['max_rate'] == 0:
#                wr_pod_dict[int(snp_id)][int(pod_id)]['max_rate'] = '-'
            
            #if float(feature.GetField("acre_feet_")) != -999: 
            if float(feature.GetField("acre_feet_")) != 0:
                wr_pod_dict[int(snp_id)][int(pod_id)]['acre_feet'] = float(feature.GetField("acre_feet_"))
            #else:
            #    wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['acre_feet'] = float('nan')
                
            wr_pod_dict[int(snp_id)][int(pod_id)]['source'] = feature.GetField("source_typ").strip()
            wr_pod_dict[int(snp_id)][int(pod_id)]['source_name'] = feature.GetField("source").strip()
            
            if feature.GetField("tributary_"):
                wr_pod_dict[int(snp_id)][int(pod_id)]['tributary'] = feature.GetField("tributary_").strip()
            else:
                wr_pod_dict[int(snp_id)][int(pod_id)]['tributary'] = 'Not specified'
                
            
            CR_as_source = ['COLUMBIA R','COLUMBIA RIVER']
            #CR_as_tributary = ['COLUMBIA R','COLUMBIA RIVER','COLUMBIA RIVER BASIN']
            
            if wr_pod_dict[int(snp_id)][int(pod_id)]['wr_type'] == 'SW': #it should be noted that there are a few ST with CR as the tributary
                #if wr_pod_dict[int(snp_id)][int(pod_id)]['source_name'] in  CR_as_source or wr_pod_dict[int(snp_id)][int(pod_id)]['tributary'] in CR_as_tributary:
                if wr_pod_dict[int(snp_id)][int(pod_id)]['source_name'] in  CR_as_source:
                    wr_pod_dict[int(snp_id)][int(pod_id)]['wr_type_2'] = 'CR'
            
            RES_as_source = ['COLD SPRINGS RESERVOIR']
            #CSRes_as_tributary = ['COLUMBIA R','COLUMBIA RIVER','COLUMBIA RIVER BASIN']
            
            if wr_pod_dict[int(snp_id)][int(pod_id)]['wr_type'] == 'SW': #it should be noted that there are a few ST with CR as the tributary
                #if wr_pod_dict[int(snp_id)][int(pod_id)]['source_name'] in  CSRes_as_source or wr_pod_dict[int(snp_id)][int(pod_id)]['tributary'] in CSRes_as_tributary:
                if wr_pod_dict[int(snp_id)][int(pod_id)]['source_name'] in  RES_as_source:
                    #wr_pod_dict[int(snp_id)][int(pod_id)]['wr_type_2'] = 'RES'
                    wr_pod_dict[int(snp_id)][int(pod_id)]['wr_type_2'] = 'CSRes'
                  
                    
            wr_pod_dict[int(snp_id)][int(pod_id)]['begin_mon'] = int(feature.GetField("begin_mont"))
            if int(feature.GetField("begin_day")) == 0:
                wr_pod_dict[int(snp_id)][int(pod_id)]['begin_day'] = 1
            else:
                wr_pod_dict[int(snp_id)][int(pod_id)]['begin_day'] = int(feature.GetField("begin_day"))
                
            
            wr_pod_dict[int(snp_id)][int(pod_id)]['end_mon'] = int(feature.GetField("end_month"))
            if int(feature.GetField("end_day")) == 0:
                wr_pod_dict[int(snp_id)][int(pod_id)]['end_day'] = 1
                
            elif int(feature.GetField("end_day")) == 29 and wr_pod_dict[int(snp_id)][int(pod_id)]['end_mon'] == 2:
                wr_pod_dict[int(snp_id)][int(pod_id)]['end_day'] = 28
                
            else:
                wr_pod_dict[int(snp_id)][int(pod_id)]['end_day'] = int(feature.GetField("end_day"))
            
            if wr_pod_dict[int(snp_id)][int(pod_id)]['end_mon'] == 0:
                wr_pod_dict[int(snp_id)][int(pod_id)]['end_mon'] = 12
            
            if wr_pod_dict[int(snp_id)][int(pod_id)]['end_mon'] == 0 and wr_pod_dict[int(snp_id)][int(pod_id)]['end_day'] == 0:
                wr_pod_dict[int(snp_id)][int(pod_id)]['end_mon'] = 12
                wr_pod_dict[int(snp_id)][int(pod_id)]['end_day'] = 31
            
            wr_pod_dict[int(snp_id)][int(pod_id)]['pump_start'] = datetime.datetime(2010,wr_pod_dict[int(snp_id)][int(pod_id)]['begin_mon'],wr_pod_dict[int(snp_id)][int(pod_id)]['begin_day']).timetuple().tm_yday
            wr_pod_dict[int(snp_id)][int(pod_id)]['pump_end'] = datetime.datetime(2010,wr_pod_dict[int(snp_id)][int(pod_id)]['end_mon'],wr_pod_dict[int(snp_id)][int(pod_id)]['end_day']).timetuple().tm_yday
            
            wr_pod_dict[int(snp_id)][int(pod_id)]['subbasin'] = sub_pod_dict[str(snp_id)][str(pod_id)]
            
    return wr_pod_dict


#%%
def WR_Parser(wr_pou_dict, wr_pod_dict):
    num_pou_recs=0
    num_pou_nan = []
    for i in wr_pou_dict.keys():    
        for ii in wr_pou_dict[i].keys():
            if np.isnan(wr_pou_dict[i][ii]['acre_feet']) == 1:
                num_pou_nan.append((i,ii,0))
            else:
                num_pou_nan.append((i,ii,1))  
                
            num_pou_recs = num_pou_recs + 1        

    wris_per_use = dict()
    empty_acres = []
    non_empty_acres = []
    
    for i in wr_pou_dict.keys():    
        for ii in wr_pou_dict[i].keys():
            for iii in range(len(wr_pou_dict[i][ii]['Use'])):
            
                if wr_pou_dict[i][ii]['Use'][iii] not in wris_per_use.keys():
                    wris_per_use[wr_pou_dict[i][ii]['Use'][iii]] = dict()
                    
                if wr_pou_dict[i][ii]['Source'] not in wris_per_use[wr_pou_dict[i][ii]['Use'][iii]].keys():
                    wris_per_use[wr_pou_dict[i][ii]['Use'][iii]][wr_pou_dict[i][ii]['Source']] = 0
                
                comment = ''
                if np.isnan(wr_pou_dict[i][ii]['acre_feet']) == 1:
                    if i in wr_pod_dict.keys():
                        temp_sum = 0
                        num_pods = 0
                        pods_temp = []
                        pods_time_temp = []
                        find_use_case_match = 0
                        
                        for jjj in wr_pod_dict[i].keys():
                            use_case_match = 0
                            if wr_pou_dict[i][ii]['Use'][iii] == wr_pod_dict[i][jjj]['use_code'] or wr_pou_dict[i][ii]['Unmod_Use'] == wr_pod_dict[i][jjj]['Unmod_Use']:
                                use_case_match = 1
                            if wr_pou_dict[i][ii]['Use'][iii] == 'STORAGE' and wr_pod_dict[i][jjj]['use_code'] == 'IRRIGATION':
                                use_case_match = 1
                                
                            if wr_pou_dict[i][ii]['Use'][iii] == 'DOMESTIC' and wr_pod_dict[i][jjj]['Unmod_Use'] == 'IRRIGATION AND DOMESTIC':
                                use_case_match = 1
                                
                            if use_case_match == 1:
                                find_use_case_match = 1
                                
                            if ('acre_feet' in wr_pod_dict[i][jjj].keys() and wr_pou_dict[i][ii]['Source'] == wr_pod_dict[i][jjj]['wr_type'] and use_case_match == 1):
                                temp_sum = temp_sum + wr_pod_dict[i][jjj]['acre_feet']
                                num_pods = num_pods + 1
                                pods_temp.append(jjj)
                                pods_time_temp.append([wr_pod_dict[i][jjj]['pump_start'],wr_pod_dict[i][jjj]['pump_end']])

                                
                        if temp_sum != 0:                                    
                            wris_per_use[wr_pou_dict[i][ii]['Use'][iii]][wr_pou_dict[i][ii]['Source']] = wris_per_use[wr_pou_dict[i][ii]['Use'][iii]][wr_pou_dict[i][ii]['Source']] + temp_sum
                            pods_time_temp = np.asarray(pods_time_temp)
                            
                            if len(np.unique(pods_time_temp[:,0])) > 1 or len(np.unique(pods_time_temp[:,1])) > 1:
                                start_date = np.min(pods_time_temp[:,0])
                                end_date = np.max(pods_time_temp[:,1])
                            else:
                                start_date = pods_time_temp[0,0]
                                end_date = pods_time_temp[0,1]
                            
                            jjj = pods_temp[len(pods_temp)-1]
                            if num_pods > 1:
                                if len(comment) == 0:
                                    comment = 'Multiple PODs'
                                else:
                                    comment = comment + ', Multiple PODs'
                                
                                if len(np.unique(pods_time_temp[:,0])) > 1 or len(np.unique(pods_time_temp[:,1])) > 1:
                                    comment = comment + ', Multiple Pumping Dates'
                                
                                sub_pods = dict()
                                for subpodid in pods_temp:
                                    if wr_pod_dict[i][subpodid]['subbasin'] not in sub_pods:  
                                        sub_pods[wr_pod_dict[i][subpodid]['subbasin']] = subpodid
                                        
                                if len(sub_pods.keys()) > 1:
                                    for sub_pod_id in sub_pods.keys():
                                        print i, jjj, sub_pod_id
                                        jjj = sub_pods[sub_pod_id]
                                        
                                        non_empty_acres.append((2,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii], wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['priority'],
                                                    '-',num_pods,wr_pod_dict[i][jjj]['Unmod_Use'], wr_pod_dict[i][jjj]['use_code'],wr_pod_dict[i][jjj]['source_type'],wr_pod_dict[i][jjj]['wr_type_2'],
                                                    '-',temp_sum,wr_pod_dict[i][jjj]['max_rate'],start_date,end_date,sub_pod_id,comment)) 
                                    
                                else:
                                    non_empty_acres.append((2,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii], wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['priority'],
                                                    '-',num_pods,wr_pod_dict[i][jjj]['Unmod_Use'], wr_pod_dict[i][jjj]['use_code'],wr_pod_dict[i][jjj]['source_type'],wr_pod_dict[i][jjj]['wr_type_2'],
                                                    '-',temp_sum,wr_pod_dict[i][jjj]['max_rate'],start_date,end_date,wr_pod_dict[i][jjj]['subbasin'],comment)) 
                            else:
                                
                                non_empty_acres.append((2,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii], wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['priority'],
                                                    '-',pods_temp[0],wr_pod_dict[i][jjj]['Unmod_Use'], wr_pod_dict[i][jjj]['use_code'],wr_pod_dict[i][jjj]['source_type'],wr_pod_dict[i][jjj]['wr_type_2'],
                                                    '-',temp_sum,wr_pod_dict[i][jjj]['max_rate'],start_date,end_date,wr_pod_dict[i][jjj]['subbasin'],comment)) 
                            
                                
                        else:
                            if find_use_case_match == 0:
                                empty_acres.append((2,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],'-', 'No POD with same USE CODE OR WR Source'))
                            else:
                                num_pods = 0
                                pods_temp = []
                                pods_time_temp = []
                                
                                for jjj in wr_pod_dict[i].keys():
                                    use_case_match = 0
                                    if wr_pou_dict[i][ii]['Use'][iii] == wr_pod_dict[i][jjj]['use_code'] or wr_pou_dict[i][ii]['Unmod_Use'] == wr_pod_dict[i][jjj]['Unmod_Use']:
                                        use_case_match = 1
                                    
                                    if wr_pou_dict[i][ii]['Use'][iii] == 'STORAGE' and wr_pod_dict[i][jjj]['use_code'] == 'IRRIGATION':
                                        use_case_match = 1
                                    
                                    if wr_pou_dict[i][ii]['Use'][iii] == 'DOMESTIC' and wr_pod_dict[i][jjj]['Unmod_Use'] == 'IRRIGATION AND DOMESTIC':
                                        use_case_match = 1
                                
                                    if (wr_pou_dict[i][ii]['Source'] == wr_pod_dict[i][jjj]['wr_type'] and use_case_match == 1):
                                        num_pods = num_pods + 1
                                        pods_temp.append(jjj)
                                        pods_time_temp.append([wr_pod_dict[i][jjj]['pump_start'],wr_pod_dict[i][jjj]['pump_end']])
                                
                                pods_time_temp = np.asarray(pods_time_temp)
                                if len(np.unique(pods_time_temp[:,0])) > 1 or len(np.unique(pods_time_temp[:,1])) > 1:
                                    start_date = np.min(pods_time_temp[:,0])
                                    end_date = np.max(pods_time_temp[:,1])
                                else:
                                    start_date = pods_time_temp[0,0]
                                    end_date = pods_time_temp[0,1]
                                
                                jjj = pods_temp[len(pods_temp)-1]
                                if len(comment) == 0:
                                    comment = 'Max rate (cfs)'
                                else:
                                    comment = comment + ', Max rate (cfs)'
                                
                                if len(np.unique(pods_time_temp[:,0])) > 1 or len(np.unique(pods_time_temp[:,1])) > 1:
                                    comment = comment + ', Multiple Pumping Dates'
                                    
                                non_empty_acres.append((2,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii], wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['priority'],
                                                    '-',pods_temp[0],wr_pod_dict[i][jjj]['Unmod_Use'], wr_pod_dict[i][jjj]['use_code'],wr_pod_dict[i][jjj]['source_type'],wr_pod_dict[i][jjj]['wr_type_2'],
                                                    '-','-',wr_pod_dict[i][jjj]['max_rate'],start_date,end_date,wr_pod_dict[i][jjj]['subbasin'],comment))
                                    
                    else:
                        empty_acres.append((2,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],'-','No POD Match'))
                        
                else:
                    
                    if i in wr_pod_dict.keys():
                        old_duty = -999
                        num_pods = 0
                        duty_total = 0
                        pods_temp = []
                        pods_time_temp = []
                        
                        for jjj in wr_pod_dict[i].keys():
                            use_case_match = 0
                            if wr_pou_dict[i][ii]['Use'][iii] == wr_pod_dict[i][jjj]['use_code'] or wr_pou_dict[i][ii]['Unmod_Use'] == wr_pod_dict[i][jjj]['Unmod_Use']:
                                use_case_match = 1
                            if wr_pou_dict[i][ii]['Use'][iii] == 'STORAGE' and wr_pod_dict[i][jjj]['use_code'] == 'IRRIGATION':
                                use_case_match = 1
                            if wr_pou_dict[i][ii]['Use'][iii] == 'SUPPLEMENTAL IRRIGATION' and wr_pod_dict[i][jjj]['use_code'] == 'IRRIGATION':
                                use_case_match = 1
                            if wr_pou_dict[i][ii]['Use'][iii] == 'IRRIGATION' and wr_pod_dict[i][jjj]['use_code'] == 'SUPPLEMENTAL IRRIGATION':
                                use_case_match = 1
                                
                            if wr_pou_dict[i][ii]['Use'][iii] == 'DOMESTIC' and wr_pod_dict[i][jjj]['Unmod_Use'] == 'IRRIGATION AND DOMESTIC':
                                use_case_match = 1
    
                            
                            if (wr_pou_dict[i][ii]['Source'] == wr_pod_dict[i][jjj]['wr_type'] and use_case_match == 1):
                                if 'duty' in wr_pod_dict[i][jjj].keys():
                                    
                                    temp = wr_pod_dict[i][jjj]['duty']
                                    if old_duty != temp and num_pods > 1:
                                        print i, ii, wr_pod_dict[i][jjj]['use_code'], temp
                                        comment = 'Mutiple values for duty'
                                        
                                    if num_pods == 0:
                                        old_duty = temp
                                    duty_total = duty_total + temp
                                    
                                num_pods = num_pods + 1
                                pods_temp.append(jjj)
                                pods_time_temp.append([wr_pod_dict[i][jjj]['pump_start'],wr_pod_dict[i][jjj]['pump_end']])
                        
                        if len(pods_temp) > 0:
                            jjj = pods_temp[len(pods_temp)-1]
                        
                        if old_duty != -999:
                            wris_per_use[wr_pou_dict[i][ii]['Use'][iii]][wr_pou_dict[i][ii]['Source']] = wris_per_use[wr_pou_dict[i][ii]['Use'][iii]][wr_pou_dict[i][ii]['Source']] + wr_pou_dict[i][ii]['acre_feet']*old_duty
                            pods_time_temp = np.asarray(pods_time_temp)
                            
                            if len(np.unique(pods_time_temp[:,0])) > 1:
                                start_date = np.min(pods_time_temp[:,0])
                                end_date = np.max(pods_time_temp[:,1])
                            else:
                                start_date = pods_time_temp[0,0]
                                end_date = pods_time_temp[0,1]
                            
                            if num_pods > 1:
                                if len(comment) == 0:
                                    comment = 'Multiple PODs'
                                else:
                                    comment = comment + ', Multiple PODs'
                                
                                if len(np.unique(pods_time_temp[:,0])) > 1 or len(np.unique(pods_time_temp[:,1])) > 1:
                                    comment = comment + ', Multiple Pumping Dates'
                                    
                                sub_pods = dict()
                                for subpodid in pods_temp:
                                    if wr_pod_dict[i][subpodid]['subbasin'] not in sub_pods:  
                                        sub_pods[wr_pod_dict[i][subpodid]['subbasin']] = subpodid
                                if len(sub_pods.keys()) > 1:
                                    for sub_pod_id in sub_pods.keys():
                                        print i, jjj, sub_pod_id
                                        jjj = sub_pods[sub_pod_id]
                                        non_empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['priority'],
                                                    wr_pou_dict[i][ii]['acre_feet'],num_pods,wr_pod_dict[i][jjj]['Unmod_Use'],wr_pod_dict[i][jjj]['use_code'],wr_pod_dict[i][jjj]['source_type'],wr_pod_dict[i][jjj]['wr_type_2'],
                                                    old_duty,wr_pou_dict[i][ii]['acre_feet']*old_duty,wr_pod_dict[i][jjj]['max_rate'],start_date,end_date,sub_pod_id,comment))
                                else:    
                                    non_empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['priority'],
                                                    wr_pou_dict[i][ii]['acre_feet'],num_pods,wr_pod_dict[i][jjj]['Unmod_Use'],wr_pod_dict[i][jjj]['use_code'],wr_pod_dict[i][jjj]['source_type'],wr_pod_dict[i][jjj]['wr_type_2'],
                                                    old_duty,wr_pou_dict[i][ii]['acre_feet']*old_duty,wr_pod_dict[i][jjj]['max_rate'],start_date,end_date,wr_pod_dict[i][jjj]['subbasin'],comment))
                            else:
                                non_empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['priority'],
                                                    wr_pou_dict[i][ii]['acre_feet'],pods_temp[0],wr_pod_dict[i][jjj]['Unmod_Use'],wr_pod_dict[i][jjj]['use_code'],wr_pod_dict[i][jjj]['source_type'],wr_pod_dict[i][jjj]['wr_type_2'],
                                                    old_duty,wr_pou_dict[i][ii]['acre_feet']*old_duty,wr_pod_dict[i][jjj]['max_rate'],start_date,end_date,wr_pod_dict[i][jjj]['subbasin'],comment))
                        else:
                            if len(pods_temp) > 0:
    
                                if 'IRRIGATION' in wr_pou_dict[i][ii]['Use'][iii] or 'SUPPLEMENTAL IRRIGATION' in wr_pou_dict[i][ii]['Use'][iii]:
                                    if len(comment) == 0:
                                        comment = 'Assumed Duty of 3.5 acre-ft'
                                    else:
                                        comment = comment + ', Assumed Duty of 3.5 acre-ft'
                                    
                                    non_empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['priority'],
                                                    wr_pou_dict[i][ii]['acre_feet'],pods_temp[0],wr_pod_dict[i][jjj]['Unmod_Use'],wr_pod_dict[i][jjj]['use_code'],wr_pod_dict[i][jjj]['source_type'],wr_pod_dict[i][jjj]['wr_type_2'],
                                                    3.5,wr_pou_dict[i][ii]['acre_feet']*3.5,wr_pod_dict[i][jjj]['max_rate'],start_date,end_date,wr_pod_dict[i][jjj]['subbasin'],comment))
                                else:

                                    if len(comment) == 0:
                                        comment = 'No Duty'
                                    else:
                                        comment = comment + ', No Duty'
                                    non_empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['priority'],
                                                        wr_pou_dict[i][ii]['acre_feet'],pods_temp[0],wr_pod_dict[i][jjj]['Unmod_Use'],wr_pod_dict[i][jjj]['use_code'],wr_pod_dict[i][jjj]['source_type'],wr_pod_dict[i][jjj]['wr_type_2'],
                                                        '-','-',wr_pod_dict[i][jjj]['max_rate'],start_date,end_date,wr_pod_dict[i][jjj]['subbasin'],comment))
                            else:
                                empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['acre_feet'],'No POD with same USE CODE OR WR Source'))
                    else:
                        empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['acre_feet'],'No POD Match'))
    
    return non_empty_acres, empty_acres

#%%
    
def GetInstream_WR(wr_pod_dict, sub_pod_dict):
    
    instream_wrs = dict()
    for snpid in wr_pod_dict.keys():
        for podid in wr_pod_dict[snpid].keys():
             if 'instream' in wr_pod_dict[snpid][podid]['use_code'].lower(): 
            #if sub_pod_dict[str(snpid)][str(podid)] != 'nan' and 'instream' in wr_pod_dict[snpid][podid]['use_code'].lower(): 
                if snpid not in instream_wrs.keys():
                    instream_wrs[snpid] = dict()
                
                instream_wrs[snpid][podid] = wr_pod_dict[snpid][podid]
                
    return instream_wrs

