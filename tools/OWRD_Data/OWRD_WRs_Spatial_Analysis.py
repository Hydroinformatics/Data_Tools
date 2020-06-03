# -*- coding: utf-8 -*-

import sys, datetime, gdal, os, shutil, zipfile, json, shapely
from statistics import mode
from scipy import stats, ndimage
qswat_utils_path = r'C:\Users\sammy\Documents\GitHub\InterACTWEL\src'
sys.path.append(qswat_utils_path)

from qswat import QSWAT_utils

from osgeo import ogr
from shapely.geometry import Point

import matplotlib.pyplot as plt
import numpy as np

import rasterio
from rasterio.mask import mask

#%%

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

#%%

def EliminateDuplicates(root_path, input_file_name, output_file_name):    
    
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(root_path + '/' + input_file_name, 0)
    wr_layer = dataSource.GetLayer()
    
    spatialRef = wr_layer.GetSpatialRef()
    
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
                    
                    if fname == 'use_code_d' and field_val == 'SUPPLEMENTAL IRRIGATION':
                        field_val = 'IRRIGATION'
                    if fname == 'use_code' and field_val == 'IS':
                        field_val = 'IR'
                        
                    if field_val == feature_dict[snp_id][objs][fname]: 
                        feat_sim_count = feat_sim_count + 1
                    else:
                        non_equal_feat.append(fname)
                
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
 
#################################################################### 
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
def GetWR_POD_Data(file_name, root_path, pou_keys):
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
            if feature.GetField("tributary_"):
                wr_pod_dict[int(snp_id)][int(pod_id)]['tributary'] = feature.GetField("tributary_").strip()
            else:
                wr_pod_dict[int(snp_id)][int(pod_id)]['tributary'] = 'Not specified'
             
            
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
            #wr_pod_dict[int(snp_id)][int(pod_id)]['pump_end'] = 
            
    return wr_pod_dict


#%%
def WR_Parser(wr_pou_dict, wr_pod_dict):
    num_pou_recs=0
    num_pou_nan = []
    for i in wr_pou_dict.keys():    
        for ii in wr_pou_dict[i].keys():
            if np.isnan(wr_pou_dict[i][ii]['acre_feet']) == 1:
                #num_pou_nan = num_pou_nan + 1
                num_pou_nan.append((i,ii,0))
            else:
                num_pou_nan.append((i,ii,1))
            
            num_pou_recs = num_pou_recs + 1
            
#    bb = np.asarray(num_pou_nan)

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
                    #empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Use']))

                    if i in wr_pod_dict.keys():
                        temp_sum = 0
                        num_pods = 0
                        pods_temp = []
                        pods_time_temp = []
                        find_use_case_match = 0
                        #temp_sum = dict()
                        
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
                                #wr_pou_dict[i][ii]['Source'] == pod_source_dict[wr_pod_dict[i][jjj]['source_type']] 
                                
                                #if pod_source_dict[wr_pod_dict[i][iii]['source_type']] not in temp_sum.keys():
                                #    temp_sum[pod_source_dict[wr_pod_dict[i][iii]['source_type']]] = 0
                                
                                temp_sum = temp_sum + wr_pod_dict[i][jjj]['acre_feet']
                                num_pods = num_pods + 1
                                pods_temp.append(jjj)
                                pods_time_temp.append([wr_pod_dict[i][jjj]['pump_start'],wr_pod_dict[i][jjj]['pump_end']])
                                #temp_sum[pod_source_dict[wr_pod_dict[i][iii]['source_type']]] = temp_sum[pod_source_dict[wr_pod_dict[i][iii]['source_type']]] + wr_pod_dict[i][iii]['acre_feet']
                                
                        if temp_sum != 0:
                        #if len(temp_sum.keys()) != 0:
                        #    for zi in temp_sum.keys():
                        #        if zi not in wris_per_use[wr_pou_dict[i][ii]['Use']].keys():
                        #            wris_per_use[wr_pou_dict[i][ii]['Use']][zi] = 0
                                    
                            wris_per_use[wr_pou_dict[i][ii]['Use'][iii]][wr_pou_dict[i][ii]['Source']] = wris_per_use[wr_pou_dict[i][ii]['Use'][iii]][wr_pou_dict[i][ii]['Source']] + temp_sum
                            #non_empty_acres.append((i,ii,wr_pou_dict[i][ii]['Use']))
                            
                            pods_time_temp = np.asarray(pods_time_temp)
                            if len(np.unique(pods_time_temp[:,0])) > 1:
                                start_date = np.min(pods_time_temp[:,0])
                                end_date = np.max(pods_time_temp[:,1])
#                            elif len(pods_time_temp[:,0]) > 1:
#                                start_date = mode(pods_time_temp[:,0])
#                                end_date = mode(pods_time_temp[:,1])
                            else:
                                start_date = pods_time_temp[0,0]
                                end_date = pods_time_temp[0,1]
                            
                            if num_pods > 1:
                                jjj = pods_temp[len(pods_temp)-1]
                                if len(comment) == 0:
                                    comment = 'Multiple PODs'
                                else:
                                    comment = comment + ', Multiple PODs'
                                
                                if len(np.unique(pods_time_temp[:,0])) > 1 or len(np.unique(pods_time_temp[:,1])) > 1:
                                    comment = comment + ', Multiple Pumping Dates'
                                
                                non_empty_acres.append((2,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii], wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['priority'],
                                                    '-',num_pods,wr_pod_dict[i][jjj]['Unmod_Use'], wr_pod_dict[i][jjj]['use_code'],wr_pod_dict[i][jjj]['source_type'],wr_pod_dict[i][jjj]['wr_type'],
                                                    '-',temp_sum,wr_pod_dict[i][jjj]['max_rate'],start_date,end_date,comment)) 
                            else:
                                non_empty_acres.append((2,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii], wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['priority'],
                                                    '-',pods_temp[0],wr_pod_dict[i][jjj]['Unmod_Use'], wr_pod_dict[i][jjj]['use_code'],wr_pod_dict[i][jjj]['source_type'],wr_pod_dict[i][jjj]['wr_type'],
                                                    '-',temp_sum,wr_pod_dict[i][jjj]['max_rate'],start_date,end_date,comment)) 
                            
                                
                            
                        else:
                            if find_use_case_match == 0:
                                empty_acres.append((2,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],'-', 'No POD with same USE CODE OR WR Source'))
                            else:
                                #empty_acres.append((2,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],'-','No acre-feet'))
    #                            if num_pods > 1:
    #                                jjj = pods_temp[len(pods_temp)-1]
    #                                if len(comment) == 0:
    #                                    comment = 'Multiple PODs'
    #                                else:
    #                                    comment = comment + ', Multiple PODs'
                                        
    #                                non_empty_acres.append((2,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii], wr_pou_dict[i][ii]['Source'],
    #                                                'None',num_pods,wr_pod_dict[i][jjj]['Unmod_Use'], wr_pod_dict[i][jjj]['use_code'],wr_pod_dict[i][jjj]['source_type'],wr_pod_dict[i][jjj]['wr_type'],
    #                                                'None',temp_sum,comment)) 
    #                            else:
                                
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
                                if len(np.unique(pods_time_temp[:,0])) > 1:
                                    start_date = np.min(pods_time_temp[:,0])
                                    end_date = np.max(pods_time_temp[:,1])
#                                elif len(pods_time_temp[:,0]) > 1:
#                                    start_date = mode(pods_time_temp[:,0])
#                                    end_date = mode(pods_time_temp[:,1])
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
                                                    '-',pods_temp[0],wr_pod_dict[i][jjj]['Unmod_Use'], wr_pod_dict[i][jjj]['use_code'],wr_pod_dict[i][jjj]['source_type'],wr_pod_dict[i][jjj]['wr_type'],
                                                    '-','-',wr_pod_dict[i][jjj]['max_rate'],start_date,end_date,comment))
                                    
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
                                    
                        if old_duty != -999:
                            #wris_per_use[wr_pou_dict[i][ii]['Use'][iii]][wr_pou_dict[i][ii]['Source']] = wris_per_use[wr_pou_dict[i][ii]['Use'][iii]][wr_pou_dict[i][ii]['Source']] + wr_pou_dict[i][ii]['acre_feet']*temp
                            wris_per_use[wr_pou_dict[i][ii]['Use'][iii]][wr_pou_dict[i][ii]['Source']] = wris_per_use[wr_pou_dict[i][ii]['Use'][iii]][wr_pou_dict[i][ii]['Source']] + wr_pou_dict[i][ii]['acre_feet']*old_duty
                            #non_empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Use'],wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['acre_feet']))
                            pods_time_temp = np.asarray(pods_time_temp)
                            
                            if len(np.unique(pods_time_temp[:,0])) > 1:
                                start_date = np.min(pods_time_temp[:,0])
                                end_date = np.max(pods_time_temp[:,1])
#                            elif len(pods_time_temp[:,0]) > 1:
#                                start_date = mode(pods_time_temp[:,0])
#                                end_date = mode(pods_time_temp[:,1])
                            else:
                                start_date = pods_time_temp[0,0]
                                end_date = pods_time_temp[0,1]
                            
                            if num_pods > 1:
                                jjj = pods_temp[len(pods_temp)-1]
                                if len(comment) == 0:
                                    comment = 'Multiple PODs'
                                else:
                                    comment = comment + ', Multiple PODs'
                                
                                if len(np.unique(pods_time_temp[:,0])) > 1 or len(np.unique(pods_time_temp[:,1])) > 1:
                                    comment = comment + ', Multiple Pumping Dates'
                                    
                                non_empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['priority'],
                                                    wr_pou_dict[i][ii]['acre_feet'],num_pods,wr_pod_dict[i][jjj]['Unmod_Use'],wr_pod_dict[i][jjj]['use_code'],wr_pod_dict[i][jjj]['source_type'],wr_pod_dict[i][jjj]['wr_type'],
                                                    old_duty,wr_pou_dict[i][ii]['acre_feet']*old_duty,wr_pod_dict[i][jjj]['max_rate'],start_date,end_date,comment))
                            else:
                                non_empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['priority'],
                                                    wr_pou_dict[i][ii]['acre_feet'],pods_temp[0],wr_pod_dict[i][pods_temp[0]]['Unmod_Use'],wr_pod_dict[i][pods_temp[0]]['use_code'],wr_pod_dict[i][pods_temp[0]]['source_type'],wr_pod_dict[i][pods_temp[0]]['wr_type'],
                                                    old_duty,wr_pou_dict[i][ii]['acre_feet']*old_duty,wr_pod_dict[i][jjj]['max_rate'],start_date,end_date,comment))
                        else:
                            if len(pods_temp) > 0:
    
                                if 'IRRIGATION' in wr_pou_dict[i][ii]['Use'][iii] or 'SUPPLEMENTAL IRRIGATION' in wr_pou_dict[i][ii]['Use'][iii]:
                                    if len(comment) == 0:
                                        comment = 'Assumed Duty of 3 acre-ft'
                                    else:
                                        comment = comment + ', Assumed Duty of 3 acre-ft'
                                    
                                    non_empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['priority'],
                                                    wr_pou_dict[i][ii]['acre_feet'],pods_temp[0],wr_pod_dict[i][pods_temp[0]]['Unmod_Use'],wr_pod_dict[i][pods_temp[0]]['use_code'],wr_pod_dict[i][pods_temp[0]]['source_type'],wr_pod_dict[i][pods_temp[0]]['wr_type'],
                                                    3,wr_pou_dict[i][ii]['acre_feet']*3,wr_pod_dict[i][jjj]['max_rate'],start_date,end_date,comment))
                                else:
                                    if len(comment) == 0:
                                        comment = 'No Duty'
                                    else:
                                        comment = comment + ', No Duty'
                                    non_empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['priority'],
                                                        wr_pou_dict[i][ii]['acre_feet'],pods_temp[0],wr_pod_dict[i][pods_temp[0]]['Unmod_Use'],wr_pod_dict[i][pods_temp[0]]['use_code'],wr_pod_dict[i][pods_temp[0]]['source_type'],wr_pod_dict[i][pods_temp[0]]['wr_type'],
                                                        '-','-',wr_pod_dict[i][jjj]['max_rate'],start_date,end_date,comment))
                            else:
                                empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['acre_feet'],'No POD with same USE CODE OR WR Source'))
                    else:
                        empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['acre_feet'],'No POD Match'))
    
    return non_empty_acres, empty_acres

#%%
#root_path = r'Z:\Projects\INFEWS\Modeling\GIS_Data\Data\wr_state'
root_path = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\HRUs_Meghna'
pou_file = 'wr_v_pou_public_proj_REGION.shp'
pod_file = 'wr_v_pod_public_proj_REGIONv2.shp'
pou_file_noduplicate = 'wr_v_pou_public_proj_REGION_NoDuplicatesPy_V2'

#root_path = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\DMW'
#pou_file = 'pou_dmw.shp'
#pod_file = 'pod_dmw.shp'
#pou_file_noduplicate = 'pou_dmw_NoDuplicatesPy_V2'

feature_dict, remove_features = EliminateDuplicates(root_path, pou_file, pou_file_noduplicate)

pou_file = pou_file_noduplicate + '.shp'

wr_pou_dict = GetWR_POU_Data(pou_file, root_path)
wr_pod_dict = GetWR_POD_Data(pod_file, root_path, wr_pou_dict.keys())

non_empty_acres, empty_acres = WR_Parser(wr_pou_dict, wr_pod_dict)

non_empty_acres_arr = np.asarray(non_empty_acres)
empty_acres_arr = np.asarray(empty_acres)

hru_wr_area_tresh = 5

#%%
cdl_path = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\Umatilla_Input_data\CDL_Data\CDL_Data_Filter'

cdl_file = 'CDL_2017_clip_20180821151150_mfv2.tif'
        
#%%

driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(root_path + '/' + pou_file, 0)
wr_layer = dataSource.GetLayer()

spatialRef = wr_layer.GetSpatialRef()

#pou_uses = dict()
#for feature in wr_layer:
#    use_cd = feature.GetField("use_code_d")
#    wr_type = feature.GetField("wr_type")
#    
#    if use_cd not in pou_uses.keys():
#        pou_uses[use_cd] = dict()
#    
#    if wr_type not in pou_uses[use_cd].keys():
#        pou_uses[use_cd][wr_type] = []
#    
#    snp_id = feature.GetField("snp_id")
#    pou_uses[use_cd][wr_type].append(snp_id)
#    
#wr_layer.ResetReading()

pou_uses = dict()
for i in range(len(non_empty_acres)):
    use_cd = non_empty_acres[i][4]
    wr_type = non_empty_acres[i][5]
    
    if use_cd not in pou_uses.keys():
        pou_uses[use_cd] = dict()
    
    if wr_type not in pou_uses[use_cd].keys():
        pou_uses[use_cd][wr_type] = []
    
    snp_id = non_empty_acres[i][1]
    pou_uses[use_cd][wr_type].append(snp_id)
    
#%%
outPath = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\HRUs_Meghna\Zonal_Stats\WR_Shapefiles'

clayer = 0
#temp_use = ['IRRIGATION','SUPPLEMENTAL IRRIGATION']
temp_use = ['IRRIGATION']
temp_wr_type = ['SW']

total_wr_vol = dict()

#for use_code in pou_uses.keys():
for use_code in temp_use:
    if use_code not in total_wr_vol.keys():
        total_wr_vol[use_code] = dict()
        
    for wr_type in pou_uses[use_code].keys():
    #for wr_type in temp_wr_type:
        if wr_type not in total_wr_vol[use_code].keys():
            total_wr_vol[use_code][wr_type] = dict()
        
        outFileName = 'wr_v_pou_public_proj_REGION_' + str(clayer) + '_' + str(wr_type) + 'v2.shp'
        print outFileName
        
        outShapefile = outPath + '/' + outFileName
        outDriver = ogr.GetDriverByName("ESRI Shapefile")
    
        # Remove output shapefile if it already exists
        if os.path.exists(outShapefile):
            outDriver.DeleteDataSource(outShapefile)
        
        proj_outFileName = '\wr_v_pou_public_proj_REGION_' + str(clayer) + '_' + str(wr_type) + 'v2.prj'
        with open(outPath + proj_outFileName, 'w') as prj_file:
            prj_file.write(str(spatialRef.ExportToWkt()))
        #prj_file.write(spatialRef)
        prj_file.close()
    
        # Create the output shapefile
        outDataSource = outDriver.CreateDataSource(outShapefile)
        outLayer = outDataSource.CreateLayer("wr_v_pou_public", geom_type = ogr.wkbPolygon)
    
        # Add an ID field
        idField = ogr.FieldDefn("id", ogr.OFTInteger)
        outLayer.CreateField(idField)
        
        idField = ogr.FieldDefn("snp_id", ogr.OFTInteger)
        outLayer.CreateField(idField)
        
        idField = ogr.FieldDefn("pou_use_id", ogr.OFTInteger)
        outLayer.CreateField(idField)
        
        idField = ogr.FieldDefn("wr_type", ogr.OFTString)
        outLayer.CreateField(idField)
        
        idField = ogr.FieldDefn("use_code_d", ogr.OFTString)
        outLayer.CreateField(idField)
        
        idField = ogr.FieldDefn("duty", ogr.OFTReal)
        outLayer.CreateField(idField)
        
        idField = ogr.FieldDefn("wr_vol", ogr.OFTReal)
        outLayer.CreateField(idField)
        
        idField = ogr.FieldDefn("wris_acres", ogr.OFTReal)
        outLayer.CreateField(idField)
        
        idField = ogr.FieldDefn("poly_area", ogr.OFTReal)
        outLayer.CreateField(idField)
        
        
        ft_idc = 0
        for feature in wr_layer:
            snp_id = feature.GetField("snp_id")
            pou_id = feature.GetField("pou_use_id")
            use_cd = feature.GetField("use_code_d")
            wr_typ = feature.GetField("wr_type")
            wr_acres = feature.GetField("wris_acres")
            
            if snp_id in pou_uses[use_code][wr_type] and use_cd == use_code and wr_typ == wr_type:    
                points = feature.GetGeometryRef()
                
                # Create the feature and set values
                featureDefn = outLayer.GetLayerDefn()
                featureo = ogr.Feature(featureDefn)
        
                featureo.SetField("id", ft_idc)
                featureo.SetField("snp_id", snp_id)
                featureo.SetField("pou_use_id", pou_id)
                featureo.SetField("wr_type", wr_typ)
                featureo.SetField("use_code_d", use_cd)
                featureo.SetField("wris_acres", wr_acres)
                
                wr_vol = -999.0
                wr_duty = -999.0
                for i in range(len(non_empty_acres)):
                    if non_empty_acres[i][1] == snp_id and non_empty_acres[i][4] == use_cd and non_empty_acres[i][5] == wr_type:
                        wr_duty = non_empty_acres[i][13]
                        wr_vol = non_empty_acres[i][14]
                
                #print snp_id, use_cd, wr_type, wr_vol
                featureo.SetField("duty", wr_duty)
                featureo.SetField("wr_vol", wr_vol)
                
                if snp_id not in total_wr_vol[use_code][wr_type].keys():
                    total_wr_vol[use_code][wr_type][snp_id] = 0
                
                total_wr_vol[use_code][wr_type][snp_id] = total_wr_vol[use_code][wr_type][snp_id] + wr_vol
                
                temp_area = []
                diam_bool = 0
                len_pts = 0
                for pt in points:
                    temp_area.append(pt.GetArea())
                    len_pts = len_pts + 1
                
                if np.mean(temp_area) > 220.0 and np.mean(temp_area) < 225.0:
                    diam_bool = 1
                    if len_pts > 1:
                        wkt = 'MULTIPOLYGON ('
                    else:
                        wkt = 'POLYGON '
                
                #print points
                pts = []
                area = []
                cpoly = 1
                if diam_bool == 1:
                    for pt in points:
                        if not pt.GetPoints():
                            ptsc = pt.GetGeometryRef(0)
                            ptsc = np.asarray(ptsc.GetPoints())
                        else:
                            ptsc = np.asarray(pt.GetPoints())
                        
                        xc, yc = centeroidnp(ptsc)
                        #wkt = 'POINT (' + str(xc) + ' ' + str(yc) + ')'
                        centroid_pt = Point(xc,yc)
                        
                        xlen = 407.0/2
                        ylen = 405.0/2
#                        
#                        #wkt = 'LINEARRING ('
                        wkt = 'POLYGON '
                        wkt = wkt + '(('
                        wkt = wkt + str(xc+xlen) + ' ' + str(yc+ylen) + ','
                        wkt = wkt + str(xc+xlen) + ' ' + str(yc-ylen) + ','
                        wkt = wkt + str(xc-xlen) + ' ' + str(yc-ylen) + ','
                        wkt = wkt + str(xc-xlen) + ' ' + str(yc+ylen) + ','
                        wkt = wkt + str(xc+xlen) + ' ' + str(yc+ylen) + '))'
                        poly = ogr.CreateGeometryFromWkt(wkt)
#                    
#                    if len_pts > 1:
#                        wkt = wkt + ')'
                        ##bufferDistance = 405.0
                        ##ptsc = ogr.CreateGeometryFromWkt(wkt)
                        ##poly = ptsc.Buffer(bufferDistance,cap_style=3)
                        
                        #bufferDistance = 405.0/2
                        #buf = centroid_pt.buffer(bufferDistance, cap_style=3)
                        #poly = ogr.CreateGeometryFromWkt(buf.wkt)
                        
                        if cpoly == 1:
                            poly2 = poly
                            #poly2 = ogr.Geometry(ogr.wkbMultiPolygon)
                            cpoly = 0
                        else:
                            poly2 = poly2.Union(poly)
                        
                    #print wkt
                    #ptsc = ogr.CreateGeometryFromWkt(wkt)
                    #poly = ptsc.Buffer(bufferDistance)
                    #featureo.SetGeometry(ptsc)
                    area_poly = poly2.GetArea()/4046.86
                    featureo.SetField("poly_area", area_poly)
                    featureo.SetGeometry(poly2)
                else:
                    area_poly = points.GetArea()/4046.86
                    featureo.SetField("poly_area", area_poly)
                    featureo.SetGeometry(points)
                    
                outLayer.CreateFeature(featureo)
                
                ft_idc = ft_idc + 1
                
        wr_layer.ResetReading()
        
        feature = None
        
        # Save and close DataSource
        outDataSource = None
    clayer = clayer + 1
   
#%%
#fnames = os.listdir(cdl_path)  
#no_crop_ids = Get_No_Crop_Ids()
#    
#raster_file = cdl_path + '\\' + cdl_file
#region, region_NoData, region_obj = QSWAT_utils.Read_Raster(raster_file)
#region_prj = region_obj.GetProjection()
#
#output_path = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\HRUs_Meghna\Zonal_Stats\WR_Tiffs'
#
#hru_wr = dict()
#
#clayer = 0
##for use_code in pou_uses.keys():
#for use_code in temp_use:
#    if use_code not in hru_wr.keys():
#        hru_wr[use_code] = dict()
#        
#    for wr_type in pou_uses[use_code].keys():
#    #for wr_type in temp_wr_type:
#        if wr_type not in hru_wr[use_code].keys():
#            hru_wr[use_code][wr_type] = dict()
#            
#        new_wrs = np.zeros(region.shape, dtype=np.int64)
#        new_wrs_duty = np.zeros(region.shape, dtype=np.int64)
#        
#        outFileName = '\wr_v_pou_public_proj_REGION_' + str(clayer) + '_' + str(wr_type) + 'v2.shp'        
#        boundary_lyr = outPath + outFileName
#        
#        dataSource = driver.Open(boundary_lyr, 0)
#        temp_layer = dataSource.GetLayer()
#        
#        cfeat = 1
#        for feature in temp_layer:
#            outFileName = '/temp.shp'
#            outShapefile = output_path + outFileName
#            outDriver = ogr.GetDriverByName("ESRI Shapefile")
#        
#            # Remove output shapefile if it already exists
#            if os.path.exists(outShapefile):
#                outDriver.DeleteDataSource(outShapefile)
#                
#            proj_outFileName = '/temp.prj'
#            with open(output_path + proj_outFileName, 'w') as prj_file:
#                #prj_file.write(str(spatialRef.ExportToWkt()))
#                prj_file.write(str(region_prj))
#            prj_file.close()
#        
#        
#            # Create the output shapefile
#            outDataSource = outDriver.CreateDataSource(outShapefile)
#            outLayer = outDataSource.CreateLayer("temp", geom_type = ogr.wkbPolygon)
#        
#            # Add an ID field
#            idField = ogr.FieldDefn("id", ogr.OFTInteger)
#            outLayer.CreateField(idField)
#            
#            featureDefn = outLayer.GetLayerDefn()
#            featureo = ogr.Feature(featureDefn)
#                
#            featureo.SetField("id", 1)
#            points = feature.GetGeometryRef()
#            featureo.SetGeometry(points)
#            
#            feat_id = feature.GetField("id")
#            snp_id = feature.GetField("snp_id")
#            pou_id = feature.GetField("pou_use_id")
#            wr_vol = feature.GetField("wr_vol")
#            wr_duty = feature.GetField("duty")
#            
##            total_wr_vol[use_code][wr_type] = total_wr_vol[use_code][wr_type] + feature.GetField("wr_vol")
#                            
#            outLayer.CreateFeature(featureo)
#            outDataSource = None
#            
#            BoundaryRaster = output_path + '/temp_boundary.tif'
#            QSWAT_utils.CreateTiff(outShapefile, BoundaryRaster, raster_file,'',1,'','tiff')
#            
#            region_temp, region_NoData_temp, region_obj_temp = QSWAT_utils.Read_Raster(BoundaryRaster)
#            
#            if len(np.where(region_temp > 0)) == 0:
#                break
#            
#            # CHECK WHY THIS IS NEEDED
#            region_temp = np.asarray(np.flip(region_temp,0),dtype=float)
#            
#            if region_NoData_temp != None:
#                region_temp[region_temp == region_NoData_temp] = 0.0
#            
#            result, num_features = ndimage.measurements.label(region_temp)
#            labeled_area = ndimage.measurements.sum(region_temp,result,range(1,num_features+1))
#            
#            for bi in range(len(labeled_area)):
#                if labeled_area[bi] <= hru_wr_area_tresh and bi != np.argmax(labeled_area):
#                    region_temp[result == bi] = 0
#        
#            #print cfeat
#            if cfeat == 1:
#                new_wrs[region_temp > 0.0] = cfeat
#                hru_wr[use_code][wr_type][cfeat] = [(snp_id, pou_id)]
#                cfeat = cfeat + 1
#                
#            else:
#                print cfeat
##                if cfeat == 1459:
##                    stpo = 0
#                old_new_wrs = new_wrs.copy()
#                new_wrs[region_temp > 0.0] = cfeat
#                hru_wr[use_code][wr_type][cfeat] = [(snp_id, pou_id)]
#                
#                uhru = np.unique(old_new_wrs[region_temp > 0.0])
#                uhru = uhru[uhru != 0]
#                uhru = uhru[uhru != cfeat]
#                
#                cfeat = cfeat + 1
#                
#                if len(uhru) > 0:
#                    for u in uhru:
#                        temp = []
#                        temp.append((snp_id, pou_id))
#                        #if len(hru_wr[use_code][wr_type][u]):
#                        for usnp_id, upou_id in hru_wr[use_code][wr_type][u]:
#                            temp.append((usnp_id,upou_id))
#                            
#                        temp_cond = (region_temp > 0.0) & (old_new_wrs == u)
#                        temp_cond = np.asarray(temp_cond, dtype=int)
#                        result, num_features = ndimage.measurements.label(temp_cond)
#                        labeled_area = ndimage.measurements.sum(temp_cond,result,range(1,num_features+1))
#            
#                        for bi in range(len(labeled_area)):
#                            if labeled_area[bi] <= hru_wr_area_tresh:
#                                temp_cond[result == bi] = 0
#                    
#                        if np.sum(np.asarray(temp_cond.flatten(), dtype=int)) > hru_wr_area_tresh:
#                            hru_wr[use_code][wr_type][cfeat] = temp
#                            new_wrs[temp_cond == True] = cfeat
#                            cfeat = cfeat + 1
#                    
##                else:
##                    new_wrs[region_temp > 0.0] = cfeat
##                    hru_wr[use_code][wr_type][cfeat] = [snp_id]
##                    cfeat = cfeat + 1
#            
#                    
#        outRasterName = 'WR_POU_proj_REGION_' + str(clayer) + '_' + str(wr_type) + 'v2.tif'
#        print outRasterName
#        newRaster_name = output_path + '/' + outRasterName
#        QSWAT_utils.Save_NewRaster(new_wrs, region_obj, region, newRaster_name, -999.0)
#        
#    clayer = clayer + 1
#
#with open('WR_IRR_DICT_HRU_Masterv2.json', 'w') as fp:
#    json.dump(hru_wr, fp)
    
#%%
with open('WR_IRR_DICT_HRU_Masterv2.json') as json_file:
    hru_wr = json.load(json_file)
    
output_path = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\HRUs_Meghna\Zonal_Stats\WR_Tiffs'

output_path = output_path.replace('\\','/')

lu_raster_file = output_path + '/CDL_2017_clip_20180821151150_mfv2_clp.tif'
landuseFile_clip = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\Umatilla_Output_data_TEST\HRUs_LU\UHRUs_Region_Latest.tif'
landuseFile = output_path + '/UHRUs_Region_Latest_Clip.tif'
QSWAT_utils.Clipraster(landuseFile_clip, landuseFile, lu_raster_file, gdal.GRA_Mode)

lu_region, lu_region_NoData, lu_region_obj = QSWAT_utils.Read_Raster(landuseFile)
lu_region_prj = lu_region_obj.GetProjection()

if lu_region_NoData != None:
    lu_region[lu_region == lu_region_NoData] = 0.0

lu_raster_flat = lu_region.flatten()

#BoundaryRaster = output_path + '/hrus.tif'
#outShapefile = root_path + '\Correct_HRUs_11534_ExemptLU_333.shp'
#QSWAT_utils.CreateTiff(outShapefile, BoundaryRaster, lu_raster_file,'HRUGIS','','','tiff',gdal.GDT_Float64)
#hru_raster, hru_NoData, hru_obj_temp = QSWAT_utils.Read_Raster(BoundaryRaster)
#hru_raster = np.asarray(np.flip(hru_raster,0),dtype=float)

GW_WR_Raster = output_path + '/WR_POU_proj_REGION_0_GWv2.tif'
GW_WR_Raster_clip = output_path + '/WR_POU_GW_clip.tif'
QSWAT_utils.Clipraster(GW_WR_Raster, GW_WR_Raster_clip, landuseFile, gdal.GRA_Mode)
gwwr_raster, gwwr_NoData, gwwr_obj_temp = QSWAT_utils.Read_Raster(GW_WR_Raster_clip)

#hru_raster = np.asarray(np.flip(hru_raster,0),dtype=float)
#plt.matshow(hru_raster)
if gwwr_NoData != None:
    gwwr_raster[gwwr_raster == gwwr_NoData] = 0.0   
gwwr_flat = gwwr_raster.flatten()

SW_WR_Raster = output_path + '/WR_POU_proj_REGION_0_SWv2.tif'
SW_WR_Raster_clip = output_path + '/WR_POU_SW_Clip.tif'
QSWAT_utils.Clipraster(SW_WR_Raster, SW_WR_Raster_clip, landuseFile, gdal.GRA_Mode)
swwr_raster, swwr_NoData, swwr_obj_temp = QSWAT_utils.Read_Raster(SW_WR_Raster_clip)

#hru_raster = np.asarray(np.flip(hru_raster,0),dtype=float)
#plt.matshow(hru_raster)
if swwr_NoData != None:
    swwr_raster[swwr_raster == swwr_NoData] = 0.0
swwr_flat = swwr_raster.flatten()

#%%

swwr_ids = np.asarray(np.unique(swwr_flat, axis=0),dtype = int)
gwwr_ids = np.asarray(np.unique(gwwr_flat, axis=0),dtype = int)

unique_sw_snp_ids = []
for objid in swwr_ids:
    if str(objid) in hru_wr['IRRIGATION']['SW'].keys():
        for snp_ids, pouid in  hru_wr['IRRIGATION']['SW'][str(objid)]:
            if [int(snp_ids),int(pouid),1] not in unique_sw_snp_ids:
                unique_sw_snp_ids.append([int(snp_ids),int(pouid),1])
                #unique_sw_snp_ids.append([int(objid),int(snp_ids),int(pouid),1])

unique_gw_snp_ids = []              
for objid in gwwr_ids:
    if str(objid) in hru_wr['IRRIGATION']['GW'].keys():
        for snp_ids, pouid in  hru_wr['IRRIGATION']['GW'][str(objid)]:
            if [int(snp_ids),int(pouid),2] not in unique_gw_snp_ids:
                unique_gw_snp_ids.append([int(snp_ids),int(pouid),2])
                #unique_gw_snp_ids.append([int(objid),int(snp_ids),int(pouid),2])

unique_sw_snp_ids = np.asarray(unique_sw_snp_ids)
#unique_sw_snp_ids = unique_sw_snp_ids[:,[1,2,3]]

unique_gw_snp_ids = np.asarray(unique_gw_snp_ids)
#unique_gw_snp_ids = unique_gw_snp_ids[:,[1,2,3]]

unique_snp_ids = np.concatenate((unique_sw_snp_ids, unique_gw_snp_ids))

if len(unique_snp_ids) != (len(unique_sw_snp_ids[:,0]) + len(unique_gw_snp_ids[:,0])):
    print 'Warning: You might have GW and SW water rights sharing a snap_ID'

unique_snp_idsr, usnp_ids_counts = np.unique(unique_snp_ids, axis=0, return_counts = 1)

#%%

use_cd = 'IRRIGATION'
swat_to_wrtype = dict()
swat_to_wrtype['SW'] = 1
swat_to_wrtype['GW'] = 3
swat_to_wrtype['ST'] = 5

wr_swat_file = []
wr_snp_ids = np.asarray(non_empty_acres_arr[:,[1,2]], dtype=int)
wr_id_count = 1

for usnpid, pouid, wrtype in unique_snp_ids:
    row_id = (usnpid == wr_snp_ids[:,0]) & (pouid == wr_snp_ids[:,1])
    row_id = np.where(row_id == True)[0]
    for ri in row_id:
        if non_empty_acres[ri][4] == use_cd:
            
            if non_empty_acres_arr[ri][6] == '-': # Only true for Umatilla based on a manual review of the WRs records
                non_empty_acres_arr[ri][6] = '1985-01-01'
                
            prior_date = datetime.datetime.strptime(non_empty_acres_arr[ri][6], '%Y-%m-%d').date()
            wr_swat_file.append([wr_id_count, usnpid, pouid, swat_to_wrtype[non_empty_acres[ri][5]], prior_date,round(non_empty_acres[ri][14]),non_empty_acres[ri][16],non_empty_acres[ri][17]]) 
            wr_id_count = wr_id_count + 1

dates_list = [wr_swat_file[ri][4] for ri in range(len(wr_swat_file))]
prior_dates, sort_index = np.unique(dates_list,return_inverse = 1)

for ri in range(len(wr_swat_file)):
    wr_swat_file[ri][4] = sort_index[ri] + 1

#%%

csv_file = output_path + '/wr_swat_filev3.csv'
filein = open(csv_file,'w')

atxt = 'WR_ID, SNAP_ID, POU_ID, WR_TYPE, PRIOR, VOL, START_DATE, END_DATE'
filein.write(atxt + '\n')

for i in range(len(wr_swat_file)):
    c=0
    for val in wr_swat_file[i]:
        if c == 0:
            atxt = str(val)
            c = 1
        else:
            atxt = atxt + ',' + str(val)      
    filein.write(atxt + '\n') 
filein.close()


#%%
# Version 3 - Only WRs
temp_lu = np.ones(lu_region.shape)
temp_lu[lu_region == 0] = 0
lu_region = temp_lu
lu_raster_flat = lu_region.flatten()

## Version 4 - WRs + CDL fields

## Version 5 - WRs + Small fields from CDL
#uni_lnduse = np.unique(lu_region.flatten())
#uni_lnduse = uni_lnduse[uni_lnduse > 0]
#
#temp_lu = lu_region.copy()
#for lnduse in uni_lnduse:
#    if len(np.where(lu_region.flatten() == lnduse)[0]) > 700.0:
#        temp_lu[lu_region == lnduse] = 0
#
#lu_region = temp_lu
#lu_raster_flat = lu_region.flatten()

combo_arr = np.transpose(np.vstack((lu_raster_flat, swwr_flat, gwwr_flat)))

unique_rows = np.unique(combo_arr, axis=0)
no_wrs_data = np.where(np.sum(unique_rows[:,[1,2]],axis=1) == 0)[0]

#unique_rows = unique_rows[no_wrs_data != 0,:]

hru_wr_dict = []

new_hrus = np.zeros(lu_region.shape, dtype=np.float)

hru_counter = 1
per_complt = 0.

for i in range(len(unique_rows)):
    write_bool = 0
    if i not in no_wrs_data:
        temp_wrs = []
        twrc = 0
        if str(int(unique_rows[i,1])) in hru_wr['IRRIGATION']['SW'].keys():
            for snp_ids, pouid in  hru_wr['IRRIGATION']['SW'][str(int(unique_rows[i,1]))]:
                for ii in range(len(wr_swat_file)):
                    if wr_swat_file[ii][1] == snp_ids and wr_swat_file[ii][2] == pouid and wr_swat_file[ii][3] == 1:
                        #hru_wr_dict.append([hru_counter,wr_swat_file[ii][0],1])
                        temp_wrs.append([hru_counter,wr_swat_file[ii][0],wr_swat_file[ii][4],wr_swat_file[ii][3],twrc,0])
                        #twrc = twrc + 1
        
        if str(int(unique_rows[i,2])) in hru_wr['IRRIGATION']['GW'].keys():
            for snp_ids, pouid in  hru_wr['IRRIGATION']['GW'][str(int(unique_rows[i,2]))]:
                for ii in range(len(wr_swat_file)):
                    if wr_swat_file[ii][1] == snp_ids and wr_swat_file[ii][2] == pouid and wr_swat_file[ii][3] == 3:
                        #hru_wr_dict.append([hru_counter,wr_swat_file[ii][0],3])
                        temp_wrs.append([hru_counter,wr_swat_file[ii][0],wr_swat_file[ii][4],wr_swat_file[ii][3],twrc,0])
                        #twrc = twrc + 1
        
        temp_wrs = np.asarray(temp_wrs)
        utemp_wrs, temp_wrs_order, temp_wrs_counts = np.unique(temp_wrs[:,2], return_counts = 1, return_inverse = 1)
        
        if len(np.where(temp_wrs_counts > 1)[0]) == 0:
            temp_wrs[:,5] = temp_wrs_order + 1
        else:
            temp_wrs_order = np.zeros(np.shape(temp_wrs_order))
            counter = 1
            for j in range(len(utemp_wrs)):
                row_ids = np.where(temp_wrs[:,2] == utemp_wrs[j])[0]
               
                if temp_wrs_counts[j] == 1:
                    temp_wrs_order[row_ids] = counter
                    counter = counter + 1
                    
                else:
                    uwr_type, torder, tcounts = np.unique(temp_wrs[row_ids,3], return_counts = 1, return_inverse = 1)
                    for wrt in uwr_type:
                        trow_id = np.where(temp_wrs[row_ids,3] == wrt)[0]
                        for jj in trow_id:
                            temp_wrs_order[row_ids[jj]] = counter
                            counter = counter + 1
            
            temp_wrs[:,5] = temp_wrs_order
        
        
        ids_bool = (lu_region == unique_rows[i,0]) & (swwr_raster == unique_rows[i,1]) & (gwwr_raster == unique_rows[i,2])
        
        if len(temp_wrs) > 1 and np.sum(np.asarray(ids_bool.flatten(), dtype=int)) <= hru_wr_area_tresh*2.0:
            write_bool = 1
      
        if write_bool == 0:
            for r in range(len(temp_wrs)):
                hru_wr_dict.append([temp_wrs[r,0],temp_wrs[r,1],temp_wrs[r,2],temp_wrs[r,3],temp_wrs[r,4],temp_wrs[r,5]])
          
            # Create new raster of CDL + WRs
            new_hrus[ids_bool == True] = hru_counter
            
            hru_counter = hru_counter + 1
        
    else:
        
        if unique_rows[i,0] > 0:
            hru_wr_dict.append([hru_counter,9999,9999,0,0,0])
            ids_bool = (lu_region == unique_rows[i,0]) & (swwr_raster == unique_rows[i,1]) & (gwwr_raster == unique_rows[i,2])
            new_hrus[ids_bool == True] = hru_counter
            
            hru_counter = hru_counter + 1
        
    if (i/float(len(unique_rows)))*100.0 > per_complt:
        print 'Completed %: ' + str(per_complt)
        per_complt = per_complt + 10.0

    
new_hrus[new_hrus == 0.0] = gwwr_NoData
newRaster_name = output_path + '/UHRUs_WR_Region_Latest_Clipv3.tif'
print ('Writing new HRU raster')
QSWAT_utils.Save_NewRaster(new_hrus, gwwr_obj_temp, gwwr_raster, newRaster_name, gwwr_NoData)


#%%
#hru_wr_dict_arr = np.asarray(hru_wr_dict)
#prior_dates, sort_index = np.unique(hru_wr_dict_arr[:,2],return_inverse = 1)

#csv_file = output_path + '/hru_wr_swat_filev2.csv'
#filein = open(csv_file,'w')
#
#atxt = 'HRU_ID, WR_ID, PRIOR, HRU_PRIOR'
#filein.write(atxt + '\n')
#
#for i in range(len(hru_wr_dict)):
#    c=0
#    for val in hru_wr_dict[i]:
#        if c == 0:
#            atxt = str(val)
#        elif c != 4:
#            atxt = atxt + ',' + str(val)
#        c = c + 1
#            
#    filein.write(atxt + '\n') 
#filein.close()

#%%
#hru_wr_dict = np.asarray(hru_wr_dict)
#hru_wr_file = np.zeros((len(hru_wr_dict),4))
#for rid in range(len(hru_wr_dict)):
#    hru_wr_dict
#
#
#
##%%
#swwr_raster_bin = np.zeros(np.shape(swwr_raster))
#gwwr_raster_bin = np.zeros(np.shape(gwwr_raster))
#
#swwr_raster_bin[swwr_raster > 0] = 1
#gwwr_raster_bin[gwwr_raster > 0] = 1
#
#total_wr_area = np.sum(swwr_raster_bin.flatten()) + np.sum(gwwr_raster_bin.flatten())
#
#dual_wr_area = swwr_raster_bin.flatten() + gwwr_raster_bin.flatten()
#dual_wr_area = len(dual_wr_area[dual_wr_area > 1])
#
#print str(round((dual_wr_area/total_wr_area)*100)) + '%'


#%%
#outShapefile = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT\Umatilla_Meghna_Test\Watershed\Shapes\hru1.shp'    
#
#output_path = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\HRUs_Meghna\Zonal_Stats\WR_Tiffs' 
#raster_file = output_path + '/WR_POU_SW_Clip.tif'
#BoundaryRaster = output_path + '/temp_boundary.tif'
#
#QSWAT_utils.CreateTiff(outShapefile, BoundaryRaster, raster_file,'HRUGIS','','','tiff',gdal.GDT_Float64)
#hru_raster, hru_NoData, hru_obj_temp = QSWAT_utils.Read_Raster(BoundaryRaster)
#hru_raster = np.asarray(np.flip(hru_raster,0),dtype=float)
#
#newRaster_name = output_path + '/UHRUs_WR_Region_Latest_Clipv4.tif'
#region_temp, region_NoData_temp, region_obj_temp = QSWAT_utils.Read_Raster(newRaster_name)
##region_temp = np.asarray(np.flip(region_temp,0),dtype=float)
#
#if region_NoData_temp != None:
#    region_temp[region_temp == region_NoData_temp] = 0.0
#
##%%
#hru_flat = hru_raster.flatten()
#uhrus, uhrus_counts = np.unique(hru_flat, return_counts = 1)
#uhrus = uhrus[uhrus > 0]
#
#region_flat = region_temp.flatten()
#
#combo_arr = np.transpose(np.vstack((region_flat, hru_flat)))
#
#unique_rows = np.unique(combo_arr, axis=0)
##unique_rows, rows_counts = np.unique(combo_arr, axis=0, return_counts = 1)
##rows_counts = rows_counts[np.where(unique_rows[:,1] > 0)[0],:] 
#unique_rows = unique_rows[np.where(unique_rows[:,1] > 0)[0],:] 
#unique_rows = unique_rows[np.where(unique_rows[:,0] > 0)[0],:] 
#unique_hrus = np.unique(unique_rows[:,1])
#unique_wr = np.unique(unique_rows[:,0])
#
##for uh in uhrus:
##    uwrs = np.unique(region_flat[hru_flat == uh])
##    if len(uwrs) > 1:
##        print uh, uwrs
##        break

#%%
#fnames = os.listdir(cdl_path)  
#no_crop_ids = Get_No_Crop_Ids()
#    
#raster_file = cdl_path + '\\' + cdl_file
#region, region_NoData, region_obj = QSWAT_utils.Read_Raster(raster_file)
#region_prj = region_obj.GetProjection()
#
##centroids_X, centroids_Y = GetPixelCentroids(region_obj)
##
##top_row, last_row = QSWAT_utils.Raster_row_boundaries(region)
##left_col, right_col = QSWAT_utils.Raster_col_boundaries(region)
#
#output_path = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\HRUs_Meghna\Zonal_Stats\WR_Tiffs'
#
#BoundaryRaster = output_path + '/hrus.tif'
#outShapefile = root_path + '\Correct_HRUs_11534_ExemptLU_333.shp'
#QSWAT_utils.CreateTiff(outShapefile, BoundaryRaster, raster_file,'HRUGIS','','','tiff',gdal.GDT_Float64)
#hru_raster, hru_NoData, hru_obj_temp = QSWAT_utils.Read_Raster(BoundaryRaster)
#hru_raster = np.asarray(np.flip(hru_raster,0),dtype=float)
#
##plt.matshow(hru_raster)
#if hru_NoData != None:
#    hru_raster[hru_raster == hru_NoData] = 0.0
#
##hru_flat = hru_raster.flatten()
##uhru = np.unique(hru_raster.flatten())
##uhru = uhru(np.where(uhru!=0))
#
#print 'Rasterized HRUs'
#
##new_wrs_vol = dict()
##total_wr_vol = dict()
#
#hru_wr = dict()
#
#clayer = 0
##for use_code in pou_uses.keys():
#for use_code in temp_use:
#    if use_code not in hru_wr.keys():
#        hru_wr[use_code] = dict()
##        new_wrs_vol[use_code] = dict()
##        total_wr_vol[use_code] = dict()
#        
#    for wr_type in pou_uses[use_code].keys():
#        if wr_type not in hru_wr[use_code].keys():
#            hru_wr[use_code][wr_type] = dict()
##            new_wrs_vol[use_code][wr_type] = []
##            total_wr_vol[use_code][wr_type] = 0
#            
#        new_wrs = np.zeros(region.shape, dtype=np.int64)
#        new_wrs_duty = np.zeros(region.shape, dtype=np.int64)
#        
#        outFileName = '\wr_v_pou_public_proj_REGION_' + str(clayer) + '_' + str(wr_type) + '.shp'
#        #boundary_lyr = outPath + '\wr_v_pou_public_proj_REGION_0_GW.shp'
#        
#        boundary_lyr = outPath + outFileName
#        
#        
#        dataSource = driver.Open(boundary_lyr, 0)
#        temp_layer = dataSource.GetLayer()
#        
#        for feature in temp_layer:
#            outFileName = '/temp.shp'
#            outShapefile = output_path + outFileName
#            outDriver = ogr.GetDriverByName("ESRI Shapefile")
#        
#            # Remove output shapefile if it already exists
#            if os.path.exists(outShapefile):
#                outDriver.DeleteDataSource(outShapefile)
#                
#            proj_outFileName = '/temp.prj'
#            with open(output_path + proj_outFileName, 'w') as prj_file:
#                #prj_file.write(str(spatialRef.ExportToWkt()))
#                prj_file.write(str(region_prj))
#            prj_file.close()
#        
#        
#            # Create the output shapefile
#            outDataSource = outDriver.CreateDataSource(outShapefile)
#            outLayer = outDataSource.CreateLayer("temp", geom_type = ogr.wkbPolygon)
#        
#            # Add an ID field
#            idField = ogr.FieldDefn("id", ogr.OFTInteger)
#            outLayer.CreateField(idField)
#            
#            featureDefn = outLayer.GetLayerDefn()
#            featureo = ogr.Feature(featureDefn)
#                
#            featureo.SetField("id", 1)
#            points = feature.GetGeometryRef()
#            featureo.SetGeometry(points)
#            
#            snp_id = feature.GetField("snp_id")
#            wr_vol = feature.GetField("wr_vol")
#            wr_duty = feature.GetField("duty")
#            
##            total_wr_vol[use_code][wr_type] = total_wr_vol[use_code][wr_type] + feature.GetField("wr_vol")
#                            
#            outLayer.CreateFeature(featureo)
#            outDataSource = None
#            
#            BoundaryRaster = output_path + '/temp_boundary.tif'
#            QSWAT_utils.CreateTiff(outShapefile, BoundaryRaster, raster_file,'',1,'','tiff')
#            
#            region_temp, region_NoData_temp, region_obj_temp = QSWAT_utils.Read_Raster(BoundaryRaster)
#            
#            if len(np.where(region_temp > 0)) == 0:
#                break
#            
#            region_temp = np.asarray(np.flip(region_temp,0),dtype=float)
#            
#            if region_NoData_temp != None:
#                region_temp[region_temp == region_NoData_temp] = 0.0
#            
#            uhru = np.unique(hru_raster[region_temp > 0.0])
#            uhru = uhru[uhru != 0]
#            
#            if len(uhru) > 0:
#                #if snp_id not in hru_wr[use_code][wr_type].keys():
#                #    hru_wr[use_code][wr_type][snp_id] = []
#                for u in uhru:
#                    if u not in hru_wr[use_code][wr_type].keys():
#                        hru_wr[use_code][wr_type][int(u)] = []
#                    #hru_wr[use_code][wr_type][int(u)].append(snp_id)
#                    hru_wr[use_code][wr_type][int(u)].append(wr_vol)
#            
##            new_wrs_vol[use_code][wr_type].append([snp_id,np.sum(region_temp.flatten())])
#            
#            new_wrs = new_wrs + region_temp
##            new_wrs_duty = new_wrs_duty + region_temp
#            
#            
##            new_wrs_duty = new_wrs_duty + region_temp
#            
##            id_raster = np.where(new_wrs_duty.flatten() > 1)
##            if len(id_raster[0]) > 0:
##                new_wrs_flat = new_wrs_duty.flatten()
##                for idr in id_raster[0]:
##                    if idr not in new_wrs_vol[use_code][wr_type].keys():
##                        new_wrs_vol[use_code][wr_type][idr] = []
##                    
##                    new_wrs_vol[use_code][wr_type][idr].append(wr_vol)
##            
##            new_wrs_duty[new_wrs_duty > 1.0] = 1.0
#                    
#        outRasterName = 'WR_POU_proj_REGION_' + str(clayer) + '_' + str(wr_type) + '.tif'
#        print outRasterName
#        newRaster_name = output_path + '/' + outRasterName
#        QSWAT_utils.Save_NewRaster(new_wrs, region_obj, region, newRaster_name, -999.0)
#        
#    clayer = clayer + 1
#
###ft_idc = 0
###poly_area = dict()
###pts_len =[]
###
###for feature in layer:
###    geom = feature.geometry()
###    points = feature.GetGeometryRef()
###    ft_id = feature.GetField("OBJECTID")
###    snp_id = feature.GetField("snp_id")
###    
###    # Create the feature and set values
###    featureDefn = outLayer.GetLayerDefn()
###    featureo = ogr.Feature(featureDefn)
###    
###    temp_area = []
###    diam_bool = 0
###    for pt in points:
####        if pt.GetArea() > 220.0 and pt.GetArea() < 225.0:
####            diam_bool = 1
###        temp_area.append(pt.GetArea())
###        
###    if np.mean(temp_area) > 220.0 and np.mean(temp_area) < 225.0:
###        diam_bool = 1
###        
###    pts = []
###    area = []
###    cpoly = 1
###    if diam_bool == 1:
###        for pt in points:
###    
###            featureo.SetField("id", ft_idc)
###            featureo.SetField("obj_wrid", ft_id)
###            featureo.SetField("snp_id", snp_id)
###        
###            #diam_bool = 1
###            bufferDistance = 500
###            
###            if not pt.GetPoints():
###                ptsc = pt.GetGeometryRef(0)
###                ptsc = np.asarray(ptsc.GetPoints())
###            else:
###                ptsc = np.asarray(pt.GetPoints())
###            
###            
###            xc, yc = centeroidnp(ptsc)
###            wkt = 'POINT (' + str(xc) + ' ' + str(yc) + ')'
###            ptsc = ogr.CreateGeometryFromWkt(wkt)
###            poly = ptsc.Buffer(bufferDistance)
###            ##print poly
###            ##poly = poly.GetGeometryRef(0)
###            
###            #featureo.SetGeometry(poly)
###            
###            ##pts.append(poly.GetPoints())
###            ##area.append(pt.GetArea())
###            #outLayer.CreateFeature(featureo)
###            #ft_idc = ft_idc + 1
###            
###            
###            if cpoly == 1:
###                poly2 = poly
###                #poly2 = ogr.Geometry(ogr.wkbMultiPolygon)
###                cpoly = 0
###            else:
###                poly2 = poly2.Union(poly)
###            
###            #poly2.AddGeometry(poly)
###                
###        
###        featureo.SetGeometry(poly2)
###        
###        #convexhull = poly2.ConvexHull()
###        #featureo.SetGeometry(convexhull)
###        
###        outLayer.CreateFeature(featureo)
###        ft_idc = ft_idc + 1