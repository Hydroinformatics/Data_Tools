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


import owrd_utils as owrd_utils

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
    
    # Requested by Meghna
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

feature_dict, remove_features = owrd_utils.EliminateDuplicates(root_path, pou_file, pou_file_noduplicate)

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
#
#temp_total_wr_vol = dict()
#missed_total_wr_vol = dict()
#
##for use_code in pou_uses.keys():
#for use_code in temp_use:
#    if use_code not in hru_wr.keys():
#        hru_wr[use_code] = dict()
#        temp_total_wr_vol[use_code] = dict()
#        missed_total_wr_vol[use_code] = dict()
#        
#    for wr_type in pou_uses[use_code].keys():
#    #for wr_type in temp_wr_type:
#        if wr_type not in hru_wr[use_code].keys():
#            hru_wr[use_code][wr_type] = dict()
#            temp_total_wr_vol[use_code][wr_type] = 0
#            missed_total_wr_vol[use_code][wr_type] = 0
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
#            temp_total_wr_vol[use_code][wr_type] = temp_total_wr_vol[use_code][wr_type] + feature.GetField("wr_vol")
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
#                if labeled_area[bi] <= hru_wr_area_tresh and bi != np.argmax(labeled_area): # version < 5
#                #if labeled_area[bi] <= hru_wr_area_tresh*10: # version >= 5
#                    region_temp[result == bi+1] = 0
#        
#            #print cfeat
#            #if cfeat == 1:
#            if cfeat == 1 and np.sum(region_temp.flatten()) > 0:
#                new_wrs[region_temp > 0.0] = cfeat
#                hru_wr[use_code][wr_type][cfeat] = [(snp_id, pou_id)]
#                cfeat = cfeat + 1
#            
#            #else:
#            elif np.sum(region_temp.flatten()) > 0:
#                print cfeat
##                if cfeat == 173:
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
#            else:
#                missed_total_wr_vol[use_code][wr_type] = missed_total_wr_vol[use_code][wr_type] + feature.GetField("wr_vol")
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
            wr_swat_file.append([wr_id_count, usnpid, pouid, swat_to_wrtype[non_empty_acres[ri][5]], prior_date,round(non_empty_acres[ri][14]),non_empty_acres[ri][16],non_empty_acres[ri][17],
                                 non_empty_acres[ri][17] - non_empty_acres[ri][16] + 1]) 
            wr_id_count = wr_id_count + 1

dates_list = [wr_swat_file[ri][4] for ri in range(len(wr_swat_file))]
prior_dates, sort_index = np.unique(dates_list,return_inverse = 1)

for ri in range(len(wr_swat_file)):
    wr_swat_file[ri][4] = sort_index[ri] + 1

#No WR data
wr_swat_file.append([int('9'*len(str(ri))), int('9'*len(str(ri))), int('9'*len(str(ri))), 0, 0, 0, 1, 365, 365])

#%%

csv_file = output_path + '/Final/wr_swat_file_v3.csv'
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


num_year_sim = 12
csv_file = output_path + '/Final/wrdata.dat'
filein = open(csv_file,'w')

for yr in range(0,num_year_sim):
    for i in range(len(wr_swat_file)):
        atxt = str(yr+1).rjust(4) + ''.rjust(3)
        atxt = atxt + str(wr_swat_file[i][0]).rjust(5) + ''.rjust(3)
        atxt = atxt + str(wr_swat_file[i][3]).rjust(4)  + ''.rjust(3)
        atxt = atxt + str(int(wr_swat_file[i][5])).rjust(6)  + ''.rjust(3)
        atxt = atxt + str(wr_swat_file[i][6]).rjust(4)  + ''.rjust(3)
        atxt = atxt + str(wr_swat_file[i][7]).rjust(4)
        #if yr == num_year_sim-1 and i == len(wr_swat_file)
        filein.write(atxt + '\n') 
filein.close()


#%%
# Version 3 - Only WRs
#temp_lu = np.zeros(lu_region.shape)
#temp_lu[lu_region == 0] = 0
#lu_region = temp_lu
#lu_raster_flat = lu_region.flatten()

## Version 4 - WRs + CDL fields

## Version 5 - WRs + Small fields from CDL
uni_lnduse = np.unique(lu_region.flatten())
uni_lnduse = uni_lnduse[uni_lnduse > 0]

temp_lu = lu_region.copy()
for lnduse in uni_lnduse:
    if len(np.where(lu_region.flatten() == lnduse)[0]) > 700.0:
        temp_lu[lu_region == lnduse] = 0

lu_region = temp_lu
lu_raster_flat = lu_region.flatten()

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
                        temp_wrs.append([hru_counter,wr_swat_file[ii][0],wr_swat_file[ii][4],wr_swat_file[ii][6],wr_swat_file[ii][8],wr_swat_file[ii][3],twrc,0])
                        #twrc = twrc + 1
        
        if str(int(unique_rows[i,2])) in hru_wr['IRRIGATION']['GW'].keys():
            for snp_ids, pouid in  hru_wr['IRRIGATION']['GW'][str(int(unique_rows[i,2]))]:
                for ii in range(len(wr_swat_file)):
                    if wr_swat_file[ii][1] == snp_ids and wr_swat_file[ii][2] == pouid and wr_swat_file[ii][3] == 3:
                        #hru_wr_dict.append([hru_counter,wr_swat_file[ii][0],3])
                        #temp_wrs.append([hru_counter,wr_swat_file[ii][0],wr_swat_file[ii][4],wr_swat_file[ii][3],twrc,0])
                        temp_wrs.append([hru_counter,wr_swat_file[ii][0],wr_swat_file[ii][4],wr_swat_file[ii][6],wr_swat_file[ii][8],wr_swat_file[ii][3],twrc,0])
                        #twrc = twrc + 1
        
        #if len(temp_wrs) > 2:
        #    print i
        
        temp_wrs = np.asarray(sorted(np.asarray(temp_wrs), key=lambda x: (-x[4], x[3], x[2], x[5])))
        temp_wrs_order = 1
        for i in range(len(temp_wrs)):
            temp_wrs[i,6] = temp_wrs_order
            temp_wrs_order = temp_wrs_order + 1
        
#        temp_wrs = np.asarray(temp_wrs)
#        utemp_wrs, temp_wrs_order, temp_wrs_counts = np.unique(temp_wrs[:,2], return_counts = 1, return_inverse = 1)
#        
#        if len(np.where(temp_wrs_counts > 1)[0]) == 0:
#            temp_wrs[:,7] = temp_wrs_order + 1
#        else:
#            temp_wrs_order = np.zeros(np.shape(temp_wrs_order))
#            counter = 1
#            for j in range(len(utemp_wrs)):
#                row_ids = np.where(temp_wrs[:,2] == utemp_wrs[j])[0]
#               
#                if temp_wrs_counts[j] == 1:
#                    temp_wrs_order[row_ids] = counter
#                    counter = counter + 1
#                    
#                else:
#                    uwr_type, torder, tcounts = np.unique(temp_wrs[row_ids,3], return_counts = 1, return_inverse = 1)
#                    for wrt in uwr_type:
#                        trow_id = np.where(temp_wrs[row_ids,3] == wrt)[0]
#                        for jj in trow_id:
#                            temp_wrs_order[row_ids[jj]] = counter
#                            counter = counter + 1
#            
#            temp_wrs[:,5] = temp_wrs_order
        
        
        ids_bool = (lu_region == unique_rows[i,0]) & (swwr_raster == unique_rows[i,1]) & (gwwr_raster == unique_rows[i,2])
        
        if len(temp_wrs) > 1 and np.sum(np.asarray(ids_bool.flatten(), dtype=int)) <= hru_wr_area_tresh*2.0:
            write_bool = 1
      
        if write_bool == 0:
            for r in range(len(temp_wrs)):
                #hru_wr_dict.append([temp_wrs[r,0],temp_wrs[r,1],temp_wrs[r,2],temp_wrs[r,3],temp_wrs[r,4],temp_wrs[r,5]])
                hru_wr_dict.append([temp_wrs[r,0],temp_wrs[r,1],temp_wrs[r,2],temp_wrs[r,3],temp_wrs[r,4],temp_wrs[r,5],temp_wrs[r,6]])
                
            # Create new raster of CDL + WRs
            new_hrus[ids_bool == True] = hru_counter
            
            hru_counter = hru_counter + 1
        
    else:
        
        if unique_rows[i,0] > 0:
            no_wr_id = wr_swat_file[len(wr_swat_file)-1][0]
            hru_wr_dict.append([hru_counter,no_wr_id,no_wr_id,0,0,0,1])
            ids_bool = (lu_region == unique_rows[i,0]) & (swwr_raster == unique_rows[i,1]) & (gwwr_raster == unique_rows[i,2])
            new_hrus[ids_bool == True] = hru_counter
            
            hru_counter = hru_counter + 1
        
    if (i/float(len(unique_rows)))*100.0 > per_complt:
        print 'Completed %: ' + str(per_complt)
        per_complt = per_complt + 10.0

    
#new_hrus[new_hrus == 0.0] = gwwr_NoData
#newRaster_name = output_path + '/UHRUs_WR_Region_Latest_Clip_V3.tif'
#print ('Writing new HRU raster')
#QSWAT_utils.Save_NewRaster(new_hrus, gwwr_obj_temp, gwwr_raster, newRaster_name, gwwr_NoData)

#%%
import pyodbc
db_path = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT\Umatilla_InterACTWEL_QSWATv4\Umatilla_InterACTWEL_QSWATv4\Umatilla_InterACTWEL_QSWATv4.mdb'    
conn_str = (
    r'DRIVER={Microsoft Access Driver (*.mdb, *.accdb)};'
    r'DBQ=' + db_path + ';'
    )
cnxn = pyodbc.connect(conn_str)
crsr = cnxn.cursor()

crsr.execute('select * from hrus')
swat_hruid = dict()
for row in crsr.fetchall():
    swat_hruid[int(row[12])] = row[11]

#%%
        
outShapefile = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT\Umatilla_InterACTWEL_QSWATv4\Umatilla_InterACTWEL_QSWATv4\Watershed\Shapes\hru2.shp'    

output_path = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\HRUs_Meghna\Zonal_Stats\WR_Tiffs' 
raster_file = output_path + '/WR_POU_SW_Clip.tif'
BoundaryRaster = output_path + '/temp_boundary.tif'

QSWAT_utils.CreateTiff(outShapefile, BoundaryRaster, raster_file,'HRUGIS','','','tiff',gdal.GDT_Float64)
hru_raster, hru_NoData, hru_obj_temp = QSWAT_utils.Read_Raster(BoundaryRaster)
hru_raster = np.asarray(np.flip(hru_raster,0),dtype=float)

wrRaster_name = output_path + '/UHRUs_WR_Region_Latest_Clip_V3.tif'
wrregion_temp, wrregion_NoData_temp, wrregion_obj_temp = QSWAT_utils.Read_Raster(wrRaster_name)
#region_temp = np.asarray(np.flip(region_temp,0),dtype=float)

if wrregion_NoData_temp != None:
    wrregion_temp[wrregion_temp == wrregion_NoData_temp] = 0.0

hru_flat = hru_raster.flatten()
uhrus, uhrus_counts = np.unique(hru_flat, return_counts = 1)
uhrus = uhrus[uhrus > 0]

wr_region_flat = wrregion_temp.flatten()

combo_arr = np.transpose(np.vstack((wr_region_flat, hru_flat)))

unique_rows = np.unique(combo_arr, axis=0)
#unique_rows, rows_counts = np.unique(combo_arr, axis=0, return_counts = 1)
#rows_counts = rows_counts[np.where(unique_rows[:,1] > 0)[0],:] 
unique_rows = unique_rows[np.where(unique_rows[:,1] > 0)[0],:] 
unique_rows = unique_rows[np.where(unique_rows[:,0] > 0)[0],:]
 
unique_hrus = np.unique(unique_rows[:,1])
unique_wr = np.unique(unique_rows[:,0])

#for uh in uhrus:
#    uwrs = np.unique(region_flat[hru_flat == uh])
#    if len(uwrs) > 1:
#        print uh, uwrs
#        break

#%%
hru_wr_dict_arr = np.asarray(hru_wr_dict)

hrwwr_dat = []
for wrid, hruid in unique_rows:
    temp_ids = np.where(wrid == hru_wr_dict_arr[:,0])[0]
    for tid in temp_ids:
        hrwwr_dat.append(np.hstack((swat_hruid[int(hruid)],hru_wr_dict_arr[tid,:])))
        
miss_uhrus = np.setdiff1d(uhrus,unique_rows[:,1])
for hruid in miss_uhrus:
    hrwwr_dat.append(np.hstack((swat_hruid[int(hruid)],[0,9999,9999,0,0,0,1])))
        

#hrwwr_dat_arr = np.asarray(hrwwr_dat)
#%%

#hru_wr_dict_arr = np.asarray(sorted(np.asarray(hru_wr_dict), key=lambda x: (x[6], x[2])))
#hru_wr_dict_arr = np.asarray(sorted(np.asarray(hru_wr_dict), key=lambda x: (x[2])))
hrwwr_dat_arr = np.asarray(sorted(np.asarray(hrwwr_dat), key=lambda x: (x[3])))

#prior_dates, sort_index = np.unique(hru_wr_dict_arr[:,2],return_inverse = 1)

#csv_file = output_path + '/hruwr_v2.csv'
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

#data = [2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,20,21,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,59,60,61,63,64,65,66,75,76,77,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,146,147,152,153,154,155,156,157,158,159,160,161,162,171,172,173,174,175,176,177,178,179,180,181,189,190,191,196,197,198,202,203,204,205,206,207,208,209,210,211,212,213,214,220,221,222,223,224,234,256,257,258,259,260,261,262,263,264,265,281,282,283,284,285,286,287,288,289,291,292,293,294,295,296,297,298,299,300,301,302,312,313,314,315,316,318,319,320,321,322,323,352,353,354,355,357,358,359,360,361,362,363,364,365,366,367,368,369,370,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,396,397,398,399,400,401,402,407,408,409,410,411,421,422,423,424,425,426,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,493,494,495,496,497,498,499,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,519,523,524,525,526,527,528,529,530,531,532,533,534,535,536,537,538,539,540,541,542,543,544,545,546,547,553,554,555,556,557,558,559,560,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,604,605,606,607,608,609,610,611,612,613,614,615,616,618,619,620,626,627,628,629,630,631,642,643,644,645,646,647,648,649,650,651,652,653,654,655,656,657,658,659,660,661,662,663,664,665,667,668,669,670,671,672,673,674,675,676,681,682,683,684,685,686,687,688,689,690,691,692,693,694,695,698,699,701,702,703,704,705,706,707,708,709,710,711,712,713,716,717,718,719,720,721,722,725,726,748,749,750,808,809,810,811,812,860,861,862,863,864,866,867,868,869,948,949,950,958,959,960,961,962,963,965,966,967,968,982,983,984,993,994,995,999,1000,1001,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,1022,1088,1089,1090,1091,1092,1093,1094,1105,1106,1107,1108,1109,1114,1115,1116,1117,1118,1119,1120,1121,1122,1123,1124,1125,1126,1135,1136,1137,1138,1139,1140,1141,1142,1143,1144,1145,1146,1149,1150,1151,1152,1153,1154,1167,1168,1169,1170,1171,1172,1173,1174,1175,1176,1197,1198,1199,1200,1201,1202,1203,1204,1205,1206,1207,1212,1213,1214,1215,1216,1217,1221,1222,1223,1224,1225,1235,1236,1238,1239,1240,1241,1242,1243,1244,1245,1246,1247,1248,1249,1250,1251,1252,1253,1254,1255,1256,1259,1260,1261,1262,1263,1264,1265,1266,1267,1269,1270,1271,1272,1273,1274,1275,1276,1277,1282,1283,1284,1285,1286,1287,1288,1289,1290,1291,1292,1293,1294,1295,1296,1297,1298,1299,1306,1307,1308,1309,1310,1312,1313,1314,1315,1316,1317,1318,1319,1323,1324,1327,1336,1337,1341,1342,1343,1344,1345,1346,1347,1348,1349,1350,1351,1352,1353,1365,1366,1367,1368,1369,1370,1371,1372,1373,1374,1375,1376,1377,1378,1379,1380,1381,1382,1383,1384,1385,1386,1387,1388,1389,1390,1391,1392,1395,1396,1397,1398,1399,1400,1401,1402,1403,1404,1405,1406,1407,1408,1411,1412,1414,1415,1416,1417,1418,1419,1420,1421,1422,1423,1425,1426,1427,1428,1429,1430,1431,1432,1433,1434,1435,1457,1458,1459,1460,1461,1473,1474,1475,1476,1497,1536,1537,1538,1539,1540,1541,1542,1543,1544,1545,1546,1556,1557,1558,1559,1560,1561,1562,1563,1564,1565,1566,1567,1590,1591,1592,1593,1594,1595,1596,1597,1598,1599,1600,1606,1607,1608,1609,1614,1615,1628,1629,1630,1631,1632,1633,1643,1644,1645,1646,1647,1648,1649,1650,1677,1678,1679,1680,1681,1707,1708,1709,1710,1711,1712,1713,1714,1715,1716,1717,1718,1724,1725,1726,1727,1728,1729,1730,1731,1732,1756,1757,1758,1759,1760,1761,1762,1763,1764,1765,1766,1767,1768,1769,1770,1771,1772,1776,1777,1846,1847,1848,1849,1850,1851,1852,1854,1855,1856,1857,1859,1860,1861,1862,1863,1864,1865,1866,1867,1868,1869,1870,1871,1872,1873,1874,1875,1876,1877,1878,1879,1880,1881,1882,1883,1884,1885,1886,1887,1888,1889,1890,1891,1892,1893,1894,1895,1896,1897,1898,1899,1901,1902,1903,1904,1905,1906,1907,1908,1912,1913,1914,1915,1916,1928,1929,1930,1931,1932,1933,1934,1935,1936,1937,1938,2072,2073,2074,2075,2076,2082,2083,2084,2085,2086,2087,2088,2089,2090,2091,2092,2093,2094,2095,2096,2097,2098,2099,2108,2109,2110,2111,2112,2113,2114,2115,2116,2117,2118,2137,2138,2139,2140,2141,2142,2143,2144,2145,2146,2147,2209,2210,2211,2212,2213,2214,2215,2216,2217,2218,2219,2220,2221,2222,2223,2224,2225,2226,2227,2228,2229,2230,2231,2232,2233,2234,2235,2236,2273,2274,2275,2276,2277,2278,2279,2280,2281,2282,2283,2284,2285,2286,2287,2288,2289,2290,2291,2292,2293,2294,2295,2298,2299,2300,2301,2302,2303,2304,2305,2306,2324,2325,2326,2327,2328,2329,2330,2331,2332,2334,2335,2336,2337,2338,2339,2340,2341,2342,2343,2344,2355,2377,2378,2379,2380,2381,2382,2383,2384,2385,2386,2387,2388,2389,2390,2391,2392,2393,2395,2396,2397,2398,2399,2400,2401,2402,2403,2408,2409,2410,2411,2412,2413,2414,2415,2416,2426,2427,2428,2429,2495,2496,2497,2498,2499,2500,2501,2502,2503,2504,2505,2506,2507,2508,2509,2518,2519,2520,2521,2522,2523,2524,2525,2526,2527,2528,2529,2530,2560,2561,2562,2563,2564,2565,2566,2567,2568,2569,2570,2571,2593,2594,2595,2596,2597,2604,2674,2675,2676,2677,2678,2679,2680,2681,2682,2683,2684,2685,2686,2687,2688,2689,2690,2691,2692,2693,2694,2695,2696,2697,2698,2699,2700,2701,2702,2703,2704,2705,2706,2707,2708,2709,2710,2711,2712,2713,2718,2719,2720,2721,2722,2723,2724,2728,2729,2730,2731,2732,2733,2734,2735,2736,2737,2738,2740,2741,2742,2785,2786,2787,2788,2789,2790,2791,2792,2793,2794,2795,2810,2820,2821,2822,2823,2824,2825,2826,2827,2828,2829,2830,2855,2856,2857,2886,2887,2888,2889,2890,2891,2892,2893,2894,2895,2896,2906,2907,2908,2915,2916,2917,2918,2919,2920,2921,2922,2923,2924,2925,2926,2927,2928,2929,2931,2932,2933,2934,2935,2936,2937,2938,2959,2963,2964,2965,2966,2967,2968,2969,2970,2971,2972,2973,2974,2975,2976,2977,2978,2979,2984,2985,2986,2987,2988,2989,2990,2991,2992,2993,2994,2995,2996,2997,2998,2999,3000,3013,3014,3015,3016,3017,3018,3019,3020,3021,3027,3028,3029,3044,3045,3046,3047,3048,3049,3050,3081,3082,3083,3084,3085,3086,3087,3088,3089,3090,3091,3092,3093,3094,3095,3096,3097,3098,3099,3100,3101,3102,3103,3104,3105,3106,3107,3108,3111,3112,3113,3114,3115,3116,3117,3118,3119,3120,3121,3122,3123,3124,3125,3138,3139,3213,3214,3217,3218,3219,3220,3221,3222,3223,3224,3225,3226,3229,3230,3231,3232,3233,3234,3235,3236,3237,3238,3239,3240,3241,3242,3246,3247,3250,3251,3252,3253,3254,3256,3257,3258,3259,3260,3261,3264,3265,3266,3267,3268,3269,3270,3271,3272,3273,3274,3275,3276,3277,3278,3279,3280,3281,3282,3283,3284,3285,3286,3287,3288,3289,3290,3298,3299,3300,3301,3302,3303,3304,3305,3306,3307,3308,3316,3317,3318,3319,3320,3321,3335,3336,3340,3341,3342,3343,3364,3365,3366,3367,3368,3369,3370,3371,3372,3373,3374,3379,3380,3381,3382,3383,3384,3385,3389,3390,3391,3392,3393,3394,3395,3396,3397,3398,3399,3400,3401,3417,3418,3419,3420,3421,3422,3423,3424,3425,3426,3427,3428,3429,3430,3431,3432,3433,3434,3435,3438,3439,3440,3441,3442,3443,3452,3453,3454,3455,3456,3457,3459,3460,3461,3462,3466,3467,3468,3469,3470,3471,3472,3473,3474,3475,3490,3491,3492,3497,3498,3499,3500,3501,3502,3503,3504,3505,3506,3507,3508,3520,3521,3522,3523,3524,3525,3526,3527,3528,3529,3530,3531,3532,3534,3535,3536,3537,3538,3539,3545,3546,3547,3548,3549,3550,3551,3552,3553,3555,3556,3557,3558,3559,3560,3561,3562,3563,3564,3565,3566,3571,3572,3573,3574,3575,3576,3598,3599,3600,3601,3602,3603,3604,3605,3606,3607,3608,3609,3610,3611,3612,3613,3614,3615,3616,3617,3618,3619,3620,3621,3622,3623,3624,3625,3626,3627,3628,3629,3630,3631,3632,3633,3634,3636,3637,3638,3639,3650,3651,3652,3653,3654,3655,3656,3657,3658,3659,3660,3661,3662,3663,3664,3665,3666,3667,3668,3669,3674,3675,3676,3677,3678,3679,3680,3681,3682,3683,3684,3685,3686,3687,3688,3689,3690,3691,3692,3694,3695,3696,3699,3700,3701,3702,3703,3706,3707,3708,3709,3710,3711,3712,3713,3714,3715,3716,3717,3719,3720,3721,3722,3723,3724,3725,3728,3729,3730,3731,3732,3733,3734,3735,3736,3740,3741,3742,3743,3744,3749,3750,3751,3752,3753,3754,3755,3756,3757,3759,3760,3761,3762,3763,3764,3770,3771,3772,3773,3774,3775,3776,3777,3778,3779,3780,3781,3782,3783,3785,3786,3787,3788,3789,3790,3792,3793,3794,3795,3796,3797,3798,3799,3800,3801,3802,3803,3804,3805,3806,3807,3808,3809,3810,3811,3812,3814,3815,3816,3819,3820,3821,3822,3823,3824,3826,3827,3828,3829,3830,3831,3832,3833,3834,3835,3836,3837,3838,3839,3840,3841,3842,3843,3844,3845,3846,3848,3849,3851,3852,3853,3854,3855,3868,3869,3870,3871,3872,3873,3874,3875,3880,3881,3882,3883,3885,3886,3887,3888,3889,3890,3891,3917,3918,3919,3920,3923,3924,3925,3926,3927,3928,3929,3930,3931,3932,3933,3934,3935,3938,3939,3940,3941,3942,3943,3944,3951,3952,3953,3954,3955,3956,3957,3958,3959,3960,3961,3962,4000,4001,4002,4003,4004,4005,4006,4007,4008,4009,4010,4011,4012,4013,4021,4022,4023,4024,4025,4026,4027,4028,4029,4030,4031,4032,4033,4048,4049,4050,4051,4052,4053,4055,4056,4057,4058,4059,4060,4061,4062,4063,4064,4065,4066,4067,4068,4108,4109,4110,4146,4147,4148,4149,4150,4154,4155,4156,4157,4158,4159,4160,4161,4162,4163,4164,4165,4166,4168,4169,4170,4171,4172,4173,4174,4175,4176,4177,4178,4179,4180,4181,4182,4183,4184,4185,4186,4187,4188,4189,4190,4191,4192,4193,4194,4195,4196,4197,4198,4199,4200,4201,4202,4203,4204,4205,4206,4207,4208,4209,4211,4212,4213,4214,4215,4216,4217,4219,4220,4221,4222,4223,4224,4225,4226,4234,4235,4236,4237,4238,4239,4240,4241,4242,4243,4245,4246,4247,4248,4249,4250,4251,4252,4260,4261,4262,4263,4273,4274,4275,4276,4277,4278,4279,4280,4281,4282,4283,4284,4285,4286,4287,4288,4289,4290,4291,4292,4325,4326,4334,4335,4336,4337,4338,4339,4340,4341,4342,4343,4344,4345,4346,4347,4348,4349,4350,4351,4352,4364,4365,4366,4367,4368,4369,4370,4371,4372,4373,4374,4375,4376,4377,4378,4379,4380,4381,4382,4383,4384,4385,4386,4387,4388,4389,4390,4391,4392,4393,4394,4395,4396,4397,4398,4399,4408,4409,4410,4411,4412,4413,4414,4421,4422,4423,4424,4426,4427,4428,4429,4430,4431,4432,4433,4434,4435,4436,4437,4438,4439,4440,4441,4442,4447,4448,4449,4450,4451,4452,4453,4454,4455,4478,4479,4480,4481,4495,4496,4497,4498,4499,4500,4501,4502,4503,4504,4505,4506,4507,4508,4509,4624,4625,4626,4627,4628,4629,4630,4631,4632,4633,4634,4635,4643,4644,4645,4646,4647,4651,4652,4653,4654,4655,4656,4667,4668,4669,4692,4693,4694,4695,4696,4697,4698,4699,4700,4701,4719,4720,4721,4722,4723,4724,4725,4726,4727,4728,4729,4730,4731,4743,4744,4747,4748,4749,4750,4751,4784,4785,4786,4787,4788,4810,4811,4812,4813,4814,4815,4817,4818,4819,4822,4823,4824,4825,4826,4827,4828,4829,4830,4831,4832,4833,4834,4835,4858,4859,4871,4872,4873,4874,4882,4883,4884,4885,4886,4910,4911,4912,4913,4914,4915,4916,4920,4921,4922,4923,4924,4925,4926,4927,4932,4933,4934,4966,4967,4968,4969,4970,4971,4972,4973,4974,4975,4976,4977,4978,4979,4980,4981,4982,4983,4984,4985,4986,4987,4988,4989,4990,4991,4992,4993,5007,5008,5009,5010,5011,5012,5019,5020,5021,5022,5023,5024,5025,5026,5027,5028,5029,5030,5035,5036,5039,5040,5041,5042,5043,5044,5045,5046,5047,5066,5067,5069,5070,5071,5072,5073,5074,5075,5076,5077,5078,5079,5080,5081,5082,5083,5084,5085,5086,5088,5089,5090,5093,5094,5095,5096,5097,5098,5099,5100,5101,5102,5103,5104,5105,5106,5107,5108,5109,5110,5111,5112,5113,5114,5115,5116,5117,5118,5119,5141,5142,5143,5144,5145,5146,5147,5148,5149,5150,5151,5152,5153,5154,5155,5156,5157,5158,5159,5160,5171,5172,5173,5174,5175,5176,5177,5178,5179,5195,5196,5197,5198,5199,5200,5201,5202,5205,5206,5207,5208,5209,5223,5229,5230,5231,5232,5233,5234,5235,5236,5289,5295,5296,5297,5298,5299,5300,5301,5302,5303,5304,5305,5306,5307,5308,5309,5310,5311,5312,5313,5314,5315,5316,5317,5318,5319,5447,5448,5449,5450,5451,5585,5586,5587,5588,5589,5590,5591,5592,5593,5594,5595,5596,5597,5598,5599,5600,5601,5602,5604,5605,5606,5607,5608,5609,5610,5611,5612,5613,5614,5617,5618,5619,5620,5621,5622,5623,5624,5625,5626,5627,5635,5636,5637,5638,5639,5640,5641,5652,5653,5654,5655,5664,5665,5666,5667,5668,5676,5677,5678,5679,5680,5686,5687,5688,5689,5690,5691,5692,5693,5694,5695,5698,5699,5700,5701,5702,5703,5704,5705,5729,5730,5731,5732,5733,5734,5747,5748,5749,5750,5751,5752,5753,5754,5755,5756,5757,5769,5770,5771,5772,5773,5774,5789,5790,5791,5792,5793,5794,5795,5796,5797,5798,5799,5800,5819,5820,5845,5846,5847,5848,5849,5850,5851,5852,5853,5854,5855,5856,6058,6059,6060,6061,6062,6063,6064,6065,6067,6068,6069,6070,6071,6282,6283,6284,6285,6286,6287,6288,6289,6290,6291,6292,6293,6294,6295,6296,6297,6298,6299,6300,6301,6302,6303,6304,6305,6306,6307,6308,6309,6310,6312,6313,6314,6315,6316,6317,6318,6319,6320,6321,6322,6325,6326,6327,6328,6329,6330,6331,6332,6333,6334,6335,6336,6337,6338,6339,6340,6343,6344,6345,6346,6347,6348,6349,6350,6351,6352,6353,6354,6355,6356,6357,6363,6364,6365,6366,6367,6368,6369,6370,6371,6379,6380,6381,6382,6383,6384,6392,6393,6394,6395,6396,6402,6403,6404,6405,6406,6407,6408,6409,6443,6444,6445,6446,6447,6448,6449,6450,6455,6456,6457,6458,6459,6460,6461,6462,6463,6464,6465,6469,6470,6471,6472,6473,6474,6475,6476,6477,6478,6479,6480,6481,6482,6483,6484,6485,6486,6487,6488,6489,6490,6491,6492,6524,6525,6526,6527,6528,6529,6530,6531,6538,6539,6540,6541,6542,6543,6544,6545,6546,6547,6548,6549,6550,6551,6552,6553,6554,6555,6556,6557,6558,6559,6560,6561,6628,6629,6630,6641,6642,6643,6644,6645,6646,6708,6709,6710,6711,6712,6713,6714,6715,6716,6717,6718,6719,6720,6721,6722,6723,6726,6727,6728,6729,6730,6731,6732,6733,6734,6735,6736,6737,6738,6739,6740,6741,6742,6754,6755,6756,6757,6758,6759,6760,6761,6762,6763,6764,6765,6766,6767,6768,6769,6770,6771,6772,6773,6774,6786,6787,6788,6789,6790,6791,6792,6793,6794,6795,6796,6797,6806,6807,6808,6809,6810,6811,6812,6813,6814,6830,6831,6832,6833,6834,6835,6836,6837,6838,6839,6840,6841,6842,6843,6844,6845,6846,6847,7326,7327,7328,7329,7330,7331,7332,7333,7334,7335,7336,7337,7338,7339,7340,7341,7342,7343,7344,7345,7346,7347,7348,7349,7350,7353,7354,7355,7356,7357,7358,7359,7360,7361,7362,7363,7364,7365,7366,7367,7368,7369,7370,7371,7374,7375,7376,7377,7378,7379,7380,7381,7382,7383,7386,7387,7388,7389,7390,7391,7402,7403,7404,7405,7406,7407,7410,7411,7412,7413,7414,7415,7416,7417,7418,7419,7420,7421,7422,7423,7425,7426,7427,7428,7429,7430,7431,7432,7433,7434,7435,7436,7437,7439,7440,7441,7442,7443,7444,7445,7446,7447,7448,7449,7450,7451,7452,7454,7455,7456,7457,7458,7459,7464,7471,7472,7473,7541,7542,7543,7544,7545,7546,7547,7548,7555,7556,7557,7558,7570,7571,7572,7573,7574,7576,7577,7578,7579,7580,7581,7582,7583,7584,7590,7591,7592,7593,7594,7595,7596,7597,7598,7603,7604,7605,7606,7607,7608,7609,7610,7611,7612,7613,7614,7615,7616,7617,7618,7619,7620,7621,7622,7627,7628,7629,7630,7631,7632,7633,7634,7641,7642,7643,7644,7645,7646,7667,7668,7669,7670,7671,7672,7677,7678,7679,7680,7681,7682,7688,7689,7690,7691,7692,7693,7694,7701,7702,7703,7705,7706,7707,7708,7709,7718,7719,7720,7721,7722,7723,7724,7725,7726,7727,7728,7729,7730,7731,7732,7733,7734,7744,7765,7766,7767,7801,7876,7877,7878,7879,7892,7893,7894,7895,7896,7897,7898,7899,7900,7901,7902,7904,7905,7906,7908,7909,7910,7911,7912,7913,7914,7915,7916,7936,7937,7938,7939,7940,7941,7942,7943,7944,7945,7946,7947,7948,7949,7951,7952,7953,7954,7955,7959,7960,7963,7964,7965,7966,7967,7968,7969,7970,7971,7972,7975,7976,7997,7998,7999,8000,8001,8002,8003,8004,8005,8006,8007,8008,8009,8010,8011,8012,8013,8014,8015,8016,8017,8018,8019,8020,8021,8022,8023,8024,8025,8026,8027,8028,8029,8030,8031,8032,8033,8034,8035,8036,8037,8038,8039,8042,8043,8044,8045,8046,8047,8048,8049,8050,8051,8052,8053,8054,8055,8056,8057,8058,8059,8060,8061,8062,8067,8068,8069,8070,8071,8072,8073,8145,8146,8147,8148,8149,8150,8151,8152,8153,8158,8159,8160,8161,8162,8163,8164,8165,8166,8167,8168,8169,8175,8176,8177,8178,8179,8180,8181,8182,8183,8184,8185,8239,8240,8241,8330,8338,8339,8340,8341,8342,8343,8344,8345,8346,8347,8359,8360,8424,8425,8426,8427,8428,8429,8430,8431,8432,8433,8434,8435,8450,8451,8452,8453,8454,8455,8456,8457,8458,8459,8460,8461,8462,8463,8464,8465,8466,8467,8468,8469,8470,8471,8472,8473,8474,8475,8476,8477,8478,8479,8480,8481,8482,8483,8484,8485,8486,8487,8488,8489,8490,8491,8492,8493,8494,8495,8528,8529,8530,8531,8532,8533,8534,8535,8536,8537,8538,8539,8540,8541,8542,8543,8544,8555,8556,8557,8558,8559,8560,8561,8625,8626,8627,8628,8629,8630,8632,8633,8634,8635,8636,8637,8638,8639,8640,8641,8642,8643,8644,8645,8646,8647,8648,8649,8650,8651,8652,8653,8654,8655,8656,8657,8658,8659,8660,8661,8662,8663,8669,8670,8671,8672,8673,8674,8676,8680,8681,8682,8683,8684,8689,8690,8691,8692,8693,8694,8695,8696,8697,8698,8699,8700,8701,8702,8703,8704,8705,8706,8707,8708,8712,8713,8715,8716,8717,8718,8719,8720,8721,8722,8724,8733,8734,8735,8736,8737,8757,8758,8759,8760,8761,8762,8763,8764,8765,8766,8767,8768,8769,8770,8771,8772,8773,8895,8896,8897,8898,8899,8900,8901,8902,8903,8904,8905,8906,8907,8908,8909,8910,8911,8912,8913,8914,8915,8917,8918,8919,8920,8921,8922,8923,8924,8925,8926,8927,8928,8929,8930,8931,8932,8933,8934,8935,8936,8941,8942,8943,8944,8945,8946,8947,8948,8949,8950,8951,8953,8954,8955,8957,9010,9011,9012,9013,9014,9015,9016,9017,9018,9019,9020,9021,9022,9023,9024,9025,9026,9027,9028,9029,9030,9031,9032,9033,9035,9036,9131,9132,9133,9134,9135,9136,9137,9138,9139,9140,9141,9142,9143,9144,9159,9160,9161,9162,9163,9164,9165,9166,9167,9168,9169,9170,9171,9172,9186,9187,9188,9189,9190,9191,9192,9193,9194,9196,9197,9198,9199,9200,9201,9202,9203,9204,9205,9206,9207,9208,9209,9210,9228,9229,9230,9231,9232,9233,9234,9235,9328,9329,9330,9430,9447,9448,9449,9450,9451,9458,9468,9842,9843,9844,9845,9846,9847,9850,9851,9852,9853,9854,9855,9856,9857,9861,9862,9863,9864,9865,9866,9867,9868,9869,9870,9871,9872,9873,9874,9875,9876,9877,9878,9879,9880,9881,9882,9883,9884,9885,9886,9887,9888,9889,9890,9891,9892,9893,9894,9895,9896,9897,9898,9902,9903,9905,9906,9907,9908,9915,9916,9917,9918,9919,9920,9921,9922,9923,9924,9925,9926,9927,9928,9929,9930,9932,9933,9934,9935,9936,9937,9941,9942,9943,9944,9945,9946,9947,9948,9949,9950,9951,9952,9953,9954,9955,9956,9957,9958,9959,9960,9961,9962,9963,9964,9965,9966,9967,9968,9969,9970,9971,9972,9973,9974,9975,9978,9979,9980,9981,9984,9985,9986,9987,9988,9989,9990,9991,9992,9993,9994,9995,9996,9997,9998,9999,10000,10001,10004,10005,10006,10007,10008,10009,10010,10011,10012,10016,10017,10018,10019,10020,10021,10025,10026,10027,10028,10029,10030,10031,10032,10033,10034,10060,10061,10062,10063,10064,10067,10068,10069,10070,10071,10072,10073,10074,10075,10080,10081,10082,10083,10084,10085,10086,10087,10088,10089,10096,10097,10098,10099,10100,10117,10118,10119,10120,10122,10123,10124,10125,10126,10127,10128,10129,10130,10131,10132,10133,10134,10135,10136,10137,10138,10139,10140,10141,10142,10154,10155,10156,10157,10158,10159,10160,10161,10163,10164,10165,10166,10167,10168,10169,10171,10172,10173,10179,10180,10181,10182,10183,10184,10185,10186,10187,10188,10276,10277,10278,10279,10280,10281,10282,10283,10284,10285,10286,10287,10290,10291,10292,10293,10294,10295,10296,10297,10298,10299,10300,10301,10302,10303,10304,10305,10306,10307,10308,10309,10310,10311,10312,10313,10314,10315,10316,10317,10318,10319,10320,10321,10322,10329,10330,10331,10332,10333,10334,10335,10336,10337,10338,10342,10343,10345,10346,10347,10348,10349,10350,10351,10352,10354,10358,10482,10483,10484,10485,10486,10487,10488,10489,10490,10491,10492,10494,10495,10496,10497,10498,10499,10500,10501,10502,10503,10504,10505,10506,10507,10508,10509,10510,10511,10512,10513,10514,10515,10516,10519,10520,10521,10522,10523,10524,10525,10526,10527,10528,10529,10530,10531,10532,10542,10543,10544,10545,10546,10547,10548,10549,10550,10551,10552,10554,10555,10556,10557,10558,10559,10560,10561,10562,10569,10570,10571,10572,10573,10574,10575,10576,10577,10578,10579,10580,10581,10582,10583,10584,10704,10705,10706,10707,10708,10709,10732,10733,10734,10735,10736,10737,10738,10739,10740,10741,10742,10749,10750,10751,10758,10775,10776,10777,10785,10786,10800,10801,10811,10818,10819,10820,10821,10822,10823,10824,10825,10835,10836,10837,10838,10839,10840,10841,10842,10843,10844,10845,10846,10847,10848,10849,10889,10890,10891,10892,10893,10894,10895,10896,10899,10900,10922,10923,10924,10925,10926,10931,10937,10938,10939,10940,10941,10942,10948,10957,10958,10959,10960,10961,10962,10963,10964,10965,10966,10967,10968,10969,10970,10971,10972,10973,10974,10975,10980,10981,10982,10983,10984,10985,10986,10987,10988,10989,10990,10991,10992,10993,10994,10995,10996,10997,11001,11002,11007,11008,11009,11010,11011,11012,11013,11014,11015,11016,11017,11018,11019,11020,11021,11022,11023,11024,11025,11026,11027,11028,11029,11030,11031,11032,11033,11034,11035,11036,11037,11038,11039,11040,11041,11042,11043,11044,11045,11046,11047,11048,11049,11050,11055,11062,11063,11064,11065,11066,11071,11072,11073,11074,11075,11076,11077,11078,11079,11080,11081,11082,11083,11084,11085,11086,11087,11088,11089,11090,11091,11092,11093,11094,11095,11096,11097,11098,11099,11228,11239,11240,11241,11242,11270,11271,11272,11273,11274,11275,11276,11277,11278,11287,11288,11289,11290,11291,11292,11293,11294,11295,11296,11297,11298,11299,11300,11301,11302,11303,11304,11305,11306,11307,11308,11309,11310,11311,11312,11313,11314,11315,11316,11317,11318,11319,11320,11321,11322,11323,11324,11325,11326,11327,11328,11329,11330,11331,11334,11335,11336,11337,11339,11340,11341,11342,11343]
#
#
#csv_file = 'hruwr_missing.dat'
#filein = open(csv_file,'w')
#
##atxt = 'HRU_ID, WR_ID, PRIOR, HRU_PRIOR'
#
#for i in range(len(data)):
#    atxt = str(data[i]).rjust(6) + ''.rjust(3)
#    atxt = atxt + str(9999).rjust(6) + ''.rjust(3)
#    atxt = atxt + str(0).rjust(4) + ''.rjust(3)
#    atxt = atxt + str(1).rjust(4)
#
#    filein.write(atxt + '\n') 
#filein.close()

#%%

csv_file = output_path + '/Final/hruwr.dat'
filein = open(csv_file,'w')

#atxt = 'HRU_ID, WR_ID, PRIOR, HRU_PRIOR'

for i in range(len(hrwwr_dat_arr)):
    atxt = str(hrwwr_dat_arr[i,0]).rjust(6) + ''.rjust(3)
    atxt = atxt + str(hrwwr_dat_arr[i,2]).rjust(6) + ''.rjust(3)
    atxt = atxt + str(hrwwr_dat_arr[i,6]).rjust(4) + ''.rjust(3)
    atxt = atxt + str(hrwwr_dat_arr[i,7]).rjust(4)

    filein.write(atxt + '\n') 
filein.close()


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
#outShapefile = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT\Umatilla_InterACTWEL_QSWATv4\Umatilla_InterACTWEL_QSWATv4\Watershed\Shapes\hru2.shp'    
#
#output_path = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\HRUs_Meghna\Zonal_Stats\WR_Tiffs' 
#raster_file = output_path + '/WR_POU_SW_Clip.tif'
#BoundaryRaster = output_path + '/temp_boundary.tif'
#
#QSWAT_utils.CreateTiff(outShapefile, BoundaryRaster, raster_file,'HRUGIS','','','tiff',gdal.GDT_Float64)
#hru_raster, hru_NoData, hru_obj_temp = QSWAT_utils.Read_Raster(BoundaryRaster)
#hru_raster = np.asarray(np.flip(hru_raster,0),dtype=float)
#
#newRaster_name = output_path + '/UHRUs_WR_Region_Latest_Clip_V1.tif'
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