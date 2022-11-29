# -*- coding: utf-8 -*-

import sys, datetime, gdal, os, shutil, zipfile, json, shapely,pyodbc
from statistics import mode
from scipy import stats, ndimage
from osgeo import ogr
from shapely.geometry import Point
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from rasterio.mask import mask

qswat_utils_path = r'C:\Users\sammy\Documents\GitHub\InterACTWEL\src'
sys.path.append(qswat_utils_path)

from qswat import QSWAT_utils
import owrd_utils

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
#root_path = r'Z:\Projects\INFEWS\Modeling\GIS_Data\Data\wr_state'
root_path = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\HRUs_Meghna'
pou_file = 'wr_v_pou_public_proj_REGION.shp'
pod_file = 'wr_v_pod_public_proj_REGIONv2.shp'
subbasins_file = 'Sub_basins_147/subs1.shp'
pou_file_noduplicate = 'wr_v_pou_public_proj_REGION_NoDuplicatesPy_V2_cr'
output_path = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\HRUs_Meghna\Zonal_Stats\WR_Tiffs'

build_wr_hru_dict = 0
build_sub_pod_dict = 0
build_sub_pou_dict = 0
hru_wr_area_tresh = 5

###############################################################################
#%%

if build_sub_pod_dict == 1:
    sub_pod_dict = owrd_utils.GetSWATBasinForPODs(root_path, subbasins_file, pod_file)
    with open(output_path + '/Final/SUB_POD_DICT_CR.json', 'w') as fp:
        json.dump(sub_pod_dict, fp)
else:
    with open(output_path + '/Final/SUB_POD_DICT_CR.json') as json_file:
        sub_pod_dict = json.load(json_file)
        
        
if build_sub_pou_dict == 1:
    sub_pou_dict = owrd_utils.GetSWATBasinForPOUs(root_path, subbasins_file, pou_file)
    with open(output_path + '/Final/SUB_POU_DICT_CR.json', 'w') as fp:
        json.dump(sub_pou_dict, fp)
else:
    with open(output_path + '/Final/SUB_POU_DICT_CR.json') as json_file:
        sub_pou_dict = json.load(json_file)

feature_dict, remove_features = owrd_utils.Eliminate_Duplicate_Records(root_path, pou_file, pou_file_noduplicate)
pou_file = pou_file_noduplicate + '.shp'

wr_pou_dict = owrd_utils.GetWR_POU_Data(pou_file, root_path)
wr_pod_dict = owrd_utils.GetWR_POD_Data(pod_file, root_path, wr_pou_dict.keys(), sub_pod_dict)

#%% Specific to Umatilla Cold Spring Res - Check with Meghna
for snpid in wr_pod_dict.keys():
    for podid in wr_pod_dict[int(snpid)].keys():
        if wr_pod_dict[int(snpid)][int(podid)]['wr_type_2'] == 'CSRes':
            wr_pod_dict[int(snpid)][int(podid)]['subbasin'] = 122
    
#%%
non_empty_acres, empty_acres = owrd_utils.WR_Parser(wr_pou_dict, wr_pod_dict)
non_empty_acres_arr = np.asarray(non_empty_acres)
empty_acres_arr = np.asarray(empty_acres)

#%% Specific to Umatilla Cold Spring Res - Check with Meghna

for irow in range(len(non_empty_acres_arr)):
    if int(non_empty_acres_arr[irow][1]) == 199215 and int(non_empty_acres_arr[irow][2]) == 257459:
        non_empty_acres_arr[irow][12] = 'CSRes'
        non_empty_acres_arr[irow][18] = 122

#%%
instream_wrs = owrd_utils.GetInstream_WR(wr_pod_dict, sub_pod_dict)

instream_wrs_arr = []
for ik in instream_wrs.keys():
    temp_subbasins = []
    for poui in sub_pou_dict[str(ik)].keys():
        if sub_pou_dict[str(ik)][poui] not in temp_subbasins:
            temp_subbasins.append(sub_pou_dict[str(ik)][poui])
    
    if len(temp_subbasins[0]) > 0:
        for subid in range(len(temp_subbasins[0])):
            for podi in instream_wrs[ik].keys():
                #instream_wrs_arr.append((instream_wrs[ik][podi]['subbasin'], round(instream_wrs[ik][podi]['max_rate']/35.31466621266132,2), instream_wrs[ik][podi]['pump_start'], instream_wrs[ik][podi]['pump_end']))
                instream_wrs_arr.append((ik,temp_subbasins[0][subid], round(instream_wrs[ik][podi]['max_rate']/35.31466621266132,2), instream_wrs[ik][podi]['pump_start'], instream_wrs[ik][podi]['pump_end']))
    else:
        for podi in instream_wrs[ik].keys():
            instream_wrs_arr.append((ik,instream_wrs[ik][podi]['subbasin'], round(instream_wrs[ik][podi]['max_rate']/35.31466621266132,2), instream_wrs[ik][podi]['pump_start'], instream_wrs[ik][podi]['pump_end']))

#%%
subbasin =[sb for sb in owrd_utils.records(root_path + '/' + subbasins_file)]
instream_wrs_temp = np.asarray(instream_wrs_arr)

for sid in range(len(subbasin)):
    sub_id = subbasin[sid]['properties']['Subbasin']
    if sub_id not in np.unique(instream_wrs_temp[:,1]):
        instream_wrs_arr.append((-999, int(sub_id), -1, 1, 365))

instream_wrs_sorted = np.asarray(sorted(np.asarray(instream_wrs_arr), key=lambda x: (x[1], x[3])))

sub_instwr_dict = dict()
for ik in range(len(instream_wrs_sorted)):
    if instream_wrs_sorted[ik,1] not in sub_instwr_dict.keys():
        sub_instwr_dict[int(instream_wrs_sorted[ik,1])] = []
    if instream_wrs_sorted[ik,0] not in sub_instwr_dict[instream_wrs_sorted[ik,1]]:
        sub_instwr_dict[int(instream_wrs_sorted[ik,1])].append(instream_wrs_sorted[ik,0])
        
sub_multiple_wrs = []
for sid in sub_instwr_dict.keys():
    if len(sub_instwr_dict[sid]) > 1:
        temp_tser = np.zeros((365,len(sub_instwr_dict[sid])))
        cc = 0
        for sid2 in sub_instwr_dict[sid]:
            temp_data = instream_wrs_sorted[(instream_wrs_sorted[:,1]==sid) & (instream_wrs_sorted[:,0] == sid2)]
            for irow in range(len(temp_data)):
                for iday in range(int(temp_data[irow,3]),int(temp_data[irow,4])+1):
                    temp_tser[iday-1,cc] = temp_data[irow,2]
            
            cc = cc + 1
        
        min_flow = np.max(temp_tser,axis=1)
        replace_rows = []
        orow = 0
        fday = 1
        for irowa in range(1,365):
            if min_flow[irowa] != min_flow[orow]:
                replace_rows.append((min_flow[orow], fday, orow+1))
                fday = irowa + 1
            
            orow = irowa
        replace_rows.append((min_flow[orow], fday, orow+1))
                 
        for irowb in range(len(replace_rows)):  
            sub_multiple_wrs.append((temp_data[irow,0],temp_data[irow,1],replace_rows[irowb][0],replace_rows[irowb][1],replace_rows[irowb][2]))
    
sub_multiple_arr = np.asarray(sub_multiple_wrs)
sub_multiple_arr = np.asarray(sorted(np.asarray(sub_multiple_arr), key=lambda x: (x[1], x[3])))

#%%
csv_file = output_path + '/Final/instrwr_CR.dat'
filein = open(csv_file,'w')

#for istr in range(len(instream_wrs_sorted)):
istr = 0
for sid in sub_instwr_dict.keys():
    if len(sub_instwr_dict[sid]) > 1:
        temp_data = sub_multiple_arr[sub_multiple_arr[:,1]==sid]
        for irow in range(len(temp_data)):
            atxt = str(int(temp_data[irow][1])).rjust(6) + ''.rjust(3)
            atxt = atxt + str(temp_data[irow][2]).rjust(6) + ''.rjust(3)
            atxt = atxt + str(int(temp_data[irow][3])).rjust(4) + ''.rjust(3)
            atxt = atxt + str(int(temp_data[irow][4])).rjust(4)
            filein.write(atxt + '\n') 
    else:
        temp_data = instream_wrs_sorted[instream_wrs_sorted[:,1]==sid]
        for irow in range(len(temp_data)):
            atxt = str(int(temp_data[irow][1])).rjust(6) + ''.rjust(3)
            atxt = atxt + str(temp_data[irow][2]).rjust(6) + ''.rjust(3)
            atxt = atxt + str(int(temp_data[irow][3])).rjust(4) + ''.rjust(3)
            atxt = atxt + str(int(temp_data[irow][4])).rjust(4)
            filein.write(atxt + '\n') 
    
filein.close()

#for istr in range(len(instream_wrs_sorted)):
#    atxt = str(int(instream_wrs_sorted[istr][0])).rjust(6) + ''.rjust(3)
#    atxt = atxt + str(instream_wrs_sorted[istr][1]).rjust(6) + ''.rjust(3)
#    atxt = atxt + str(int(instream_wrs_sorted[istr][2])).rjust(4) + ''.rjust(3)
#    atxt = atxt + str(int(instream_wrs_sorted[istr][3])).rjust(4)
#    filein.write(atxt + '\n') 
#filein.close()
       
#%%
pou_uses = dict()
pou_uses_wrtype = dict()
pod_uses_wrtype = dict()

for i in range(len(non_empty_acres)):
    use_cd_pou = non_empty_acres[i][4] # from POUs
    use_cd_pod = non_empty_acres[i][10] #from PODs
    wr_type_pou = non_empty_acres[i][5] # from POUs
    wr_type_pod = non_empty_acres[i][12] # from PODs
    
    if use_cd_pou not in pou_uses.keys():
        pou_uses[use_cd_pou] = dict()
        pou_uses_wrtype[use_cd_pou] = dict()
        
    if use_cd_pod not in pod_uses_wrtype.keys():
        pod_uses_wrtype[use_cd_pod] = dict()
    
    if wr_type_pod not in pou_uses[use_cd_pou].keys():
        pou_uses[use_cd_pou][wr_type_pod] = []
        pou_uses_wrtype[use_cd_pou][wr_type_pod] = dict()
        
    if wr_type_pod not in pod_uses_wrtype[use_cd_pod].keys():
        pod_uses_wrtype[use_cd_pod][wr_type_pod] = dict()
    
    snp_id = non_empty_acres[i][1]
    pou_uses[use_cd_pou][wr_type_pod].append(snp_id)
    pou_uses_wrtype[use_cd_pou][wr_type_pod][snp_id] = wr_type_pou
    pod_uses_wrtype[use_cd_pod][wr_type_pod][snp_id] = wr_type_pod
    
pod_uses_wrtype_raw = dict()
for snpid in wr_pod_dict.keys():
    for podid in wr_pod_dict[snpid].keys():
        use_cd_pod = wr_pod_dict[snpid][podid]['use_code']
        wr_type_pod = wr_pod_dict[snpid][podid]['wr_type_2']
        
        if use_cd_pod not in pod_uses_wrtype_raw.keys():
            pod_uses_wrtype_raw[use_cd_pod] = dict()
        
        if wr_type_pod not in pod_uses_wrtype_raw[use_cd_pod].keys():
            pod_uses_wrtype_raw[use_cd_pod][wr_type_pod] = dict()
        
        pod_uses_wrtype_raw[use_cd_pod][wr_type_pod][snp_id] = wr_type_pod
    
#%%
port_of_morrow_pou = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\HRUs_Meghna\SWAT_Port_of_Morrow_POU_SW_CR.shp'        
outPath = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\HRUs_Meghna\Zonal_Stats\WR_Shapefiles'
#temp_use = ['IRRIGATION','SUPPLEMENTAL IRRIGATION']
temp_use = ['IRRIGATION']

clayer = 0
total_wr_vol = dict()

driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(root_path + '/' + pou_file, 0)
wr_layer = dataSource.GetLayer()
spatialRef = wr_layer.GetSpatialRef()


dataSource2 = driver.Open(port_of_morrow_pou, 0)
pom_pou_layer = dataSource2.GetLayer()
spatialRef2 = pom_pou_layer.GetSpatialRef()
for featureb in pom_pou_layer:
    pom_pou_points = featureb.GetGeometryRef()

#port_morrow_subs = [1,3,144,145,146,147]
port_morrow_snpid = [193154, 198391, 198571]

#for use_code in pou_uses.keys():
for use_code in temp_use:
    if use_code not in total_wr_vol.keys():
        total_wr_vol[use_code] = dict()
        
    for wr_type in pou_uses[use_code].keys():
    #for wr_type in temp_wr_type:
        if wr_type not in total_wr_vol[use_code].keys():
            total_wr_vol[use_code][wr_type] = dict()
        
        outFileName = 'wr_v_pou_public_proj_REGION_' + str(clayer) + '_' + str(wr_type) + '_cr.shp'
        print outFileName
        
        outShapefile = outPath + '/' + outFileName
        outDriver = ogr.GetDriverByName("ESRI Shapefile")
    
        # Remove output shapefile if it already exists
        if os.path.exists(outShapefile):
            outDriver.DeleteDataSource(outShapefile)
        
        proj_outFileName = '\wr_v_pou_public_proj_REGION_' + str(clayer) + '_' + str(wr_type) + '_cr.prj'
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
            
            if snp_id in pou_uses[use_code][wr_type] and use_cd == use_code and wr_typ == pou_uses_wrtype[use_code][wr_type][snp_id]:    
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
                
                if snp_id in port_morrow_snpid and wr_typ == 'SW' :
                    points = pom_pou_points
                
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
                        
                        xc, yc = owrd_utils.centeroidnp(ptsc)
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
   
featureb = None
dataSource2 = None 

#%%
output_path = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\HRUs_Meghna\Zonal_Stats\WR_Tiffs'
cdl_path = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\Umatilla_Input_data\CDL_Data\CDL_Data_Filter'
cdl_file = 'CDL_2017_clip_20180821151150_mfv2.tif'

if build_wr_hru_dict == 1:
    fnames = os.listdir(cdl_path)  
    no_crop_ids = Get_No_Crop_Ids()
        
    raster_file = cdl_path + '\\' + cdl_file
    region, region_NoData, region_obj = QSWAT_utils.Read_Raster(raster_file)
    region_prj = region_obj.GetProjection()
    
    hru_wr = dict()
    clayer = 0
    temp_total_wr_vol = dict()
    missed_total_wr_vol = dict()
    
    #for use_code in pou_uses.keys():
    for use_code in temp_use:
        if use_code not in hru_wr.keys():
            hru_wr[use_code] = dict()
            temp_total_wr_vol[use_code] = dict()
            missed_total_wr_vol[use_code] = dict()
            
        for wr_type in pou_uses[use_code].keys():
        #for wr_type in temp_wr_type:
            if wr_type not in hru_wr[use_code].keys():
                hru_wr[use_code][wr_type] = dict()
                temp_total_wr_vol[use_code][wr_type] = 0
                missed_total_wr_vol[use_code][wr_type] = 0
                
            new_wrs = np.zeros(region.shape, dtype=np.int64)
            new_wrs_duty = np.zeros(region.shape, dtype=np.int64)
            
            outFileName = '\wr_v_pou_public_proj_REGION_' + str(clayer) + '_' + str(wr_type) + '_cr.shp'        
            boundary_lyr = outPath + outFileName
            
            dataSource = driver.Open(boundary_lyr, 0)
            temp_layer = dataSource.GetLayer()
            
            cfeat = 1
            for feature in temp_layer:
                outFileName = '/temp.shp'
                outShapefile = output_path + outFileName
                outDriver = ogr.GetDriverByName("ESRI Shapefile")
            
                # Remove output shapefile if it already exists
                if os.path.exists(outShapefile):
                    outDriver.DeleteDataSource(outShapefile)
                    
                proj_outFileName = '/temp.prj'
                with open(output_path + proj_outFileName, 'w') as prj_file:
                    #prj_file.write(str(spatialRef.ExportToWkt()))
                    prj_file.write(str(region_prj))
                prj_file.close()
            
            
                # Create the output shapefile
                outDataSource = outDriver.CreateDataSource(outShapefile)
                outLayer = outDataSource.CreateLayer("temp", geom_type = ogr.wkbPolygon)
            
                # Add an ID field
                idField = ogr.FieldDefn("id", ogr.OFTInteger)
                outLayer.CreateField(idField)
                
                featureDefn = outLayer.GetLayerDefn()
                featureo = ogr.Feature(featureDefn)
                    
                featureo.SetField("id", 1)
                points = feature.GetGeometryRef()
                featureo.SetGeometry(points)
                
                feat_id = feature.GetField("id")
                snp_id = feature.GetField("snp_id")
                pou_id = feature.GetField("pou_use_id")
                wr_vol = feature.GetField("wr_vol")
                wr_duty = feature.GetField("duty")
                
                temp_total_wr_vol[use_code][wr_type] = temp_total_wr_vol[use_code][wr_type] + feature.GetField("wr_vol")
                                
                outLayer.CreateFeature(featureo)
                outDataSource = None
                
                BoundaryRaster = output_path + '/temp_boundary.tif'
                QSWAT_utils.CreateTiff(outShapefile, BoundaryRaster, raster_file,'',1,'','tiff')
                
                region_temp, region_NoData_temp, region_obj_temp = QSWAT_utils.Read_Raster(BoundaryRaster)
                
                if len(np.where(region_temp > 0)) == 0:
                    break
                
                # CHECK WHY THIS IS NEEDED
                region_temp = np.asarray(np.flip(region_temp,0),dtype=float)
                
                if region_NoData_temp != None:
                    region_temp[region_temp == region_NoData_temp] = 0.0
                
                result, num_features = ndimage.measurements.label(region_temp)
                labeled_area = ndimage.measurements.sum(region_temp,result,range(1,num_features+1))
                
                for bi in range(len(labeled_area)):
                    if labeled_area[bi] <= hru_wr_area_tresh and bi != np.argmax(labeled_area): # version < 5
                    #if labeled_area[bi] <= hru_wr_area_tresh*10: # version >= 5
                        region_temp[result == bi+1] = 0
            
                if cfeat == 1 and np.sum(region_temp.flatten()) > 0:
                    new_wrs[region_temp > 0.0] = cfeat
                    hru_wr[use_code][wr_type][cfeat] = [(snp_id, pou_id)]
                    cfeat = cfeat + 1

                elif np.sum(region_temp.flatten()) > 0:
                    old_new_wrs = new_wrs.copy()
                    new_wrs[region_temp > 0.0] = cfeat
                    hru_wr[use_code][wr_type][cfeat] = [(snp_id, pou_id)]
                    
                    uhru = np.unique(old_new_wrs[region_temp > 0.0])
                    uhru = uhru[uhru != 0]
                    uhru = uhru[uhru != cfeat]
                    
                    cfeat = cfeat + 1
                    
                    if len(uhru) > 0:
                        for u in uhru:
                            temp = []
                            temp.append((snp_id, pou_id))
                            #if len(hru_wr[use_code][wr_type][u]):
                            for usnp_id, upou_id in hru_wr[use_code][wr_type][u]:
                                temp.append((usnp_id,upou_id))
                                
                            temp_cond = (region_temp > 0.0) & (old_new_wrs == u)
                            temp_cond = np.asarray(temp_cond, dtype=int)
                            result, num_features = ndimage.measurements.label(temp_cond)
                            labeled_area = ndimage.measurements.sum(temp_cond,result,range(1,num_features+1))
                
                            for bi in range(len(labeled_area)):
                                if labeled_area[bi] <= hru_wr_area_tresh:
                                    temp_cond[result == bi] = 0
                        
                            if np.sum(np.asarray(temp_cond.flatten(), dtype=int)) > hru_wr_area_tresh:
                                hru_wr[use_code][wr_type][cfeat] = temp
                                new_wrs[temp_cond == True] = cfeat
                                cfeat = cfeat + 1
                else:
                    missed_total_wr_vol[use_code][wr_type] = missed_total_wr_vol[use_code][wr_type] + feature.GetField("wr_vol")
                        
    #                else:
    #                    new_wrs[region_temp > 0.0] = cfeat
    #                    hru_wr[use_code][wr_type][cfeat] = [snp_id]
    #                    cfeat = cfeat + 1
                        
            outRasterName = 'WR_POU_proj_REGION_' + str(clayer) + '_' + str(wr_type) + '_cr.tif'
            print outRasterName
            newRaster_name = output_path + '/' + outRasterName
            QSWAT_utils.Save_NewRaster(new_wrs, region_obj, region, newRaster_name, -999.0)
            
        clayer = clayer + 1

    with open('WR_IRR_DICT_HRU_Master_CR.json', 'w') as fp:
        json.dump(hru_wr, fp)
        
    with open('WR_IRR_DICT_HRU_Master_CR.json') as json_file:
        hru_wr = json.load(json_file)
 
else:

    with open('WR_IRR_DICT_HRU_Master_CR.json') as json_file:
        hru_wr = json.load(json_file)
        
######################################################################################
#%%    
output_path = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\HRUs_Meghna\Zonal_Stats\WR_Tiffs'

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

GW_WR_Raster = output_path + '/WR_POU_proj_REGION_0_GW_cr.tif'
GW_WR_Raster_clip = output_path + '/WR_POU_GW_clip_cr.tif'
QSWAT_utils.Clipraster(GW_WR_Raster, GW_WR_Raster_clip, landuseFile, gdal.GRA_Mode)
gwwr_raster, gwwr_NoData, gwwr_obj_temp = QSWAT_utils.Read_Raster(GW_WR_Raster_clip)

#hru_raster = np.asarray(np.flip(hru_raster,0),dtype=float)
#plt.matshow(hru_raster)
if gwwr_NoData != None:
    gwwr_raster[gwwr_raster == gwwr_NoData] = 0.0   
gwwr_flat = gwwr_raster.flatten()


SW_WR_Raster = output_path + '/WR_POU_proj_REGION_0_SW_cr.tif'
SW_WR_Raster_clip = output_path + '/WR_POU_SW_Clip_cr.tif'
QSWAT_utils.Clipraster(SW_WR_Raster, SW_WR_Raster_clip, landuseFile, gdal.GRA_Mode)
swwr_raster, swwr_NoData, swwr_obj_temp = QSWAT_utils.Read_Raster(SW_WR_Raster_clip)

if swwr_NoData != None:
    swwr_raster[swwr_raster == swwr_NoData] = 0.0
swwr_flat = swwr_raster.flatten()


CR_WR_Raster = output_path + '/WR_POU_proj_REGION_0_CR_cr.tif'
CR_WR_Raster_clip = output_path + '/WR_POU_CR_Clip_cr.tif'
QSWAT_utils.Clipraster(CR_WR_Raster, CR_WR_Raster_clip, landuseFile, gdal.GRA_Mode)
crwr_raster, crwr_NoData, crwr_obj_temp = QSWAT_utils.Read_Raster(CR_WR_Raster_clip)

if crwr_NoData != None:
    crwr_raster[crwr_raster == crwr_NoData] = 0.0
crwr_flat = crwr_raster.flatten()


RES_WR_Raster = output_path + '/WR_POU_proj_REGION_0_CSRes_cr.tif'
RES_WR_Raster_clip = output_path + '/WR_POU_CSRes_Clip_cr.tif'
QSWAT_utils.Clipraster(RES_WR_Raster, RES_WR_Raster_clip, landuseFile, gdal.GRA_Mode)
reswr_raster, reswr_NoData, reswr_obj_temp = QSWAT_utils.Read_Raster(RES_WR_Raster_clip)

if reswr_NoData != None:
    reswr_raster[reswr_raster == reswr_NoData] = 0.0
reswr_flat = reswr_raster.flatten()

#%%

swwr_ids = np.asarray(np.unique(swwr_flat, axis=0),dtype = int)
gwwr_ids = np.asarray(np.unique(gwwr_flat, axis=0),dtype = int)
crwr_ids = np.asarray(np.unique(crwr_flat, axis=0),dtype = int)
reswr_ids = np.asarray(np.unique(reswr_flat, axis=0),dtype = int)

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
                
                
unique_cr_snp_ids = []              
for objid in crwr_ids:
    if str(objid) in hru_wr['IRRIGATION']['CR'].keys():
        for snp_ids, pouid in  hru_wr['IRRIGATION']['CR'][str(objid)]:
            if [int(snp_ids),int(pouid),3] not in unique_cr_snp_ids:
                unique_cr_snp_ids.append([int(snp_ids),int(pouid),3])
                
               
unique_res_snp_ids = []              
for objid in reswr_ids:
    if str(objid) in hru_wr['IRRIGATION']['CSRes'].keys():
        for snp_ids, pouid in  hru_wr['IRRIGATION']['CSRes'][str(objid)]:
            if [int(snp_ids),int(pouid),4] not in unique_res_snp_ids:
                unique_res_snp_ids.append([int(snp_ids),int(pouid),4])
                              
                
unique_sw_snp_ids = np.asarray(unique_sw_snp_ids)
unique_gw_snp_ids = np.asarray(unique_gw_snp_ids)
unique_cr_snp_ids = np.asarray(unique_cr_snp_ids)
unique_res_snp_ids = np.asarray(unique_res_snp_ids)

unique_snp_ids = np.concatenate((unique_sw_snp_ids, unique_gw_snp_ids, unique_cr_snp_ids, unique_res_snp_ids))

if len(unique_snp_ids) != (len(unique_sw_snp_ids[:,0]) + len(unique_gw_snp_ids[:,0]) + len(unique_cr_snp_ids[:,0]) + len(unique_res_snp_ids[:,0])):
    print 'Warning: You might have GW and SW water rights sharing a snap_ID'

unique_snp_idsr, usnp_ids_counts = np.unique(unique_snp_ids, axis=0, return_counts = 1)

#%%

use_cd = 'IRRIGATION'

swat_to_wrtype = dict()
swat_to_wrtype['SW'] = 1
swat_to_wrtype['GW'] = 3
swat_to_wrtype['ST'] = 2
swat_to_wrtype['CR'] = 5
swat_to_wrtype['CSRes'] = 2

wr_swat_file = []
wr_id_count = 1
wr_snp_ids = np.asarray(non_empty_acres_arr[:,[1,2]], dtype=int)

temp_wrtracker = dict()

for usnpid, pouid, wrtype in unique_snp_ids:
    row_id = (usnpid == wr_snp_ids[:,0]) & (pouid == wr_snp_ids[:,1])
    row_id = np.where(row_id == True)[0]
    
#    if usnpid == 199215:
#        print usnpid, pouid, wrtype, row_id
    
    for ri in row_id:
        if non_empty_acres[ri][4] == use_cd:
            if usnpid not in temp_wrtracker.keys():
                temp_wrtracker[usnpid] = []
            
            if pouid not in temp_wrtracker[usnpid]:
            
                if non_empty_acres_arr[ri][6] == '-': # Only true for Umatilla based on a manual review of the WRs records
                    non_empty_acres_arr[ri][6] = '1985-01-01'
                    
                prior_date = datetime.datetime.strptime(non_empty_acres_arr[ri][6], '%Y-%m-%d').date()
                subbasin_id = non_empty_acres[ri][18]
                if subbasin_id == 'nan':
                    subbasin_id = 0
                    
                wr_swat_file.append([wr_id_count, usnpid, pouid, swat_to_wrtype[non_empty_acres[ri][12]], prior_date, round(non_empty_acres[ri][14]), 
                                     non_empty_acres[ri][16], non_empty_acres[ri][17], non_empty_acres[ri][17] - non_empty_acres[ri][16] + 1, subbasin_id])
                    
                temp_wrtracker[usnpid].append(pouid)
                wr_id_count = wr_id_count + 1
            else:
                print usnpid, pouid, wrtype, row_id

dates_list = [wr_swat_file[ri][4] for ri in range(len(wr_swat_file))]
prior_dates, sort_index = np.unique(dates_list,return_inverse = 1)

for ri in range(len(wr_swat_file)):
    wr_swat_file[ri][4] = sort_index[ri] + 1

#No WR data
wr_swat_file.append([int('9'*len(str(ri))), int('9'*len(str(ri))), int('9'*len(str(ri))), 0, 0, 0, 365, 365, 365, 0])


csv_file = output_path + '/Final/wr_swat_file_CR.csv'
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
csv_file = output_path + '/Final/wrdata_CR.dat'
filein = open(csv_file,'w')

atxt = 'YEAR_ID, WR_ID, WR_SOURCE_ID, WR_VOL_ft-acre, WR_START_PUMPING, WR_END_PUMPING'
filein.write(atxt + '\n')
for yr in range(0,num_year_sim):
    for i in range(len(wr_swat_file)):
        atxt = str(yr+1).rjust(4) + ''.rjust(3)
        atxt = atxt + str(wr_swat_file[i][0]).rjust(5) + ''.rjust(3)
        atxt = atxt + str(wr_swat_file[i][3]).rjust(4)  + ''.rjust(3)
        atxt = atxt + str(int(wr_swat_file[i][5])).rjust(6)  + ''.rjust(3)
        atxt = atxt + str(wr_swat_file[i][6]).rjust(4)  + ''.rjust(3)
        atxt = atxt + str(wr_swat_file[i][7]).rjust(4)
        #atxt = atxt + str(wr_swat_file[i][7]).rjust(4) + ''.rjust(3)
        #atxt = atxt + str(wr_swat_file[i][9]).rjust(4)
    
        filein.write(atxt + '\n') 
filein.close()


#%%

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
    swat_hruid[int(row[12])] = [row[11],int(row[1])]

#%%
# Version 3 - Only WRs
#temp_lu = np.zeros(lu_region.shape)
#temp_lu[lu_region == 0] = 0
#lu_region = temp_lu
#lu_raster_flat = lu_region.flatten()

## Version 5 - WRs + Small fields from CDL
uni_lnduse = np.unique(lu_region.flatten())
uni_lnduse = uni_lnduse[uni_lnduse > 0]

temp_lu = lu_region.copy()
for lnduse in uni_lnduse:
    if len(np.where(lu_region.flatten() == lnduse)[0]) > 700.0:
        temp_lu[lu_region == lnduse] = 0

lu_region = temp_lu
lu_raster_flat = lu_region.flatten()

combo_arr = np.transpose(np.vstack((lu_raster_flat, swwr_flat, gwwr_flat, crwr_flat, reswr_flat)))

unique_rows = np.unique(combo_arr, axis=0)
no_wrs_data = np.where(np.sum(unique_rows[:,[1,2,3,4]],axis=1) == 0)[0]

hru_wr_dict = []
new_hrus = np.zeros(lu_region.shape, dtype=np.float)

hru_counter = 1
per_complt = 0.

wrs_priority = {1:1, 5:2, 2:3, 3:4}

for i in range(len(unique_rows)):
    
    write_bool = 0
    if i not in no_wrs_data:
        temp_wrs = []
        twrc = 0
        if str(int(unique_rows[i,1])) in hru_wr['IRRIGATION']['SW'].keys():
            for snp_ids, pouid in  hru_wr['IRRIGATION']['SW'][str(int(unique_rows[i,1]))]:
                for ii in range(len(wr_swat_file)):
                    if wr_swat_file[ii][1] == snp_ids and wr_swat_file[ii][2] == pouid and wr_swat_file[ii][3] == 1:
                        temp_wrs.append([hru_counter,wr_swat_file[ii][0],wr_swat_file[ii][4],wr_swat_file[ii][6],wr_swat_file[ii][8],wr_swat_file[ii][3],twrc,wr_swat_file[ii][9],wrs_priority[wr_swat_file[ii][3]]])
        
        if str(int(unique_rows[i,2])) in hru_wr['IRRIGATION']['GW'].keys():
            for snp_ids, pouid in  hru_wr['IRRIGATION']['GW'][str(int(unique_rows[i,2]))]:
                for ii in range(len(wr_swat_file)):
                    if wr_swat_file[ii][1] == snp_ids and wr_swat_file[ii][2] == pouid and wr_swat_file[ii][3] == 3:
                        #hru_wr_dict.append([hru_counter,wr_swat_file[ii][0],3])
                        #temp_wrs.append([hru_counter,wr_swat_file[ii][0],wr_swat_file[ii][4],wr_swat_file[ii][3],twrc,0])
                        temp_wrs.append([hru_counter,wr_swat_file[ii][0],wr_swat_file[ii][4],wr_swat_file[ii][6],wr_swat_file[ii][8],wr_swat_file[ii][3],twrc,wr_swat_file[ii][9],wrs_priority[wr_swat_file[ii][3]]])
                        #twrc = twrc + 1
        
        if str(int(unique_rows[i,3])) in hru_wr['IRRIGATION']['CR'].keys():
            for snp_ids, pouid in  hru_wr['IRRIGATION']['CR'][str(int(unique_rows[i,3]))]:
                for ii in range(len(wr_swat_file)):
                    if wr_swat_file[ii][1] == snp_ids and wr_swat_file[ii][2] == pouid and wr_swat_file[ii][3] == 5:
                        #if snp_ids in port_morrow_snpid and wr_swat_file[ii][9]  
                        temp_wrs.append([hru_counter,wr_swat_file[ii][0],wr_swat_file[ii][4],wr_swat_file[ii][6],wr_swat_file[ii][8],wr_swat_file[ii][3],twrc,wr_swat_file[ii][9],wrs_priority[wr_swat_file[ii][3]]])
                        
        
        if str(int(unique_rows[i,4])) in hru_wr['IRRIGATION']['CSRes'].keys():
            for snp_ids, pouid in  hru_wr['IRRIGATION']['CSRes'][str(int(unique_rows[i,4]))]:
                for ii in range(len(wr_swat_file)):
                    if wr_swat_file[ii][1] == snp_ids and wr_swat_file[ii][2] == pouid and wr_swat_file[ii][3] == 2:
                        temp_wrs.append([hru_counter,wr_swat_file[ii][0],wr_swat_file[ii][4],wr_swat_file[ii][6],wr_swat_file[ii][8],wr_swat_file[ii][3],twrc,wr_swat_file[ii][9],wrs_priority[wr_swat_file[ii][3]]])
        
        temp_wrs = np.asarray(sorted(np.asarray(temp_wrs), key=lambda x: (-x[4], x[3], x[2], x[8])))
        
        temp_wrs_order = 1
        for ii in range(len(temp_wrs)):
            temp_wrs[ii,6] = temp_wrs_order
            temp_wrs_order = temp_wrs_order + 1
        
        ids_bool = (lu_region == unique_rows[i,0]) & (swwr_raster == unique_rows[i,1]) & (gwwr_raster == unique_rows[i,2]) & (crwr_raster == unique_rows[i,3]) & (reswr_raster == unique_rows[i,4])
        

        
        if len(temp_wrs) > 1 and np.sum(np.asarray(ids_bool.flatten(), dtype=int)) <= hru_wr_area_tresh*2.0:
            write_bool = 1
      
        if write_bool == 0:
            for r in range(len(temp_wrs)):
                #hru_wr_dict.append([temp_wrs[r,0],temp_wrs[r,1],temp_wrs[r,2],temp_wrs[r,3],temp_wrs[r,4],temp_wrs[r,5]])
                hru_wr_dict.append([temp_wrs[r,0],temp_wrs[r,1],temp_wrs[r,2],temp_wrs[r,3],temp_wrs[r,4],temp_wrs[r,5],temp_wrs[r,6],temp_wrs[r,7]])
                
            # Create new raster of CDL + WRs
            new_hrus[ids_bool == True] = hru_counter
            
            hru_counter = hru_counter + 1
        
    else:
        
        if unique_rows[i,0] > 0:
            no_wr_id = wr_swat_file[len(wr_swat_file)-1][0]
            hru_wr_dict.append([hru_counter,no_wr_id,no_wr_id,365,365,0,1,0])
            ids_bool = (lu_region == unique_rows[i,0]) & (swwr_raster == unique_rows[i,1]) & (gwwr_raster == unique_rows[i,2]) & (crwr_raster == unique_rows[i,3]) & (reswr_raster == unique_rows[i,4])
            new_hrus[ids_bool == True] = hru_counter
            
            hru_counter = hru_counter + 1
            
    if hru_counter == 4474:
        stp=0
        
    if (i/float(len(unique_rows)))*100.0 > per_complt:
        print 'Completed %: ' + str(per_complt)
        per_complt = per_complt + 10.0

new_hrus[new_hrus == 0.0] = gwwr_NoData
newRaster_name = output_path + '/UHRUs_WR_Region_Latest_Clip_cr.tif'
print ('Writing new HRU raster')
QSWAT_utils.Save_NewRaster(new_hrus, gwwr_obj_temp, gwwr_raster, newRaster_name, gwwr_NoData)

#%%
        
outShapefile = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT\Umatilla_InterACTWEL_QSWATv4\Umatilla_InterACTWEL_QSWATv4\Watershed\Shapes\hru2.shp'    

output_path = r'C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\HRUs_Meghna\Zonal_Stats\WR_Tiffs' 
raster_file = output_path + '/WR_POU_SW_Clip_cr.tif'
BoundaryRaster = output_path + '/temp_boundary.tif'

QSWAT_utils.CreateTiff(outShapefile, BoundaryRaster, raster_file,'HRUGIS','','','tiff',gdal.GDT_Float64)
hru_raster, hru_NoData, hru_obj_temp = QSWAT_utils.Read_Raster(BoundaryRaster)
hru_raster = np.asarray(np.flip(hru_raster,0),dtype=float)

wrRaster_name = output_path + '/UHRUs_WR_Region_Latest_Clip_cr.tif'
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
#hru_wr_dict_arr = np.asarray(hru_wr_dict)
#
#hrwwr_dat = []
#for wrid, hruid in unique_rows:
#    temp_ids = np.where(wrid == hru_wr_dict_arr[:,0])[0]
#    for tid in temp_ids:
#        hrwwr_dat.append(np.hstack((swat_hruid[int(hruid)][0],hru_wr_dict_arr[tid,:])))
#        
#miss_uhrus = np.setdiff1d(uhrus,unique_rows[:,1])
#for hruid in miss_uhrus:
#    hrwwr_dat.append(np.hstack((swat_hruid[int(hruid)][0],[0,9999,9999,0,0,0,1,0])))

hru_wr_dict_arr = np.asarray(hru_wr_dict)
wr_swat_file_arr = np.array(wr_swat_file)

hrwwr_dat = []
for wrid, hruid in unique_rows:
    temp_ids = np.where(wrid == hru_wr_dict_arr[:,0])[0]
    for tid in temp_ids:
        temp_subid = np.where(hru_wr_dict_arr[tid,1] == wr_swat_file_arr[:,0])[0][0]
        if hru_wr_dict_arr[tid,1] != 9999:
            insub_bool = 0
            if wr_swat_file_arr[temp_subid,3] != 5 and swat_hruid[int(hruid)][1] != wr_swat_file_arr[temp_subid,9]:
                
                usnpid = wr_swat_file_arr[temp_subid,1]
                pouid = wr_swat_file_arr[temp_subid,2]
                row_id = (usnpid == wr_snp_ids[:,0]) & (pouid == wr_snp_ids[:,1])
                row_id = np.where(row_id == True)[0]
                
                for ri in row_id:
                    if non_empty_acres[ri][4] == use_cd and non_empty_acres[ri][18] == swat_hruid[int(hruid)][1]:
                        insub_bool = 1
        
        if swat_hruid[int(hruid)][0] == 781:
            stp = 0
                        
        if insub_bool == 1:
            hrwwr_dat.append(np.hstack((swat_hruid[int(hruid)][0], hru_wr_dict_arr[tid,:-1], swat_hruid[int(hruid)][1])))
        else:
            if (wr_swat_file_arr[temp_subid,3] == 1 or wr_swat_file_arr[temp_subid,3] == 3) and wr_swat_file_arr[temp_subid,9] == 0:
                #hrwwr_dat.append(np.hstack((swat_hruid[int(hruid)][0],hru_wr_dict_arr[tid,0:5],5,hru_wr_dict_arr[tid,6:8])))
                hrwwr_dat.append(np.hstack((swat_hruid[int(hruid)][0], hru_wr_dict_arr[tid,:-1], swat_hruid[int(hruid)][1])))
            else:
                hrwwr_dat.append(np.hstack((swat_hruid[int(hruid)][0], hru_wr_dict_arr[tid,:])))
#            hrwwr_dat.append(np.hstack((swat_hruid[int(hruid)][0], hru_wr_dict_arr[tid,:])))
        
miss_uhrus = np.setdiff1d(uhrus,unique_rows[:,1])
for hruid in miss_uhrus:
    hrwwr_dat.append(np.hstack((swat_hruid[int(hruid)][0],[0,9999,9999,0,0,0,1,0])))


#%%

hrwwr_dat_old = hrwwr_dat
hrwwr_dat_new = []

for i in hrwwr_dat_old:
    hrwwr_dat_new .append([i[0],i[2],i[6],i[8],i[7]])
    
tupled_lst = set(map(tuple, hrwwr_dat_new))
hrwwr_dat_new = map(list, tupled_lst)


wrs_priority = {1:1, 5:2, 2:3, 3:4}

hrwwr_dat_arr_b = np.asarray(sorted(np.asarray(hrwwr_dat_new), key=lambda x: (x[0], x[4])))
hrwwr_dat_arr_b = np.unique(hrwwr_dat_arr_b[:,0:5], axis=0)


for rowi in range(len(hrwwr_dat_arr_b)):
    if hrwwr_dat_arr_b[rowi,2] > 0:
        hrwwr_dat_arr_b[rowi,4] = wrs_priority[hrwwr_dat_arr_b[rowi,2]]
    
hrwwr_dat_arr_b = np.unique(hrwwr_dat_arr_b[:,0:5], axis=0)
hrwwr_dat_arr_b = np.asarray(sorted(np.asarray(hrwwr_dat_arr_b), key=lambda x: (x[0], x[4])))


unique_rows = np.unique(hrwwr_dat_arr_b[:,0:2], axis=0)

hrwwr_dat_old = hrwwr_dat_arr_b
hrwwr_dat_new = []


for hruid, wrid in unique_rows:
    #rowsid = np.where((hrwwr_dat_arr_b[:,0] == hruid) & (hrwwr_dat_arr_b[:,1] == wrid))[0]
    rowsid = np.where(hrwwr_dat_arr_b[:,0] == hruid)[0]
    counter_c = 0
    for rwid in rowsid:
        #hrwwr_dat_new[rwid,:] = hrwwr_dat_arr_b[rwid,:]
        hrwwr_dat_arr_b[rwid,4] = counter_c + 1
        counter_c = counter_c + 1

hrwwr_dat_new = hrwwr_dat_arr_b
#for hruid, wrid in unique_rows:
#    rowsid = np.where((hrwwr_dat_arr_b[:,0] == hruid) & (hrwwr_dat_arr_b[:,1] == wrid))[0]
#    if len(rowsid) > 1:
#        temp = []
#        temp_arr = []
#        temp_sort = []
#        for ri in rowsid[0]:
#            if hrwwr_dat_old[ri[0]][1] not in temp:
#                temp.append(hrwwr_dat_old[ri[0]][1])
#                temp_arr.append(hrwwr_dat_old[ri[0]])
#                temp_sort.append(wrs_priority(hrwwr_dat_old[ri[0]][2]))
#                
#        temp_arr_b = np.asarray(sorted(np.asarray(temp_arr), key=lambda x: (x[0], x[4])))
#        unique_prior = np.unique(temp_arr_b[:,4], axis=0)
#        
#        if len(unique_prior) != len(temp_arr_b[:,4]):
#            
#               
#    else:
#            
#        hrwwr_dat_new .append(hrwwr_dat_old[rowsid[0]])
#    
#tupled_lst = set(map(tuple, hrwwr_dat_new))
#hrwwr_dat_new = map(list, tupled_lst)


    
#%%

#hrwwr_dat_arr = np.asarray(sorted(np.asarray(hrwwr_dat_new), key=lambda x: (x[0], -x[7])))
##hrwwr_dat_arr = np.asarray(sorted(np.asarray(hrwwr_dat), key=lambda x: (x[3])))


#csv_file = output_path + '/Final/hruwr_CR.dat'
#filein = open(csv_file,'w')
#
##atxt = 'HRU_ID, WR_ID, PRIOR, HRU_PRIOR'
#
#for i in range(len(hrwwr_dat_arr)):
#    atxt = str(hrwwr_dat_arr[i,0]).rjust(6) + ''.rjust(3)
#    atxt = atxt + str(hrwwr_dat_arr[i,2]).rjust(6) + ''.rjust(3)
#    atxt = atxt + str(hrwwr_dat_arr[i,6]).rjust(4) + ''.rjust(3)
#    atxt = atxt + str(hrwwr_dat_arr[i,8]).rjust(4) + ''.rjust(3)
#    atxt = atxt + str(hrwwr_dat_arr[i,7]).rjust(4)
#
#    filein.write(atxt + '\n') 
#filein.close()


hrwwr_dat_arr = np.asarray(sorted(np.asarray(hrwwr_dat_new), key=lambda x: (x[0], x[4])))

csv_file = output_path + '/Final/hruwr_CR.dat'
filein = open(csv_file,'w')

#atxt = 'HRU_ID, WR_ID, PRIOR, HRU_PRIOR'

for i in range(len(hrwwr_dat_arr)):
    atxt = str(hrwwr_dat_arr[i,0]).rjust(6) + ''.rjust(3)
    atxt = atxt + str(hrwwr_dat_arr[i,1]).rjust(6) + ''.rjust(3)
    atxt = atxt + str(hrwwr_dat_arr[i,2]).rjust(4) + ''.rjust(3)
    atxt = atxt + str(hrwwr_dat_arr[i,3]).rjust(4) + ''.rjust(3)
    atxt = atxt + str(hrwwr_dat_arr[i,4]).rjust(4)

    filein.write(atxt + '\n') 
filein.close()


#%%
#temp = []
#for i in range(len(hrwwr_dat_arr)):
#    if hrwwr_dat_arr[i,6] == 1 or hrwwr_dat_arr[i,6] == 3:
#        if hrwwr_dat_arr[i,8] == 0:
#            temp.append(hrwwr_dat_arr[i,:])
#
##temp = np.asarray(temp)
##tempu = np.unique(temp[:,2])
##temp_snp = []
##for j in range(len(tempu)):
##    temp_snp.append([tempu[j]-1, wr_swat_file_arr[tempu[j]-1,0], wr_swat_file_arr[tempu[j]-1,1]])
##
##temp_snp = np.asarray(temp_snp)