# -*- coding: utf-8 -*-

import sys, datetime, gdal, os, shutil, zipfile, json, shapely, urllib2
from statistics import mode

from osgeo import ogr
from shapely.geometry import Point

import matplotlib.pyplot as plt
import numpy as np

import rasterio
from rasterio.mask import mask
from tools.utils import data_mgt_utils as dmgt
from tools.QGIS_Tools import QGIS_utils as qgis

#qswat_utils_path = r'C:\Users\sammy\Documents\GitHub\InterACTWEL\src'
#sys.path.append(qswat_utils_path)
#from qswat import QSWAT_utils

#%%


class Parser:
    
    def __init__(self):
        
        self.FGDB_url = 'https://arcgis.wrd.state.or.us/data/wr_state.zip'
        self.shp_url = 'https://arcgis.wrd.state.or.us/data/wr_state_shp.zip'
        self.shp_url_ref_map = 'http://apps.wrd.state.or.us/apps/gis/gis_map_library/gis_view_image.aspx?gis_library_image_id=326'
        self.info_url = 'https://www.oregon.gov/OWRD/access_Data/Pages/Data.aspx'
        
        self.data_mgt_obj = dmgt.DataManager()
        
        self.zip_dir = None

#%%    
    def Download(self, destdir, data_type = 'shp', zip_dir = 'OWRD_Zip_Data', del_destdir = 0):
         
        self.zip_dir = destdir.replace('\\','/') + '/' + zip_dir
        
        if os.path.isdir(self.zip_dir) and del_destdir == 1:
            shutil.rmtree(self.zip_dir)
        
        if not os.path.isdir(self.zip_dir):
            os.makedirs(self.zip_dir)
        
        if data_type == 'shp':
            filedata = urllib2.urlopen(self.shp_url)
            datatowrite = filedata.read()
            with open(self.zip_dir + '/wr_state_shp.zip', 'wb') as f:
                f.write(datatowrite)
            
            filedata = urllib2.urlopen(self.shp_url_ref_map)
            datatowrite = filedata.read()
            with open(self.zip_dir + '/wr_map_image.pdf', 'wb') as f:
                f.write(datatowrite)
                
        elif data_type.upper() == 'FGDB':
            filedata = urllib2.urlopen(self.FGDB_url)
            datatowrite = filedata.read()
            with open(self.zip_dir + '/wr_state.zip', 'wb') as f:
                f.write(datatowrite)
                
        f.close()
    
    def Unzip_Data(self , zip_dir = None, uzip_path = None):
        
        if zip_dir and uzip_path:
            self.zip_dir = zip_dir
            self.wrdata_dir = uzip_path
            self.data_mgt_obj.Unzip_Data(uzip_path, zip_dir)
            
        elif not zip_dir and uzip_path:
            self.wrdata_dir = uzip_path
            if not self.zip_dir:
                self.zip_dir = uzip_path
            self.data_mgt_obj.Unzip_Data(uzip_path, self.zip_dir)
            
        elif zip_dir and not uzip_path:
            self.zip_dir = zip_dir
            uzip_path = zip_dir
            self.wrdata_dir = zip_dir
            self.data_mgt_obj.Unzip_Data(uzip_path, zip_dir)
            
        else:
            self.wrdata_dir = self.zip_dir
            uzip_path = self.zip_dir
            self.data_mgt_obj.Unzip_Data(uzip_path, self.zip_dir)

#%%
    def Clip(self, basin_acry, bnd_shp, data_path = None):
        
        if data_path:
            self.wrdata_dir = data_path.replace('\\','/')
        
        for dtype in ['_pou','_pod']:        
            for f in os.listdir(self.wrdata_dir):
                for basin in basin_acry:
                    if basin in f and dtype in f and f[-3:] == 'shp' and '_proj' not in f and '_clip' not in f:
                        if dtype == '_pou':
                            input_layer_file = self.wrdata_dir + '/' + f
                            pou_layer_file = qgis.Clip_Vector(input_layer_file, bnd_shp)
                            self.wr_pou_file = pou_layer_file
                            
                        else:
                            input_layer_file = self.wrdata_dir + '/' + f
                            driver = ogr.GetDriverByName("ESRI Shapefile")
                            
                            dataSource = driver.Open(input_layer_file, 0)
                            input_layer = dataSource.GetLayer()
                            input_spatialRef = input_layer.GetSpatialRef()
    
                            PouDataSource = driver.Open(pou_layer_file, 0)
                            pou_layer = PouDataSource.GetLayer()
                            pou_spatialRef = pou_layer.GetSpatialRef()
                            
                            if input_spatialRef != pou_spatialRef:
                                outShapefile = os.path.dirname(input_layer_file) + '/' + os.path.basename(input_layer_file)[:-4] + '_proj.shp'
                                qgis.Reproject_Vector(input_layer_file, outShapefile, pou_spatialRef)
        
                                input_layer_file = outShapefile
                                dataSource = driver.Open(input_layer_file, 0)
                                input_layer = dataSource.GetLayer()
                                input_spatialRef = input_layer.GetSpatialRef()
                            
                            outShapefile = os.path.dirname(input_layer_file) + '/' + os.path.basename(input_layer_file)[:-4] + '_clip.shp'
   
                            # Remove output shapefile if it already exists
                            if os.path.exists(outShapefile):
                                driver.DeleteDataSource(outShapefile)
                    
                            inLayerDefn = input_layer.GetLayerDefn()
                            outDs = driver.CreateDataSource(outShapefile)
                            outLayer = outDs.CreateLayer(input_layer.GetName(), geom_type=inLayerDefn.GetGeomType())
                            print 'Clipping ' + str(outShapefile)
                            
                            idField = ogr.FieldDefn("OBJECTID", ogr.OFTInteger)
                            outLayer.CreateField(idField)
                            for i in range(0, inLayerDefn.GetFieldCount()):
                                fieldDefn = inLayerDefn.GetFieldDefn(i)
                                outLayer.CreateField(fieldDefn)
                            
                            # get the output layer's feature definition
                            outLayerDefn = outLayer.GetLayerDefn()
                            countf = 1
                            for pouFeature in pou_layer:
                                snp_id = pouFeature.GetField("snp_id")
                                #Execute SQL
                                SQL = "SELECT * FROM " + os.path.basename(input_layer_file)[:-4] + " WHERE snp_id='%s'" %(int(snp_id))
                                Filter = dataSource.ExecuteSQL(SQL)
                                
                                if Filter:
                                    for inFeat in Filter:
                                        geom = inFeat.GetGeometryRef()
                                        outFeature = ogr.Feature(outLayerDefn)
                                        # set the geometry and attribute
                                        outFeature.SetGeometry(geom)
                                        
                                        for i in range(1, outLayerDefn.GetFieldCount()):
                                            tfield_name = outLayerDefn.GetFieldDefn(i).GetNameRef()
                                            tfield_val = inFeat.GetField(tfield_name)
                                            outFeature.SetField(tfield_name, tfield_val)
                                        outFeature.SetField("OBJECTID",int(countf))
                                        
                                        # add the feature to the shapefile
                                        outLayer.CreateFeature(outFeature)
                                        # dereference the features and get the next input feature
                                        outFeature = None
                                        countf = countf + 1
                                
                            # Create projection file
                            proj_outFileName = os.path.dirname(outShapefile) + '/' + os.path.basename(outShapefile)[:-4] + '.prj'
                            with open(proj_outFileName, 'w') as prj_file:
                                prj_file.write(str(input_spatialRef.ExportToWkt()))
                            prj_file.close()
                            
                            # Save and close the shapefiles
                            dataSource = None
                            PouDataSource = None
                            
                            self.wr_pod_file = outShapefile

#%%
    def EliminateDuplicatePOUs(self, output_file, root_path=None, input_file_name=None):    
        
        if not input_file_name:
            input_file_name = self.wr_pou_file
        
        if not root_path:
            root_path = os.path.dirname(input_file_name)
            
        driver = ogr.GetDriverByName("ESRI Shapefile")
        dataSource = driver.Open(input_file_name, 0)
        wr_layer = dataSource.GetLayer()
        wr_layerDefn = wr_layer.GetLayerDefn()
        spatialRef = wr_layer.GetSpatialRef()
        
        outShapefile = root_path + '/' + output_file + '.shp'
        outDriver = ogr.GetDriverByName("ESRI Shapefile")
    
        # Remove output shapefile if it already exists
        if os.path.exists(outShapefile):
            outDriver.DeleteDataSource(outShapefile)
        
        proj_outFileName = root_path + '/' + output_file + '.prj'
        
        with open(proj_outFileName, 'w') as prj_file:
            prj_file.write(str(spatialRef.ExportToWkt()))
        prj_file.close()
    
        # Create the output shapefile
        outDataSource = outDriver.CreateDataSource(outShapefile)
        outLayer = outDataSource.CreateLayer(wr_layer.GetName(), geom_type = wr_layerDefn.GetGeomType())
    
        field_names = []
        for n in range(wr_layerDefn.GetFieldCount()):
            fdefn = wr_layerDefn.GetFieldDefn(n)
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
        self.wr_pou_file = outShapefile
        
        return feature_dict, remove_features
    
#%%
    def GetWR_POU_Data(self, file_name = None, root_path = None):
        
        if not file_name:
            file_name = self.wr_pou_file
        
        if not root_path:
            root_path = os.path.dirname(file_name)
        
        print file_name
        driver = ogr.GetDriverByName("ESRI Shapefile")
        dataSource = driver.Open(file_name, 0)
        wr_layer = dataSource.GetLayer()
        
        wr_pou_dict = dict()
        for feature in wr_layer:
            snap_id = int(feature.GetField("snp_id"))
            if snap_id not in wr_pou_dict.keys():
                wr_pou_dict[snap_id]= dict()
            
            temp_len = int(feature.GetField("pou_use_id"))
            wr_pou_dict[snap_id][temp_len] = dict()
            
            wr_pou_dict[snap_id][temp_len]['Source'] = feature.GetField("wr_type").strip()
            
            #linesplit = feature.GetField("use_code_d").strip()
            linesplit = feature.GetField("use_desc").strip()
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
            
            if feature.GetField("priority"):
                celltext = feature.GetField("priority")
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
            
        self.wr_pou_dict = wr_pou_dict
        return wr_pou_dict

#%%
    def GetWR_POD_Data(self, file_name=None, root_path=None, pou_keys=None):
        
        if not file_name:
            file_name = self.wr_pod_file
        
        if not root_path:
            root_path = os.path.dirname(file_name)
        if not pou_keys:
            pou_keys = self.wr_pou_dict.keys()
        
        driver = ogr.GetDriverByName("ESRI Shapefile")
        dataSource = driver.Open(file_name, 0)
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
                
                use_code = feature.GetField("use_code").strip()
                wr_pod_dict[int(snp_id)][int(pod_id)]['Unmod_Use'] = use_code
                
                if use_code.strip() == 'PRIMARY AND SUPPLEMENTAL IRRIGATION':
                    use_code = 'IRRIGATION'
                    
                if use_code.strip() == 'IRRIGATION AND DOMESTIC':
                    use_code = 'IRRIGATION'
                
                if use_code.strip() == 'DOMESTIC EXPANDED':
                    use_code = 'DOMESTIC'
                
                wr_pod_dict[int(snp_id)][int(pod_id)]['use_code'] = use_code
                
                if feature.GetField("priority"):
                    celltext = feature.GetField("priority")
                    
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
                if float(feature.GetField("acre_feet")) != 0:
                    wr_pod_dict[int(snp_id)][int(pod_id)]['acre_feet'] = float(feature.GetField("acre_feet"))
                #else:
                #    wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['acre_feet'] = float('nan')
                    
                wr_pod_dict[int(snp_id)][int(pod_id)]['source'] = feature.GetField("source_typ").strip()
                if feature.GetField("trib_to"):
                    wr_pod_dict[int(snp_id)][int(pod_id)]['tributary'] = feature.GetField("trib_to").strip()
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
         
        self.wr_pod_dict = wr_pod_dict
        return wr_pod_dict

#%%
    def POU_To_POD_Matcher(self,wr_pou_dict=None, wr_pod_dict=None):
        
        if not wr_pou_dict:
            wr_pou_dict = self.wr_pou_dict
        if not wr_pod_dict:
            wr_pod_dict = self.wr_pod_dict
        
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