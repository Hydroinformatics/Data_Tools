# -*- coding: utf-8 -*-
import os, shutil, csv, random, pyodbc, re, sys, json
from osgeo import ogr, gdal
from scipy import stats
import numpy as np

#sys.path.append('C:/Program Files/QGIS 2.18/apps/qgis-ltr/python/plugins')
#sys.path.append('C:/Program Files/QGIS 2.18/apps/qgis-ltr/python')
#sys.path.append('C:/Program Files/QGIS 2.18/apps/Python27/Lib')
#sys.path.append('C:/Program Files/QGIS 2.18/apps/Python27/Lib/site-packages')

#import processing 
#from qgis.analysis import QgsZonalStatistics
#from PyQt4.QtCore import QFileInfo
#from qgis.core import *

#sys.path.append('C:/Program Files/QGIS 2.18/apps/qgis-ltr/python/qgis/PyQt')

#os.chdir('../parsers')
#import QSWAT_preprocess

#%%
def StringToRaster(raster):
    # Check if string is provided

    fileInfo = QFileInfo(raster)
    path = fileInfo.filePath()
    baseName = fileInfo.baseName()

    layer = QgsRasterLayer(path, baseName)
    
    return layer

#%%
def copyshp(basefile,copyfile):
    
    shutil.copy2(basefile, copyfile)
    
    exist = os.path.isfile(basefile[:-4] +'.cpg')
    if exist:
        shutil.copy2(basefile[:-4] +'.cpg', copyfile[:-4] + '.cpg')
    
    exist = os.path.isfile(basefile[:-4] +'.dbf')
    if exist:
        shutil.copy2(basefile[:-4] +'.dbf', copyfile[:-4] + '.dbf')
    
    exist = os.path.isfile(basefile[:-4] +'.prj')
    if exist:
        shutil.copy2(basefile[:-4] +'.prj', copyfile[:-4] + '.prj')
    
    exist = os.path.isfile(basefile[:-4] +'.shx')
    if exist:
        shutil.copy2(basefile[:-4] +'.shx', copyfile[:-4] + '.shx')

        
#%%
def Zonal_Statistic_CDL(files, years):
    
    HRU_CDLdict = None
    HRU_NLCDdict = None
    hrusid_cdl = None
    
    if 'landuse_file' in files:
        print ('Zonal_Stats of InterACTWEL_Landuse')
        zonal_shape = files['output_path'] + '\zonal_shape_HRUids.shp'
        copyshp(files['hru_file'],zonal_shape)
    
        zonal_lyr = QgsVectorLayer(zonal_shape, 'zonal_shape_HRUids', 'ogr')
        landuse_lyr = files['landuse_file'].replace('\\','/')
        zoneStat = QgsZonalStatistics(zonal_lyr,landuse_lyr,"", 1, QgsZonalStatistics.Majority)
        zoneStat.calculateStatistics(None)
    
        hrusid_cdl = []
        lyr = QgsVectorLayer(zonal_shape, '', 'ogr')
        for feat in lyr.getFeatures():
            attrs = feat.attributes()
            if int(attrs[len(attrs)-1]) > 148:
                hrusid_cdl.append(attrs[6]) 
    
    if 'nlcd_file' in files:
        print ('Zonal_Stats of NLCD')
        zonal_shape = files['output_path'] + '\zonal_shape_NLCD.shp'
        copyshp(files['hru_file'], zonal_shape)
        
        zonal_lyr = QgsVectorLayer(zonal_shape, 'zonal_shape_NLCD', 'ogr')
        landuse_lyr = files['nlcd_file'].replace('\\','/')

        zoneStat = QgsZonalStatistics(zonal_lyr, landuse_lyr,"", 1, QgsZonalStatistics.Majority)
        zoneStat.calculateStatistics(None)
    
        HRU_NLCDdict = dict()
        lyr = QgsVectorLayer(zonal_shape, '', 'ogr')
        for feat in lyr.getFeatures():
            attrs = feat.attributes()
            if str(attrs[len(attrs)-1]) == 'NULL':
                attrs[len(attrs)-1] = -999
           
            HRU_NLCDdict[attrs[6]] = int(attrs[len(attrs)-1])
    
    if 'cdl_path' in files:
        # Read and create dictionary of DBF tables. Helps to validate the mapping of raster value to crop type.
        fnames = os.listdir(files['cdl_path'])
        findex = []
        c = 0
        
        for fname in fnames:
            if '.tif.' not in fname and '.tfw' not in fname and '_mfv2.tif' in fname:
                findex.append((c,int(fname[4:8])))
            c += 1    
    
        findex = np.asarray(findex)
    
        HRU_CDLdict = dict()
    
        for year in years: 
            
            print ('Zonal_Stats of ' + fnames[findex[findex[:,1]==year,0][0]])
            cdl_raster = files['cdl_path'] + '\\' + fnames[findex[findex[:,1]==year,0][0]]            
            
            zonal_shape = files['output_path'] + '\zonal_shape_' + str(year) + '.shp'
            
            copyshp(files['hru_file'],zonal_shape)
            zonal_lyr = QgsVectorLayer(zonal_shape, 'zonal_shape_' + str(year), 'ogr')
            landuse_lyr = cdl_raster.replace('\\','/')

            zoneStat = QgsZonalStatistics(zonal_lyr,landuse_lyr,"", 1, QgsZonalStatistics.Majority)
            zoneStat.calculateStatistics(None)
            
            temp_dict = dict()
            lyr = QgsVectorLayer(zonal_shape, '', 'ogr')
            for feat in lyr.getFeatures():
                attrs = feat.attributes()
                #print 'Attr' + str(attrs[len(attrs)-1])
                #print attrs[6], attrs[len(attrs)-1]
                if str(attrs[len(attrs)-1]) == 'NULL':
                    attrs[len(attrs)-1] = -999
                temp_dict[attrs[6]] = int(attrs[len(attrs)-1])
            HRU_CDLdict[year] = temp_dict 
        
    return HRU_CDLdict, HRU_NLCDdict, hrusid_cdl
 
#%%
def WR_per_HRU(hru_file, hru_dict,wr_file):
    lyr = QgsVectorLayer(hru_file, '', 'ogr')
    
    query = "HRUGIS  =  '000010007'"
    #query = '"HRUGIS"  =  '000010007'
    print query
    selection = lyr.getFeatures(QgsFeatureRequest().setFilterExpression(query))
    lyr.setSelectedFeatures([hru.id() for hru in selection])
    
    zonal_shape = files['output_path'] + '\zonal_shape_HRUids.shp'
    _writer = QgsVectorFileWriter.writeAsVectorFormat(lyr, zonal_shape, "utf-8", None, "ESRI Shapefile", onlySelected=True)
    
    
    layer1 = QgsVectorLayer(hru_file, '', 'ogr')
    layer2 = QgsVectorLayer(wr_file, '', 'ogr')
    
    intersections = []
    for a in layer1.getFeatures():
        for b in layer2.getFeatures():
            #print a.geometry().intersects(b.geometry())
            if a.geometry().intersects(b.geometry()):
                intersection = a.geometry().intersection(b.geometry())
                intersections.append(intersection.geometry().area())
    
    #copyshp(files['hru_file'],zonal_shape)
    
    return intersections




#%%
def ReplaceIncompleteSeq(HRU_CDLdict):
    
    temp_hrucdl = HRU_CDLdict
    temp_hrus_nans = []
    for hruid in HRU_CDLdict.keys():
        cnans = 0
        for val in HRU_CDLdict[hruid]:
            if val == 'nan':
                cnans = cnans + 1
        
       # if np.isnan(np.asarray(HRU_CDLdict[hruid], dtype=float)).any():
        if cnans > 0 and cnans != len(HRU_CDLdict[hruid]):
            
            temp_ids = HRU_CDLdict.keys()
            for i in range(len(HRU_CDLdict[hruid])):
                if HRU_CDLdict[hruid][i] != 'nan':
                    temp_ids = FindCropSeq(temp_hrucdl, temp_ids, i, hruid)
            
            temp_seq_array = np.empty((len(temp_ids),len(HRU_CDLdict[hruid])))
            for i in range(len(temp_ids)):
                temp_seq_array[i,:] = HRU_CDLdict[temp_ids[i]]
                
            #print temp_seq_array
            #seq, counts = np.unique(temp_seq_array, return_counts=True, axis=1)
            a = temp_seq_array
            b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
            _, idx,counts = np.unique(b, return_index=True, return_counts=True)
            seq = a[idx]
            
            #print hruid, seq, counts
            maxid = np.argmax(counts)
            temp_hrucdl[hruid] = seq[maxid,:]
            
        elif cnans == len(HRU_CDLdict[hruid]):
            temp_hrus_nans.append(hruid)
    
    return temp_hrucdl, temp_hrus_nans

def FindCropSeq(HRU_CDLdict, hru_ids, colid, target_hru):
    
    temp_ids = []
    for hruid in hru_ids:
        if hruid != target_hru and np.isnan(np.asarray(HRU_CDLdict[hruid], dtype=float)).any() == False and HRU_CDLdict[hruid][colid] == HRU_CDLdict[target_hru][colid]:
            temp_ids.append(hruid)
    if len(temp_ids) == 0:
        temp_ids = hru_ids
    return temp_ids

#%%
def FindHRUMgtFiles(swat_path):
    
    
#    swat_files  = os.listdir(swat_path + '/Scenarios/Default/TxtInOut/')
#    hru_files = [f for f in swat_files if '.hru' in f and f != 'output.hru']
#    
#    HRUFiles = dict()
#    #hrufile = []
#    for hfile in hru_files:
#        cline = 0
#        with open(swat_path + '/Scenarios/Default/TxtInOut/' + hfile) as search:
#            for line in search:
#                if cline == 0:
#                    linesplit = re.split('\s',line)
#                    for sptline in linesplit:
#                        if 'HRU' in sptline and cline == 0:
#                            sptline = re.split(':',sptline)
#                            HRUFiles[hfile.strip('.hru')] = int(sptline[1])
#                            #hrufile.append((int(sptline[1]),int(hfile.strip('.mgt'))))
#                            cline = cline + 1
    
    HRUFiles = dict()
    lyr = QgsVectorLayer(swat_path, '', 'ogr')
    for feat in lyr.getFeatures():
        attrs = feat.attributes()
        HRUFiles[attrs[6]] = int(attrs[6])
                     
    return HRUFiles

#%%
def CDL_NLCDtoSWATdict(db_path):
    
    conn_str = (
        r'DRIVER={Microsoft Access Driver (*.mdb, *.accdb)};'
        r'DBQ=' + db_path + ';'
        )
    cnxn = pyodbc.connect(conn_str)
    crsr = cnxn.cursor()
    
    crsr.execute('select * from crop')
    cropdict = dict()
    for row in crsr.fetchall():
        cropdict[str(row[2])] = row[1]
        
        
    conn_str = (
    r'DRIVER={Microsoft Access Driver (*.mdb, *.accdb)};'
    r'DBQ=' + db_path + ';'
    )
    cnxn = pyodbc.connect(conn_str)
    crsr = cnxn.cursor()

    crsr.execute('select * from nlcd2001_lu')

    nlcd_cropdict = dict()
    for row in crsr.fetchall():
        nlcd_cropdict[row[1]] = dict()
        nlcd_cropdict[row[1]]['Name'] = str(row[2])
        if str(row[2]) in cropdict.keys():
            nlcd_cropdict[row[1]]['value'] = cropdict[str(row[2])]
            
    conn_str = (
    r'DRIVER={Microsoft Access Driver (*.mdb, *.accdb)};'
    r'DBQ=' + db_path + ';'
    )
    cnxn = pyodbc.connect(conn_str)
    crsr = cnxn.cursor()

    crsr.execute('select * from CDL_lu')

    cdl_cropdict = dict()
    for row in crsr.fetchall():
        if str(row[2]) in cropdict.keys():
            cdl_cropdict[row[1]] = cropdict[str(row[2])]

    return nlcd_cropdict, cdl_cropdict

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
    
    return no_crop_ids

#%%
def ReadInputFile(file_path):
    files = dict() 
    
    with open(file_path,'rb') as search:
        for line in search:
            line = line.decode('ascii')
            if len(re.split('\s',line)) > 0:
                if 'output_path' in line:
                    linesplit = re.split('\s',line)
                    files['output_path'] = linesplit[2].replace('\\','/')
                    
                elif 'nlcd_file' in line:
                    linesplit = re.split('\s',line)
                    files['nlcd_file'] = linesplit[2].replace('\\','/')
                    
                elif 'hru_file' in line:
                    linesplit = re.split('\s',line)
                    files['hru_file'] = linesplit[2].replace('\\','/')
                    
                elif 'landuse_file' in line:
                    linesplit = re.split('\s',line)
                    files['landuse_file'] = linesplit[2].replace('\\','/')
                    
                elif 'cdl_path' in line:
                    linesplit = re.split('\s',line)
                    files['cdl_path'] = linesplit[2].replace('\\','/')
                    
                elif 'wr_pou_file' in line:
                    linesplit = re.split('\s',line)
                    files['wr_pou_file'] = linesplit[2].replace('\\','/')
                    
                elif 'wr_pod_file' in line:
                    linesplit = re.split('\s',line)
                    files['wr_pod_file'] = linesplit[2].replace('\\','/')
                    
                elif 'subs_file' in line:
                    linesplit = re.split('\s',line)
                    files['subs_file'] = linesplit[2].replace('\\','/')
     
    if 'hru_file' not in files:
        print('Error: An HRU shapefile must be provided.')
        sys.exit()
    
    if 'wr_pou_file' not in files:
        print('Error: A WR shapefile must be provided.')
        sys.exit() 
        
    if 'output_path' not in files:
        print('Warning: An output path was not provided. The results will be saved at: ' + os.getcwd())
        files['output_path'] = os.getcwd()
            
    search.close()
    return files

#%%
driver = ogr.GetDriverByName("ESRI Shapefile")

files = None
#file_path = input("Please enter the File Path with input data layers: ")
file_path ='C:\Users\sammy\Documents\Research\SWAT\QSWAT_Input_Data\Umatilla\HRUs_Meghna\Umatilla_cdl_hru.txt'
files = ReadInputFile(file_path.replace('\\','/'))

years = range(2012,2017)

#%%
##hru_name_dict = FindHRUMgtFiles(files['hru_file'][0:files['hru_file'].find('Watershed')-1].replace('\\','/'))
#hru_name_dict = FindHRUMgtFiles(files['hru_file'])
#
##intersections = WR_per_HRU(files['hru_file'],hru_name_dict,files['wr_file'])
#
#layer1 = QgsVectorLayer(files['hru_file'], '', 'ogr')
#layer2 = QgsVectorLayer(files['wr_pou_file'], '', 'ogr')
#    
#intersections = dict()
#for a in layer1.getFeatures():
#    attrs_a = a.attributes()
#    intersections[attrs_a[6]] = dict()
#    intersections[attrs_a[6]]['area_sqm'] = a.geometry().area()
#    
#    temp = dict()
#    for b in layer2.getFeatures():
#        attrs_b = b.attributes()
#        
#        if a.geometry().intersects(b.geometry()):
#            
#            intersection = a.geometry().intersection(b.geometry())
#            #intersections.append(intersection.geometry().area())
#            #intersections.append(intersection.area())
#            #print attrs_a[6], attrs_b[2], b.geometry().area(), intersection.area()
#            temp[attrs_b[2]] = [(intersection.area()/b.geometry().area())*100,intersection.area()]
#    
##    if attrs_a[6] == '000010014':
##            break
#    if len(temp) > 0:
#        intersections[str(attrs_a[6])]['wris'] = temp
#
#print('Done')


#dataSource1 = driver.Open(files['hru_file'], 0)
#layer1 = dataSource1.GetLayer()
#
#feat_a_count = layer1.GetFeatureCount()
#feat_a_counter = 1
#per_compl = 0
#
#intersections = dict()
#for feat_a in layer1:
#    poly_a = feat_a.geometry()
#    attrs_a = feat_a.GetField("HRUGIS")
#    
#    intersections[attrs_a] = dict()
#    intersections[attrs_a]['area_sqm'] = poly_a.GetArea()
#    
#    temp = dict()
#    dataSource2 = driver.Open(files['wr_pou_file'], 0)
#    layer2 = dataSource2.GetLayer()
#    for feat_b in layer2:
#        
#        poly_b = feat_b.geometry()
#        attrs_b = feat_b.GetField("snp_id")
#        #print attrs_b
#        if poly_a.Intersection(poly_b) is not None and poly_a.Intersection(poly_b).GetArea() > 0:            
#            intersection = poly_a.Intersection(poly_b).GetArea()
#            temp[attrs_b] = [(intersection/poly_b.GetArea())*100,intersection]
#    
#            
#    if len(temp) > 0:
#        intersections[str(attrs_a)]['wris'] = temp
#    
#    if (float(feat_a_counter)/feat_a_count)*100 >= per_compl:
#        print 'Completed: ' + str(round((float(feat_a_counter)/feat_a_count)*100)) + ' %'
#        per_compl = per_compl + 10
#    feat_a_counter = feat_a_counter + 1
#
#print('Done')


#%%
##csv_file = files['output_path'] + '\Umatilla_HRU_WR_POU_Match.csv'
#csv_file = files['output_path'] + '\Umatilla_HRU_WR_POU_Match_OGR.csv'
#filein = open(csv_file,'w')
#
#atxt = 'HRU_ID, SNAP_ID, OVERLAP_PER, OVERLAP_AREA'
#filein.write(atxt + '\n')
#
#for hru in intersections.keys():
#    if 'wris' in intersections[hru].keys():
#        for wr in intersections[hru]['wris'].keys():
#            atxt = str(hru) + ',' + str(wr) + ',' + str(intersections[hru]['wris'][wr][0]) + ',' + str(intersections[hru]['wris'][wr][1])
#            filein.write(atxt + '\n')
#    else:
#        atxt = str(hru) + ', No WRIS'
#        filein.write(atxt + '\n')
#        
#filein.close()
#
#print ('Writing CSV')
#print ('Finished writing CSV')

#%%

#(base,suffix) = os.path.splitext(files['wr_pod_file'])
#temp_ind = base.rfind('/')
#wr_base_name = base[temp_ind+1:]
#
#(base,suffix) = os.path.splitext(files['hru_file'])
#temp_ind = base.rfind('/')
#hru_base_name = base[temp_ind+1:]
#
#join_shp = files['output_path'] + '/JoinSHP_Test.shp'
#commandtxt = 'ogr2ogr ' + join_shp + ' ' + files['wr_pod_file'] + ' ' + files['hru_file'] + ' -dialect sqlite -sql "SELECT ST_Intersection(a.geometry, b.geometry) AS geometry FROM ' + wr_base_name + ' a, ' + hru_base_name + ' b"'
##ogr2ogr -sql "SELECT n.Name, n.Capacity, n.geometry, b.BoroName from boroughs b, nursinghomes n WHERE ST_INTERSECTS(b.geometry, n.geometry)" -dialect SQLITE output.shp input.vrt
#
#exitflag = os.system(commandtxt)

dataSource = driver.Open(files['wr_pod_file'], 0)
pod_layer = dataSource.GetLayer()

Pod_Sub_dict = dict()
Pod_dict = dict()

feat_a_count = pod_layer.GetFeatureCount()
feat_a_counter = 1
per_compl = 0

hrudataSource = driver.Open(files['hru_file'], 0)
hru_layer = hrudataSource.GetLayer()

hru_extent = hru_layer.GetExtent() #x_min, x_max, y_min, y_max =

hru_layer_dict = dict()
for hru in hru_layer:
    hru_layer_dict[hru] = dict()
    hru_layer_dict[hru]['GEOM'] = hru.GetGeometryRef()
    hru_layer_dict[hru]['SUBID'] = hru.GetField("Subbasin")
    hru_layer_dict[hru]['HRUID'] = hru.GetField("HRUGIS")


pod_dict = dict()
for pod in pod_layer:    
    podid = pod.GetField("pod_use_id")
    points_json = [json.loads(pod.GetGeometryRef().ExportToJson())]
    pod_coords = points_json[0]['coordinates']
    
    if pod_coords[0] >= hru_extent[0] and pod_coords[0] <= hru_extent[1] and pod_coords[1] >= hru_extent[2] and pod_coords[1] <= hru_extent[3]:
        pod_dict[podid] = dict()
        #pod_dict[podid]['GEOM'] = pod.geometry()
        pod_dict[podid]['wr_type'] = pod.GetField("wr_type")
        pod_dict[podid]['use_code'] = pod.GetField("use_code_d")

dataSource = driver.Open(files['wr_pod_file'], 0)
pod_layer = dataSource.GetLayer()

#for podid in pod_dict.keys():
for pod in pod_layer: 
    
    podid = pod.GetField("pod_use_id")
    #print podid
#    #points_json = [json.loads(pod.GetGeometryRef().ExportToJson())]
#    #pod_coords = points_json[0]['coordinates']
#    pod_coords = pod.geometry()
#    
#    Pod_dict[podid] = dict()
#    Pod_dict[podid]['wr_type'] = pod.GetField("wr_type")
#    Pod_dict[podid]['use_code'] = pod.GetField("use_code_d")
#    
#    #subdataSource = driver.Open(files['subs_file'], 0)
#    #subs_layer = subdataSource.GetLayer()
#    #hrudataSource = driver.Open(files['hru_file'], 0)
#    #hru_layer = hrudataSource.GetLayer()
#    
#    #for sub in subs_layer:
    
    #pod_coords = pod_dict[podid]['GEOM']
    if podid in pod_dict.keys():
        
        pod_coords = pod.geometry()
        
        for hru in hru_layer_dict.keys():
            #sub_points = sub.GetGeometryRef()
            #sub_id = sub.GetField("Subbasin")
            #hru_points = hru.GetGeometryRef()
            #subs_id = hru.GetField("Subbasin")
            #hru_id = hru.GetField("HRUGIS")
            
    #        if pod_coords.Within(sub_points):
    #            Pod_Sub_dict[podid] = sub_id
    #            break
            if pod_coords.Within(hru_layer_dict[hru]['GEOM']):
                Pod_Sub_dict[podid] = [hru_layer_dict[hru]['SUBID'],hru_layer_dict[hru]['HRUID']]
                break
        
    if (float(feat_a_counter)/feat_a_count)*100 >= per_compl:
        print 'Completed: ' + str(round((float(feat_a_counter)/feat_a_count)*100)) + ' %'
        per_compl = per_compl + 10
    feat_a_counter = feat_a_counter + 1
            

#%%        
csv_file = files['output_path'] + '\Umatilla_HRU_WR_POD_Match_OGR.csv'
filein = open(csv_file,'w')

atxt = 'POD_ID, SUB_ID, HRU_ID, WR_TYPE, USE_CODE'
filein.write(atxt + '\n')

for pod in pod_dict.keys():
    if pod in Pod_Sub_dict.keys():
        atxt = str(pod) + ',' + str(Pod_Sub_dict[pod][0]) + ',' + str(Pod_Sub_dict[pod][1])  + ',' + str(pod_dict[pod]['wr_type']) + ',' + str(pod_dict[pod]['use_code'])
        filein.write(atxt + '\n')
    else:
        atxt = str(pod) + ',-,-,' + str(pod_dict[pod]['wr_type']) + ',' + str(pod_dict[pod]['use_code'])
        filein.write(atxt + '\n')
        
filein.close()

print ('Writing CSV')
print ('Finished writing CSV')