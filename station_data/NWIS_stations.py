# -*- coding: utf-8 -*-

import os, xmltodict, csv, urllib2, re
import numpy as np


class NWIS_Station_Data(object):
    def __init__(self,base_url='',params=[]):
        
        if len(base_url) > 0:
            self.station_list_url = base_url
            self.station_info_url = 'https://waterservices.usgs.gov/nwis/site/?site='
            self.station_data_url = ['https://nwis.waterdata.usgs.gov/nwis/dv?cb_', '=on&format=rdb&site_no=', '&period=&begin_date=', '&end_date=']
            self.output_path = os.getcwd()
            self.parameters = params
            self.sites_info = dict()
            self.station_data = dict()
            
        else:
            print 'A url to a list of stations is required.'
            quit
        
#%% 
    def URLXML_to_Dict(self,request):
        try:
            temp_file = urllib2.urlopen(request)
            data = temp_file.read()
            temp_file.close()
    
            data = xmltodict.parse(data)
            return data
        
        except urllib2.HTTPError, err:
            print("Bad Request to URL: " + err.reason)
            print("URL: " + request)
            
#%%        
    def ExtractDates(self,linesplit):
        dates = []
        for lsplit in linesplit:
            if len(lsplit) > 6:
                if lsplit[4] == '-' and lsplit[7] == '-':
                    dates.append(lsplit)
        return dates
        
    def LookForPar(self,linesplit):
        linebool = 0
        for lsplit in linesplit:
            for par in self.parameters:
                if lsplit == par:
                    linebool = 1
                    break
        return linebool
    
#%%       
    def StationInfo(self,site_id):
        base_url = self.station_info_url + site_id
        
        if len(self.parameters) > 0:
            base_url = base_url + '&parameterCd='
            for par in self.parameters:
                base_url = base_url + par + ','
            base_url = base_url[:-1]
        base_url = base_url + '&seriesCatalogOutput=true'
        
        data = urllib2.urlopen(base_url)
        linebool = 0
        site_info = dict()
        clcounter = 0
        
        for line in data:
            temp_dict = dict()
            if linebool == 1 and len(line.strip()) > 0:
                linesplit = re.split('\s',line)
                
                if self.LookForPar(linesplit) == 1:
                    temp_dict['line'] = line
                    temp_dict['dates'] = self.ExtractDates(linesplit)
                    site_info[clcounter] = temp_dict
                    clcounter += 1
            if '5s	15s' in line:
                linebool = 1
        
        return site_info, base_url   
    
#%%    
    def GetStationsInfo(self, date_filter = [], csv_file = ''):
        # Find list of stations with specified search criteria (e.g., lat/long box, county) provided by user and list of parameters
        if len(self.parameters) > 0:
            base_url = self.station_list_url + '&parameterCd=' + ",".join(parameters)
        if len(date_filter) > 0:
            base_url = base_url + '&startDT=' + date_filter[0] + '&endDT=' + date_filter[1]
        base_url = base_url + '&seriesCatalogOutput=true'
        
        doc = self.URLXML_to_Dict(base_url)
        
        Sites = dict()  
        for st_type in doc['mapper'].keys():
            for st_id in range(len(doc['mapper'][st_type]['site'])):
                sno = doc['mapper']['sites']['site'][st_id]['@sno']
                
                site_info, site_url = self.StationInfo(sno)
                
                if len(site_info) > 0:
                    start_dates = []
                    end_dates = []
                    for lineid in site_info.keys():
                        if len(site_info[lineid]['dates']) > 0:
                            temp_date = site_info[lineid]['dates'][0].split('-')
                            start_dates.append(int(''.join(temp_date)))
                            temp_date = site_info[lineid]['dates'][1].split('-')
                            end_dates.append(int(''.join(temp_date)))
                        else:
                            start_dates.append(float("inf"))
                            end_dates.append(0)
                    
                    min_date = np.argmin(start_dates, axis=None)
                    start_date = site_info[min_date]['dates'][0]
                    site_info['start_date'] = start_date
                    max_date = np.argmax(end_dates, axis=None)
                    end_date = site_info[max_date]['dates'][1]
                    site_info['end_date'] = end_date
                    
                    Sites[int(sno)] = site_info 
         
        self.sites_info = Sites    
    
        if len(csv_file) > 0 and len(Sites) > 0:
            with open(self.output_path.replace('\\','/') + '/' + csv_file, mode='wb') as outputcsv:
                outputcsv_writer = csv.writer(outputcsv, delimiter = ',', quotechar = '"', quoting=csv.QUOTE_MINIMAL)
                outputcsv_writer.writerow(['STID','STtype','LAT','LNG','AGENCY','START_DATE','END_DATE'])
                
                for st_type in doc['mapper'].keys():
                    for st_id in range(len(doc['mapper'][st_type]['site'])):
                        sno = doc['mapper']['sites']['site'][st_id]['@sno']
                        cat = doc['mapper']['sites']['site'][st_id]['@cat']
                        lat = doc['mapper']['sites']['site'][st_id]['@lat']
                        lng = doc['mapper']['sites']['site'][st_id]['@lng']
                        agc = doc['mapper']['sites']['site'][st_id]['@agc']
                        outputcsv_writer.writerow([sno,cat,lat,lng,agc,Sites[int(sno)]['start_date'],Sites[int(sno)]['end_date']])
                
            outputcsv.close()

#%%
    def ReadURLData(self, base_url):
        data = urllib2.urlopen(base_url)
        if 'No sites/data found using the selection criteria specified' not in data.read(): 
            data = urllib2.urlopen(base_url)
            line_bool = 0
            temp_dict = dict()
#            old_str_day = ''
#            old_linesplit = ''
#            temp_mean_data = []
            start_date_bool = 0
            for line in data:
                if '#' not in line and line_bool == 2:
                    linesplit = "/".join(line.split())
#                    old_linesplit = linesplit
                    linesplit = re.split('/',linesplit)
                    linedate = re.split('-',linesplit[2])
#                    new_str_day = linesplit[2]
                    if len(linesplit) > 4:
                        try:
                            float(linesplit[3])
                        except:
                            linesplit[3] = '-999.0'
#                    
#                    if old_str_day == new_str_day or old_str_day == '':
#                        temp_mean_data.append(float(linesplit[3]))
#                        old_str_day = new_str_day
#                    else:
#                        old_linesplit = re.split('/',old_linesplit)
#                        old_linedate = re.split('-',old_linesplit[2])
                          
                    if start_date_bool == 0:
                        temp_dict['start_date'] = linesplit[2]
                        start_date_bool = 1
                        
                    temp_dict2 = dict()
                    if str(int(linedate[1]))+'/'+str(int(linedate[2])) not in temp_dict.keys():
                        if len(linesplit) > 4:
                            temp_dict2[int(linedate[0])] = float(linesplit[3])
                        else:
                            temp_dict2[int(linedate[0])] = -999.0
                        
                        temp_dict[str(int(linedate[1]))+'/'+str(int(linedate[2]))] = temp_dict2
                    else:
                        if len(linesplit) > 4:
                            temp_dict[str(int(linedate[1]))+'/'+str(int(linedate[2]))][int(linedate[0])] = float(linesplit[3])
                        else:
                            temp_dict[str(int(linedate[1]))+'/'+str(int(linedate[2]))][int(linedate[0])] = -999.0
                        
#                    if str(linedate[1])+'/'+str(linedate[2]) not in temp_dict.keys():
#                        if len(linesplit) > 1:
#                            temp_dict2[int(linedate[0])] = float(np.mean(temp_mean_data))
#                        else:
#                            temp_dict2[int(old_linedate[0])] = -999.0
#                        
#                        temp_dict[str(old_linedate[1])+'/'+str(old_linedate[2])] = temp_dict2
#                    else:
#                        if len(old_linesplit) > 1:
#                            temp_dict[str(old_linedate[1])+'/'+str(old_linedate[2])][int(old_linedate[0])] = float(np.mean(temp_mean_data))
#                        else:
#                            temp_dict[str(old_linedate[1])+'/'+str(old_linedate[2])][int(old_linedate[0])] = -999.0
                        
#                        temp_mean_data = []
#                        temp_mean_data.append(float(linesplit[3]))
#                        old_str_day = new_str_day
#                        old_linesplit = linesplit
                    
                elif '#' not in line:
                    line_bool = line_bool + 1
              
                                
            temp_dict['end_date'] = linesplit[2]
                        
            data.close()
            return temp_dict
        
        data.close()
        
              
#%%        
    def GetStationData(self, stations_info = dict(), csv_file=''):
        
        if len(self.sites_info) > 0 or len(stations_info.keys()) > 0:
            if len(stations_info) == 0:
                stations_ids = self.sites_info
            else:
                for stid in stations_info.keys():
                    if len(stations_info[stid].date_filter) == 0:
                        print 'Error: Date range (e.g., {"start_date":"2007-01-01","end_date": "2018-01-01"}) needed for all stations.'
                        quit()
                stations_ids = stations_info
            
            station_data = dict()
            for parm in self.parameters: 
                for stid in stations_ids:
                    base_url = self.station_data_url[0] + str(parm) + self.station_data_url[1] + str(stid) + self.station_data_url[2] + stations_ids[stid]['start_date'] + self.station_data_url[3] + stations_ids[stid]['end_date']  
                    station_data[str(stid)+'_'+str(parm)] = self.ReadURLData(base_url)
            
            self.station_data = station_data
            
            if len(csv_file) > 0:
                csv_file = self.output_path.replace('\\','/') + '/' + csv_file
                with open(csv_file, mode='wb') as outputcsv:
                    outputcsv_writer = csv.writer(outputcsv, delimiter = ',', quotechar = '"', quoting=csv.QUOTE_MINIMAL)
                    outputcsv_writer.writerow(['OWNER','STID','PARAMID','DAY','MONTH','YEAR','Mean_daily_flow_cfs'])
                    for stid in station_data.keys():
                        if station_data[stid] != None:
                            stidstr = re.split('_',stid)[0]
                            parm_id = re.split('_',stid)[1]
                            start_date = re.split('-',station_data[stid]['start_date'])
                            end_date = re.split('-',station_data[stid]['end_date'])
                            for year in range(int(start_date[0]),int(end_date[0])+1):
                                #for mon in range(int(start_date[1]),13):
                                for mon in range(1,13):
                                    #for day in range(int(start_date[2]),31):
                                    for day in range(1,32):
#                                        if len(str(mon)):
#                                            mon = '0' + str(mon)
#                                        if len(str(day)):
#                                            day = '0' + str(day)
                                        if str(mon) + '/' + str(day) in station_data[stid].keys() and year in station_data[stid][str(mon) + '/' + str(day)].keys():
                                            temp_val = station_data[stid][str(mon) + '/' + str(day)][year]
                                            outputcsv_writer.writerow(['USGS',stidstr,parm_id, day,mon,year,temp_val])
            outputcsv.close

#%%
date_filter = []
parameters = []
csv_file = []

#csv_file = 'nwis_stations.csv'

#csv_file = 'nwis_stations_ALL_streamflow.csv'
#parameters = ['00060']

#csv_file = 'nwis_stations_2007_2018_streamflow.csv'
#parameters = ['00060']
#date_filter = ['2007-01-01','2018-01-01']

#csv_file = 'nwis_stations_ALL_streamTemp.csv'
#parameters = ['00021']

#csv_file = 'nwis_stations_ALL_totalNitrogen.csv'
#parameters = ['00600']

parameters = ['00060']
base_url = 'https://nwis.waterservices.usgs.gov/nwis/site/?format=mapper&countyCd=41021,41059,41049'

nwis_data = NWIS_Station_Data(base_url,parameters)
nwis_data.GetStationsInfo(date_filter)
nwis_data.GetStationData(csv_file='USGS_MeanDailyStreamflow_Data.csv')

