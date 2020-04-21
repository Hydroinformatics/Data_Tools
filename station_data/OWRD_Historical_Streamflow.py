# -*- coding: utf-8 -*-
import csv, urllib2, re, os
import numpy as np
from selenium import webdriver

# Gather KML map with all stations
#import xmltodict
## ref website: https://apps.wrd.state.or.us/apps/sw/hydro_report/

#file_kml = 'Z:\Projects\INFEWS\Modeling\FEW_Data\USGS_Gages\owrd_gaging_station.kml'
#with open(file_kml.replace('\\','/')) as fd:
#    doc = xmltodict.parse(fd.read())
#    
#basins = []
#for b in range(len(doc['kml']['Document']['Folder'])):
#    basins.append(doc['kml']['Document']['Folder'][b]['name'])
   

class OWRD_Station_Data(object):
    def __init__(self):
        
        self.station_data_url = ['https://apps.wrd.state.or.us/apps/sw/hydro_report/data_Results.aspx?station_nbr=', '&start_date=','&end_date=','&record_type=']
        self.station_info_url = 'https://apps.wrd.state.or.us/apps/sw/hydro_report/complete_gage_information.aspx?station_nbr='
        self.record_types = {'Mean_daily_flows': 'mdfDaily_time_series'}
        self.output_path = os.getcwd()
        self.chrome_exe_path = ''
        self.data_dir = 'OWRD_Streamflow_Data'
        
        self.OWRD_stations = None
        self. MeanDailyFlowData = None
        
        self.CreateDataFolder()
        
    def CreateDataFolder(self):
        if not os.path.isdir(self.output_path.replace('\\','/') + '/' + self.data_dir):
            os.mkdir(self.output_path.replace('\\','/') + '/' + self.data_dir)

#%%
    def ParseStreamFlowData(self, file_path):
        line_bool = 0
        temp_dict = dict()
        with open(file_path,'rb') as search:
            for line in search:
                if 'Record Date' in line:
                    line = search.next()
                    #line = search.next()
                    line_bool = 1
                    
                if line_bool == 1 and '-----------' not in line:
                    #print line
                    linesplit = "-".join(line.split())
                    linesplit = re.split('-',linesplit)
                    linedate = re.split('/',linesplit[0])
                    temp_dict2 = dict()
                    
                    if str(int(linedate[0]))+'/'+str(int(linedate[1])) not in temp_dict.keys():
                        if len(linesplit) > 1:
                            temp_dict2[int(linedate[2])] = float(linesplit[1])
                        else:
                            temp_dict2[int(linedate[2])] = -999.0
                        
                        temp_dict[str(int(linedate[0]))+'/'+str(int(linedate[1]))] = temp_dict2
                    else:
                        if len(linesplit) > 1:
                            temp_dict[str(int(linedate[0]))+'/'+str(int(linedate[1]))][int(linedate[2])] = float(linesplit[1])
                        else:
                            temp_dict[str(int(linedate[0]))+'/'+str(int(linedate[1]))][int(linedate[2])] = -999.0
        search.close()
        
        return temp_dict
        
#%% If mean daily measurments are selected, a report needs to be downloaded. This fuction uses Chromedrive to
### automatically submit the post command responsible for generating the data report.
    
    def GetMeanDailyFlowData(self, csv_file =''):
        
        #base_url = 'https://apps.wrd.state.or.us/apps/sw/hydro_report/data_Results.aspx?station_nbr=STAID&start_date=01/01/1960&end_date=01/01/2019&record_type=mdfDaily_time_series'
        #desired_caps = {'prefs': {'download': {'default_directory': 'C:\Users\sammy\Documents\GitHub\Data_Tools'}}}
        
        record_type = self.record_types['Mean_daily_flows']
        MeanDailyFlowData = dict()
        
        for stid in self.OWRD_stations.keys():
            print stid
            if self.OWRD_stations[stid]['report']['mean_daily_flow'] is not None:
                start_date = self.OWRD_stations[stid]['report']['mean_daily_flow'][0]
                end_date = self.OWRD_stations[stid]['report']['mean_daily_flow'][1]
                
                serviceurl = self.station_data_url[0] + str(stid) + self.station_data_url[1] + start_date + self.station_data_url[2] + end_date + self.station_data_url[3] + record_type
            
                chromeOptions = webdriver.ChromeOptions()
                prefs = {"download.default_directory" : self.output_path.replace('\\','/') + '/' + self.data_dir}
                chromeOptions.add_experimental_option("prefs",prefs)
            
                driver = webdriver.Chrome(executable_path = self.chrome_exe_path, chrome_options=chromeOptions)
                driver.get(serviceurl)
             
                # click radio button
                python_button = driver.find_elements_by_xpath("//input[@name='ctl00$PageData$btn_download_report']")[0]
                python_button.click()
                
                file_bool = 0
                while file_bool == 0:
                    if os.path.isfile(self.output_path.replace('\\','/') + '/' + self.data_dir + '/station_' + str(stid) + '_mdfDaily_time_series.txt'):
                        driver.close()
                        file_bool = 1
                file_path = self.output_path.replace('\\','/') + '/' + self.data_dir + '/station_' + str(stid) + '_mdfDaily_time_series.txt'
                MeanDailyFlowData[stid] = self.ParseStreamFlowData(file_path)
        
        self.MeanDailyFlowData = MeanDailyFlowData
        
        if len(csv_file) > 0:
            csv_file = self.output_path.replace('\\','/') + '/' + self.data_dir + '/' + csv_file
            with open(csv_file, mode='wb') as outputcsv:
                outputcsv_writer = csv.writer(outputcsv, delimiter = ',', quotechar = '"', quoting=csv.QUOTE_MINIMAL)
                outputcsv_writer.writerow(['OWNER','STID','PARAMID','DAY','MONTH','YEAR','Mean_daily_flow_cfs'])
                
                for stid in MeanDailyFlowData.keys():
                    start_date = re.split('/',self.OWRD_stations[stid]['report']['mean_daily_flow'][0])
                    end_date = re.split('/',self.OWRD_stations[stid]['report']['mean_daily_flow'][1])
                    for year in range(int(start_date[2]),int(end_date[2])+1):
                        #for mon in range(int(start_date[0]),13):
                        for mon in range(1,13):
                            #for day in range(int(start_date[1]),31):
                            for day in range(1,32):
                                if str(mon) + '/' + str(day) in MeanDailyFlowData[stid].keys() and year in MeanDailyFlowData[stid][str(mon) + '/' + str(day)].keys():
                                    temp_val = MeanDailyFlowData[stid][str(mon) + '/' + str(day)][year]
                                    outputcsv_writer.writerow([self.OWRD_stations[stid]['owner'],stid,'MeanDailyFlow',day,mon,year,temp_val])
            outputcsv.close
            
#%% Function gathers the gage stations' information which includes: Owner, water source, data types, basin, reporting periods of each data type
##  Input: Array of stations names (as strings) or CSV with stations names in first columm (firt row is assumed to be a header)
    def GetStationInfo(self, stations_ids, csv_file =''):
        
        if not isinstance(stations_ids, list):
            file_stations = stations_ids
            stations_ids = []
            with open(file_stations,'rb') as search:
                search.next()
                for line in search:
                    linesplit = re.split(',',line)
                    stations_ids.append(linesplit[0])
            
        OWRD_stations = dict()     
        for stid in stations_ids:
            temp_dict = dict()
            #base_url = 'https://apps.wrd.state.or.us/apps/sw/hydro_report/complete_gage_information.aspx?station_nbr=' + stid
            base_url = self.station_info_url + stid
            data = urllib2.urlopen(base_url)
            
            temp_report_dict = dict()
            temp_report_dict['mean_daily_flow'] = None
            temp_report_dict['low_flow'] = None
            temp_report_dict['peak_flow'] = None
            
            for line in data:
                
                if 'ctl00_PageData_lbl_owner' in line:
                    line = line.strip()
                    line = line.replace('<span id="ctl00_PageData_lbl_owner">','')
                    linesplit = re.split('<',line)
                    temp_dict['owner'] = linesplit[0]
                    
                elif 'ctl00_PageData_lbl_source' in line:
                    line = line.strip()
                    line = line.replace('<span id="ctl00_PageData_lbl_source">','')
                    linesplit = re.split('<',line)
                    temp_dict['source'] = linesplit[0]
                
                elif 'ctl00_PageData_lbl_collected_data' in line:
                    line = line.strip()
                    line = line.replace('<span id="ctl00_PageData_lbl_collected_data"><b>','')
                    linesplit = re.split('<',line)
                    temp_dict['datatypes'] = linesplit[0]
                    
                elif 'ctl00_PageData_lbl_basin' in line:
                    line = line.strip()
                    line = line.replace('<span id="ctl00_PageData_lbl_basin">','')
                    linesplit = re.split('<',line)
                    temp_dict['basin'] = linesplit[0]
                    
                elif 'ctl00_PageData_lbl_mean_daily_flow_por' in line and 'No Mean Daily Flow Records' not in line:
                    line = line.strip()
                    line = line.replace('<span id="ctl00_PageData_lbl_mean_daily_flow_por">','')
                    linesplit = re.split('<',line)
                    temp_report_dict['mean_daily_flow'] = linesplit[0].split(' - ')
                    #temp_dict['report'] = linesplit[0].split(' - ')
                
                elif 'ctl00_PageData_lbl_low_flow_por' in line and 'No Low Flow Records Available for this Gaging Station' not in line: 
                    line = line.strip()
                    line = line.replace('<span id="ctl00_PageData_lbl_low_flow_por">','')
                    linesplit = re.split('<',line)
                    temp_report_dict['low_flow'] = linesplit[0].split(' - ')
                    #temp_dict['report'] = linesplit[0].split(' - ')
                    
                elif 'ctl00_PageData_lbl_peak_flow_por' in line: 
                    line = line.strip()
                    line = line.replace('<span id="ctl00_PageData_lbl_peak_flow_por">','')
                    linesplit = re.split('<',line)
                    temp_report_dict['peak_flow'] = linesplit[0].split(' - ')
                    #temp_dict['report'] = linesplit[0].split(' - ')
            temp_dict['report'] = temp_report_dict
            OWRD_stations[stid] = temp_dict
        
        if len(csv_file) > 0:
            csv_file = self.output_path.replace('\\','/') + '/' + self.data_dir + '/' + csv_file
            with open(csv_file, mode='wb') as outputcsv:
                outputcsv_writer = csv.writer(outputcsv, delimiter = ',', quotechar = '"', quoting=csv.QUOTE_MINIMAL)
                outputcsv_writer.writerow(['STID','OWNER','BASIN','START_DATE_Mean_Flow','END_DATE_Mean_Flow','START_DATE_Peak_Flow','END_DATE_Peak_Flow','START_DATE_Low_Flow','END_DATE_Low_Flow'])
                
                for st in OWRD_stations.keys():
                    if 'report' in OWRD_stations[st].keys():
                        outputcsv_writer.writerow([st,OWRD_stations[st]['owner'],OWRD_stations[st]['basin'],OWRD_stations[st]['report']['mean_daily_flow'][0],OWRD_stations[st]['report']['mean_daily_flow'][1],
                                                   OWRD_stations[st]['report']['peak_flow'][0],OWRD_stations[st]['report']['peak_flow'][1],OWRD_stations[st]['report']['low_flow'][0],OWRD_stations[st]['report']['low_flow'][1]])
                
            outputcsv.close()
        
        self.OWRD_stations = OWRD_stations
        #return 

#%%
OWRD_sta = OWRD_Station_Data()  
OWRD_sta.chrome_exe_path = os.path.realpath('..') + '\\chromedriver.exe'

station_file = os.path.realpath('..') + '\\example_input_files\\owrd_gaging_station_Region.csv'
OWRD_sta.GetStationInfo(station_file)

csv_file = 'OWRD_gage_station_data.csv'
#station_id = ['14019109','14019110']
#OWRD_sta.GetStationInfo(station_id)
OWRD_sta.GetMeanDailyFlowData(csv_file)

#file_path = 'C:\Users\sammy\Documents\GitHub\Data_Tools\station_14019109_mdfDaily_time_series.txt'
#temp_dict = OWRD_sta.ParseStreamFlowData(file_path.replace('\\','/'))
