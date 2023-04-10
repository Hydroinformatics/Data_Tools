#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 15:56:15 2023

@author: sammy
"""

import pandas as pd
import binascii, os, re
import numpy as np
import matplotlib.pyplot as plt


# def chunks(file, size=1024):
#     while 1:
#         startat=fh.tell()
#         print(startat) #file's object current position from the start
#         fh.seek(size,1) #offset from current postion -->1
#         data=fh.readline()
#         yield startat,fh.tell()-startat #doesnt store whole list in memory
#         if not data:
#             break
        
# fname = r'/Users/sammy/Library/CloudStorage/Box-Box/Research/SWAT/SWAT_JetStream_runs/ITERS_Results_2/output_0.hru'

# if os.path.isfile(fname):
#     try:
#         fh=open(fname,'rb') 
#     except IOError as e: #file --> permission denied
#         print ("I/O error({0}): {1}".format(e.errno, e.strerror))
#     except Exception as e1: #handle other exceptions such as attribute errors
#         print( "Unexpected error: {0}".format(e1))
        
#     for ele in chunks(fh):
#         fh.seek(ele[0])#startat
#         data=fh.read(ele[1])#endat
#         b_data=(fh.read(ele[1]))#endat This is one chunk of ascii data in binary format
#         a_data=((binascii.b2a_qp(b_data)).decode('utf-8')) #Data chunk in 'split' ascii format
#         data_chunk = (a_data.replace('=\n','').strip()) #Splitting characters removed
#         data_list = data_chunk.split('\n')  #List containing lines in chunk
#         #print(data_list,'\n')
#         #time.sleep(1)
        
        
#         for j in range(len(data_list)): #iterate through data_list to get each item 
#             #i += 1
            
#             line_of_data = data_list[j]
#             print(line_of_data)
            
#         #print (data)


varname = ["LULC","HRU","GIS","SUB","MGT","AREAkm2","BIOMt/ha","YLDt/ha","IRRmm","NAUTOkg/ha","PAUTOkg/ha", "ETmm", "SW_INITmm", "SW_ENDmm", "PERCmm", "GW_RCHGmm", "REVAPmm", "W_STRS"]
varcol = [i for i in range(0,len(varname))]

tfile = r'/Users/sammy/Library/CloudStorage/Box-Box/Research/SWAT/SWAT_JetStream_runs/ITERS_Results_2'

files_hru = ['output_org.hru', 'output_0.hru']


data = dict()

for j in range(0,2):
    
    data_array = dict()
    varbool = 0
    rowid = 0
    year = 1998
    print(tfile + '/' + files_hru[j])
    with open(tfile + '/' + files_hru[j]) as search:
        for line in search:
        
            if 'HRU'.lower() in line.lower():    
                varbool = 1
                
            elif varbool == 1:
                
                linesplit = re.split('\s',line)
                if len(linesplit[0]) > 4:
                    linesplit[2] = linesplit[1]
                    linesplit[1] = linesplit[0][4:]
                    linesplit[0] = linesplit[0][0:4]
                
                linesplit = [e for e in linesplit if e != '']
           
                if int(linesplit[1]) == 12:
                    data_array[rowid] = dict()
                    for i in range(0,len(varcol)):
                        if i != 5 and i != 0:
                            data_array[rowid][varname[i]] = float(linesplit[varcol[i]])
                        elif i != 5 and i == 0:
                            data_array[rowid][varname[i]] = linesplit[varcol[i]]
                        else:
                            
                            if int(linesplit[5].split('.')[0]) == 1:
                                year = year + 1
                            
                            data_array[rowid]['YEAR'] = year
                            data_array[rowid]['MON'] = int(linesplit[5].split('.')[0])
                            data_array[rowid][varname[i]] = float('0.'+ linesplit[5].split('.')[1])
                    
                    rowid = rowid + 1
                
    search.close()

    data[j] = pd.DataFrame.from_dict(data_array, orient='index')
    
#%%

for year in [1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006]:
    for i in range(0,2):
        df = data[i]

        #temp = np.asarray(df.loc[df['YEAR'] == 1999, 'IRRmm'])
        temp = np.asarray(df.loc[df['YEAR'] == 1999, 'YLDt/ha'])
        
        if i == 0:
            plt_data = temp
        else:    
            plt_data = np.vstack((plt_data.T,temp))
        
    plt_data =   plt_data.T
    plt.plot(np.cumsum(plt_data[:,0]),np.cumsum(plt_data[:,1]))
    
    
    