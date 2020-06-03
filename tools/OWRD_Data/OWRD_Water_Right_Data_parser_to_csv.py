# -*- coding: utf-8 -*-

import re, datetime
import numpy as np

#%%
def unique(list1): 
  
    # intilize a null list 
    unique_list = [] 
      
    # traverse for all elements 
    for x in list1: 
        # check if exists in unique_list or not 
        if x not in unique_list: 
            unique_list.append(x)
    
    return unique_list


def GetPOD_CODES(file_name, root_path):
    line_counter = 0
    pod_source_dict = dict()   
    
    with open(root_path + '/' + file_name,'rb') as search: 
        for line in search:
            
            linesplit = line.strip().split(',')

            if line_counter == 0:
                line_counter = line_counter + 1
                
            else:
                pod_source_dict[linesplit[0]] = linesplit[2].strip()
                
    search.close()   
    return pod_source_dict


def GetWR_Data(file_name, root_path):
    #with open(root_path + '/' + file_name,'rb') as search: 
#    for line in search:
#        
#        linesplit = line.strip().split(',')
#        #linesplit = line.strip().split('\t')
#        #linesplit = [i for i in linesplit if len(i) > 0]
#        
#        if line_counter == 0:
#            line_counter = line_counter + 1
#            
#        else:
#            cellcount = 0
#            cellcount_url = 0
#            cellcount_source = -999
#            slashid = []
#            snap_id = 0
#            
#            wr_use = line.strip().split('"')
#            
#            for celltext in linesplit:
#                celltext = celltext.strip(' ')
#                
#                if cellcount_url == 1:
#                    snap_id = int(celltext)
#                    wr_pou_dict[snap_id] = dict()
#                    cellcount_url = 0
#                    cellcount_source = 1
#                    if len(wr_use) > 1 and wr_use[1].isupper():
#                        wr_pou_dict[snap_id]['Use'] = wr_use[1]
#                
#                if 'http://apps.wrd.state.or.us' in celltext:
#                    cellcount_url = 1
#                    
#                if celltext == 'SW' or celltext == 'ST' or celltext == 'GW':
#                    wr_pou_dict[snap_id]['Source'] = celltext
#                    
#                slashid = [i for i in range(len(celltext)) if celltext[i] is '/' or celltext[i] is ':']
#                
#                if snap_id == 127909:
#                    stpo = 0
#                    
#                if len(slashid) > 3 and cellcount_source == 1 and cellcount < 25 :
#                    
#                    year = int(celltext[slashid[1]+1:slashid[1]+5])
#                    day = int(celltext[slashid[0]+1:slashid[1]])
#                    mon = int(celltext[0:slashid[0]])
#                    
#                    wr_pou_dict[snap_id]['Priority'] = datetime.datetime(year,mon,day)
#                    wr_pou_dict[snap_id]['Use'] = linesplit[cellcount - 1]
#                        
#                    cellcount_source = 0
#                
#                cellcount = cellcount + 1
#            
#            if 'Use'not in wr_pou_dict[snap_id].keys() and len(wr_use) <= 1:
#                wr_pou_dict[snap_id]['Use'] = linesplit[21]
#            elif 'Use'not in wr_pou_dict[snap_id].keys() and len(wr_use) > 1:
#                wr_pou_dict[snap_id]['Use'] = linesplit[23]
    line_counter = 0
    wr_pou_dict = dict()   
    
    with open(root_path + '/' + file_name,'rb') as search: 
        for line in search:
            
            linesplit = line.strip().split(',')

            if line_counter == 0:
                line_counter = line_counter + 1
                
            else:

                snap_id = int(linesplit[0])
                if snap_id not in wr_pou_dict.keys():
                    wr_pou_dict[snap_id]= dict()
                
                temp_len = int(linesplit[1])
                wr_pou_dict[snap_id][temp_len] = dict()
                
                wr_pou_dict[snap_id][temp_len]['Source'] = linesplit[2].strip()
                
#                if snap_id == 17085:
#                    sttpo = 0
                wr_pou_dict[snap_id][temp_len]['Unmod_Use'] = linesplit[5].strip()
                if 'IRRIGATION LIVESTOCK' in linesplit[5].strip():
                    linesplit[5] = linesplit[5].replace('N L','N, L')
                
                if ' AND ' in linesplit[5].strip() or ', ' in linesplit[5].strip():
                    temp_use = linesplit[5].strip().split(', ')
                    and_split =  temp_use[len(temp_use)-1].split(' AND ')
                    temp_use[len(temp_use)-1] = and_split[0]
                    temp_use.append(and_split[1])
                    wr_pou_dict[snap_id][temp_len]['Use'] = temp_use
                else:
                    wr_pou_dict[snap_id][temp_len]['Use'] = [linesplit[5].strip()]
                
                celltext = linesplit[6]
                        
                slashid = [i for i in range(len(celltext)) if celltext[i] is '/' or celltext[i] is ':']
                        
                if len(slashid) >= 3:
                        
                    year = int(celltext[slashid[1]+1:slashid[1]+5])
                    day = int(celltext[slashid[0]+1:slashid[1]])
                    mon = int(celltext[0:slashid[0]])
                    
                    wr_pou_dict[snap_id][temp_len]['Priority'] = datetime.datetime(year,mon,day)
                  
                if len(linesplit[8]) > 0:
                    wr_pou_dict[snap_id][temp_len]['acre_feet'] = float(linesplit[8])
                else:
                    wr_pou_dict[snap_id][temp_len]['acre_feet'] = float('nan')
                
#                if 'Use'not in wr_pou_dict[snap_id].keys() and len(wr_use) <= 1:
#                    wr_pou_dict[snap_id]['Use'] = linesplit[14]
#                elif 'Use'not in wr_pou_dict[snap_id].keys() and len(wr_use) > 1:
#                    wr_pou_dict[snap_id]['Use'] = linesplit[15]
    search.close()
    return wr_pou_dict

#%%
def GetWR_POD_Data(file_name, root_path, pou_keys):
    line_counter = 0
    wr_pod_dict = dict()
    with open(root_path + '/' + file_name,'rb') as search: 
        for line in search:
            if line_counter > 0:
                linesplit = line.strip().split(',')
                
#                if int(linesplit[2]) == 170304:
#                        stpo = 0
#                print linesplit[2]
                
                wr_use = line.strip().split('"')
                
                temp_list = []
                if len(wr_use) > 1:
                    for i in range(1,len(wr_use)-1):
                        cc = 0
                        for ii in range(len(linesplit)):
                            if linesplit[ii].strip('"') in wr_use[i] and len(linesplit[ii]) > 0 and '"' in linesplit[ii]:
                                if cc == 0:
                                    temp_list.append(wr_use[i])
                                    cc = 1
                                #else:
                                    
                            else:
                                temp_list.append(linesplit[ii])
                                
                    #linesplit = unique(linesplit)
                    linesplit = temp_list
                
                if int(linesplit[0]) in pou_keys:
                    if int(linesplit[0]) not in wr_pod_dict.keys():
                        wr_pod_dict[int(linesplit[0])] = dict()
                    
                    wr_pod_dict[int(linesplit[0])][int(linesplit[2])] = dict()
                    wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['pod_id'] = int(linesplit[1])
                    wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['pod_use_id'] = int(linesplit[2])
                    wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['wr_type'] = linesplit[3].strip()
                    wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['source_type'] = linesplit[8].strip()
                    
                    use_code = linesplit[11].strip()
                    wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['Unmod_Use'] = use_code
                    
                    if use_code.strip() == 'PRIMARY AND SUPPLEMENTAL IRRIGATION':
                        use_code = 'IRRIGATION'
                        
                    if use_code.strip() == 'IRRIGATION AND DOMESTIC':
                        use_code = 'IRRIGATION'
                    
                    if use_code.strip() == 'DOMESTIC EXPANDED':
                        use_code = 'DOMESTIC'
                    
                    wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['use_code'] = use_code
                    
                    if len(linesplit[12]) > 5:
                        celltext = linesplit[12]
                        
                        slashid = [i for i in range(len(celltext)) if celltext[i] is '/' or celltext[i] is ':']
                        year = int(celltext[slashid[1]+1:slashid[1]+5])
                        day = int(celltext[slashid[0]+1:slashid[1]])
                        mon = int(celltext[0:slashid[0]])
                            
                        wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['priority'] = datetime.datetime(year,mon,day)
                    
                    if float(linesplit[13]) != -999:
                        wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['duty'] = float(linesplit[13])
                    #else:
                    #    wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['duty'] = float('nan')
                        
                    if float(linesplit[14]) != -999:
                        wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['rate'] = float(linesplit[14])
                    else:
                        wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['rate'] = float('nan')
                        
                    if float(linesplit[16]) != -999:    
                        wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['max_rate'] = float(linesplit[16])
                    else:
                        wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['max_rate'] = float('nan')
                    
                    if wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['max_rate'] == 0:
                        wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['max_rate'] = '-'
                    
                    if float(linesplit[17]) != -999:    
                        wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['acre_feet'] = float(linesplit[17])
                    #else:
                    #    wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['acre_feet'] = float('nan')
                        
                    wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['source'] = linesplit[20]
                    wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['tributary'] = linesplit[21]
                    wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['begin_mon'] = int(linesplit[25])
                    wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['begin_day'] = int(linesplit[26])
                    wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['end_mon'] = int(linesplit[27])
                    wr_pod_dict[int(linesplit[0])][int(linesplit[2])]['end_day'] = int(linesplit[28])
                    
            else:
                line_counter = line_counter + 1
    search.close()
    return wr_pod_dict
    
#%%
    
root_path = r'Z:\Projects\INFEWS\Modeling\GIS_Data\Data\wr_state'

pou_file = 'WR_Region_POU.csv'
pod_file = 'WR_Region_POD.csv'
pod_source = 'POD_Source_Codes.csv'

pod_source_dict = GetPOD_CODES(pod_source,root_path)

wr_pou_dict = GetWR_Data(pou_file, root_path)
wr_pod_dict = GetWR_POD_Data(pod_file, root_path, wr_pou_dict.keys())

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
        
bb = np.asarray(num_pou_nan)

#%%
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
#                if i == 106519:
#                    stpo = 0
#                    
                if i in wr_pod_dict.keys():
                    temp_sum = 0
                    num_pods = 0
                    pods_temp = []
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
                            #temp_sum[pod_source_dict[wr_pod_dict[i][iii]['source_type']]] = temp_sum[pod_source_dict[wr_pod_dict[i][iii]['source_type']]] + wr_pod_dict[i][iii]['acre_feet']
                            
                    if temp_sum != 0:
                    #if len(temp_sum.keys()) != 0:
                    #    for zi in temp_sum.keys():
                    #        if zi not in wris_per_use[wr_pou_dict[i][ii]['Use']].keys():
                    #            wris_per_use[wr_pou_dict[i][ii]['Use']][zi] = 0
                                
                        wris_per_use[wr_pou_dict[i][ii]['Use'][iii]][wr_pou_dict[i][ii]['Source']] = wris_per_use[wr_pou_dict[i][ii]['Use'][iii]][wr_pou_dict[i][ii]['Source']] + temp_sum
                        #non_empty_acres.append((i,ii,wr_pou_dict[i][ii]['Use']))
                        
                        if num_pods > 1:
                            jjj = pods_temp[len(pods_temp)-1]
                            if len(comment) == 0:
                                comment = 'Multiple PODs'
                            else:
                                comment = comment + ', Multiple PODs'
                            non_empty_acres.append((2,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii], wr_pou_dict[i][ii]['Source'],
                                                '-',num_pods,wr_pod_dict[i][jjj]['Unmod_Use'], wr_pod_dict[i][jjj]['use_code'],wr_pod_dict[i][jjj]['source_type'],wr_pod_dict[i][jjj]['wr_type'],
                                                '-',temp_sum,wr_pod_dict[i][jjj]['max_rate'],comment)) 
                        else:
                            non_empty_acres.append((2,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii], wr_pou_dict[i][ii]['Source'],
                                                '-',pods_temp[0],wr_pod_dict[i][jjj]['Unmod_Use'], wr_pod_dict[i][jjj]['use_code'],wr_pod_dict[i][jjj]['source_type'],wr_pod_dict[i][jjj]['wr_type'],
                                                '-',temp_sum,wr_pod_dict[i][jjj]['max_rate'],comment)) 
                        
                            
                        
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
                                    
                            jjj = pods_temp[len(pods_temp)-1]
                            if len(comment) == 0:
                                comment = 'Max rate (cfs)'
                            else:
                                comment = comment + ', Max rate (cfs)'
                            non_empty_acres.append((2,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii], wr_pou_dict[i][ii]['Source'],
                                                '-',pods_temp[0],wr_pod_dict[i][jjj]['Unmod_Use'], wr_pod_dict[i][jjj]['use_code'],wr_pod_dict[i][jjj]['source_type'],wr_pod_dict[i][jjj]['wr_type'],
                                                '-','-',wr_pod_dict[i][jjj]['max_rate'],comment))
                                
                else:
                    empty_acres.append((2,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],'-','No POD Match'))
                    
            else:
                
#                if i == 174985:
#                    stpo = 0
                
                if i in wr_pod_dict.keys():
                    old_duty = -999
                    num_pods = 0
                    duty_total = 0
                    pods_temp = []
                    
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
                                    print i, ii, wr_pod_dict[i][jjj]['use_code']
                                    comment = 'Mutiple values for duty'
                                    
                                if num_pods == 0:
                                    old_duty = temp
                                duty_total = duty_total + temp
                                
                            num_pods = num_pods + 1
                            pods_temp.append(jjj)
                                
                    if old_duty != -999:
                        #wris_per_use[wr_pou_dict[i][ii]['Use'][iii]][wr_pou_dict[i][ii]['Source']] = wris_per_use[wr_pou_dict[i][ii]['Use'][iii]][wr_pou_dict[i][ii]['Source']] + wr_pou_dict[i][ii]['acre_feet']*temp
                        wris_per_use[wr_pou_dict[i][ii]['Use'][iii]][wr_pou_dict[i][ii]['Source']] = wris_per_use[wr_pou_dict[i][ii]['Use'][iii]][wr_pou_dict[i][ii]['Source']] + wr_pou_dict[i][ii]['acre_feet']*old_duty
                        #non_empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Use'],wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['acre_feet']))
                        
                        if num_pods > 1:
                            jjj = pods_temp[len(pods_temp)-1]
                            if len(comment) == 0:
                                comment = 'Multiple PODs'
                            else:
                                comment = comment + ', Multiple PODs'
                                
                            non_empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],
                                                wr_pou_dict[i][ii]['acre_feet'],num_pods,wr_pod_dict[i][jjj]['Unmod_Use'],wr_pod_dict[i][jjj]['use_code'],wr_pod_dict[i][jjj]['source_type'],wr_pod_dict[i][jjj]['wr_type'],
                                                old_duty,wr_pou_dict[i][ii]['acre_feet']*old_duty,wr_pod_dict[i][jjj]['max_rate'],comment))
                        else:
                            non_empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],
                                                wr_pou_dict[i][ii]['acre_feet'],pods_temp[0],wr_pod_dict[i][pods_temp[0]]['Unmod_Use'],wr_pod_dict[i][pods_temp[0]]['use_code'],wr_pod_dict[i][pods_temp[0]]['source_type'],wr_pod_dict[i][pods_temp[0]]['wr_type'],
                                                old_duty,wr_pou_dict[i][ii]['acre_feet']*old_duty,wr_pod_dict[i][jjj]['max_rate'],comment))
                    else:
                        if len(pods_temp) > 0:

                            if 'IRRIGATION' in wr_pou_dict[i][ii]['Use'][iii]:
                                if len(comment) == 0:
                                    comment = 'Assumed Duty of 3 acre-ft'
                                else:
                                    comment = comment + ', Assumed Duty of 3 acre-ft'
                                
                                non_empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],
                                                wr_pou_dict[i][ii]['acre_feet'],pods_temp[0],wr_pod_dict[i][pods_temp[0]]['Unmod_Use'],wr_pod_dict[i][pods_temp[0]]['use_code'],wr_pod_dict[i][pods_temp[0]]['source_type'],wr_pod_dict[i][pods_temp[0]]['wr_type'],
                                                3,wr_pou_dict[i][ii]['acre_feet']*3,wr_pod_dict[i][jjj]['max_rate'],comment))
                            else:
                                if len(comment) == 0:
                                    comment = 'No Duty'
                                else:
                                    comment = comment + ', No Duty'
                                non_empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],
                                                    wr_pou_dict[i][ii]['acre_feet'],pods_temp[0],wr_pod_dict[i][pods_temp[0]]['Unmod_Use'],wr_pod_dict[i][pods_temp[0]]['use_code'],wr_pod_dict[i][pods_temp[0]]['source_type'],wr_pod_dict[i][pods_temp[0]]['wr_type'],
                                                    '-','-',wr_pod_dict[i][jjj]['max_rate'],comment))
                        else:
                            empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['acre_feet'],'No POD with same USE CODE OR WR Source'))
                else:
                    empty_acres.append((1,i,ii,wr_pou_dict[i][ii]['Unmod_Use'],wr_pou_dict[i][ii]['Use'][iii],wr_pou_dict[i][ii]['Source'],wr_pou_dict[i][ii]['acre_feet'],'No POD Match'))

non_empty_acres_arr = np.asarray(non_empty_acres)
empty_acres_arr = np.asarray(empty_acres)

#%%
#del_list = []            
#for j in range(len(empty_acres)):  
#    i = empty_acres[j][0]
#    ii = empty_acres[j][1]
#    
#    if i in wr_pod_dict.keys():    
#        temp_sum = dict()
#        for iii in wr_pod_dict[i].keys():
#            
#            if ('acre_feet' in wr_pod_dict[i][iii].keys() and wr_pou_dict[i][ii]['Use'] == wr_pod_dict[i][iii]['use_code'] and 
#                wr_pou_dict[i][ii]['Source'].strip() == pod_source_dict[wr_pod_dict[i][iii]['source_type'].strip()]):
#                if wr_pod_dict[i][iii]['source_type'].strip() not in temp_sum.keys():
#                    temp_sum[wr_pod_dict[i][iii]['source_type'].strip()] = 0
#                temp_sum[wr_pod_dict[i][iii]['source_type'].strip()] = temp_sum[wr_pod_dict[i][iii]['source_type'].strip()] + wr_pod_dict[i][iii]['acre_feet']
#        
#        if len(temp_sum.keys()) != 0:
#            for zi in temp_sum.keys():
#                tempsource = pod_source_dict[zi.strip()].strip()
#                if tempsource not in wris_per_use[wr_pou_dict[i][ii]['Use']].keys():
#                    wris_per_use[wr_pou_dict[i][ii]['Use']][tempsource] = 0
#                
#                wris_per_use[wr_pou_dict[i][ii]['Use']][tempsource] = wris_per_use[wr_pou_dict[i][ii]['Use']][tempsource] + temp_sum[zi]
#                non_empty_acres.append((2,i,ii,wr_pou_dict[i][ii]['Use'],tempsource,temp_sum[zi]))
#                
#            del_list.append(j)
#
#temp_empty_acres = empty_acres
#empty_acres = []
#for i in range(len(temp_empty_acres)): 
#    if i not in del_list:
#        empty_acres.append(temp_empty_acres[i])
#
#del temp_empty_acres     

##%%
#wrate_per_use = dict()
#empty_rates = []
#non_empty_rates = []
#
#for i in wr_pou_dict.keys():
#    if i in wr_pod_dict.keys():
#        for ii in wr_pou_dict[i].keys():
#            if np.isnan(wr_pou_dict[i][ii]['acre_feet']) == 1:
#                for iii in wr_pod_dict[i].keys():
#    
#                    if wr_pou_dict[i][ii]['Use'] not in wrate_per_use.keys():
#                        wrate_per_use[wr_pou_dict[i][ii]['Use']] = dict()
#                        
#                    if wr_pou_dict[i][ii]['Source'] not in wrate_per_use[wr_pou_dict[i][ii]['Use']].keys():
#                        wrate_per_use[wr_pou_dict[i][ii]['Use']][wr_pou_dict[i][ii]['Source']] = 0
#                        
#                    if wr_pou_dict[i][ii]['Use'] == wr_pod_dict[i][iii]['use_code'] and np.isnan(wr_pod_dict[i][iii]['max_rate']) == 0:
#                        wrate_per_use[wr_pou_dict[i][ii]['Use']][wr_pou_dict[i][ii]['Source']] = wrate_per_use[wr_pou_dict[i][ii]['Use']][wr_pou_dict[i][ii]['Source']] + wr_pod_dict[i][iii]['max_rate']
#                        non_empty_rates.append((i,iii,wr_pod_dict[i][iii]['max_rate']))
#                        
#                    elif np.isnan(wr_pod_dict[i][iii]['max_rate']) == 1:
#                        empty_rates.append((i,iii,float('nan'),wr_pou_dict[i][ii]['Use'],wr_pod_dict[i][iii]['use_code']))
#    else:
#        empty_rates.append((i,-999, float('nan'),None,None))
#  
#              
##%%
#non_empty_acres_arr = np.asarray(non_empty_acres)
#empty_acres_bool_arr = np.asarray(empty_acres)
#
#for i in range(len(empty_acres)):
#    if empty_acres[i][0] in non_empty_rates[:,0]:
#        empty_acres_bool[i,0] = 1
#    else:
#        empty_acres_bool[i,0] = 0
#
#          
##%%
#pou_pod_source = dict()
#pou_pod_use = dict()
#missing_pou_pod = []
#
#for i in wr_pou_dict.keys():
#    if i in wr_pod_dict.keys():
#        pou_pod_source[i] = []
#        pou_pod_use[i] = []
#        
#        for ii in wr_pod_dict[i].keys():
#            if wr_pod_dict[i][ii]['source_type'] not in pou_pod_source[i]:
#                pou_pod_source[i].append(wr_pod_dict[i][ii]['source_type'])
#                
#            if wr_pod_dict[i][ii]['use_code'] not in pou_pod_use[i]:
#                pou_pod_use[i].append(wr_pod_dict[i][ii]['use_code'])
#    
#    else:
#        missing_pou_pod.append(i)
#                
##            if 'IRRIGATION' in wr_pod_dict[i][ii]['use_code']
#        