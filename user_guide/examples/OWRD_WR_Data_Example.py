# -*- coding: utf-8 -*-
"""
Created on Fri May 29 11:20:31 2020

@author: sammy
"""
import os, sys
os.chdir('..')
sys.path.append(os.getcwd() + '\tools')

from tools.OWRD_Data import OWRD_WR_Parser as wr_parser

#%% Initialize Parser
wr_parser = wr_parser.Parser()

#%% Download and Unzip data
#dest_dir = r'..\Data_Tools\examples\Download_test'
#wr_parser.Download(dest_dir)
#wr_parser.Unzip_Data()

#%% Clip POU and POD shapefile
#boundary_shp = r'...\UmatillaBoundary_proj.shp'
boundary_shp = r'...\subs_input.shp'
wrdata_dir = r'...\wr_state_sh'
wr_parser.Clip(['uma'], boundary_shp, wrdata_dir)

#%%
pou_file_noduplicate = 'wr_v_pou_public_proj_REGION_NoDuplicatesPy_V2'
feature_dict, remove_features = wr_parser.EliminateDuplicatePOUs(pou_file_noduplicate)

wr_pou_dict = wr_parser.GetWR_POU_Data()
wr_pod_dict = wr_parser.GetWR_POD_Data()

non_empty_acres, empty_acres = wr_parser.POU_To_POD_Matcher()

#non_empty_acres_arr = np.asarray(non_empty_acres)
#empty_acres_arr = np.asarray(empty_acres)
#
#hru_wr_area_tresh = 5   
    
