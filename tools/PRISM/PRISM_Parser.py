# -*- coding: utf-8 -*-

from tools.utils import data_mgt_utils as dmgt

#%%
class Parser:
    
    def __init__(self):
        self.data_mgt_obj = dmgt.DataManager()
        self.data_url = 'ftp://prism.nacse.org/'
        self.zip_dir = 'PRISM_Zip_Data'
    
    # --2. Class of data: '4km' (4km normals), '800m' (800m normals), 'monthly' 
    #      (monthly data), 'daily' (daily data), 'archive' (data_archive)
    # --3. Variables: 
    # ---- ppt (precipitation)
    # ---- tmin, tmean, tmax (min/mean/max temperature)
    # ---- tdmean ( mean dewpoint temperature), vpdmin, vpdmax (min/max vapor pressure deficit
    
    def Download(self, destdir, data_class = 'daily', var = '', spec_years = None, del_destdir = 0):
        
        if data_class != 'archive':
        
        else:
        
        self.data_mgt_obj._data_url = self.data_url + 
        self.data_mgt_obj.Download(destdir, self.zip_dir)
    
    def Unzip_Data(self):
        self.data_mgt_obj.Unzip_Data()
        


#
##%%
#if __name__ == '__main__':
###%% Parse Path to DEM File
#    parser = argparse.ArgumentParser(description='Inputs File')
#    parser.add_argument('path', metavar='-p', type=str, nargs='+',
#                        help='Path to text file with list of input parameters')
#    
#    args = parser.parse_args()
#    main_path = args.path[0].replace('\\','/')
#    csv_file_name = args.path[1].replace('\\','/')
#    system_boundary = args.path[2].replace('\\','/')
#    begin_year = int(args.path[3])
#    end_year = int(args.path[4])
#    
#    #main_path = r'Z:\Projects\INFEWS\Modeling\FEW_Data\PRISM_RecentYears'
#    #system_boundary = 'System_Boundary_16km.tif'
#    #begin_year = 1981
#    #end_year = 2019
#
#    pathunzip = main_path + '/temp'
#    
#    #%%
#    path = main_path + '/' + system_boundary
#    region, region_NoData, region_obj = Read_Raster(path)
#    
#    centroids_X, centroids_Y = GetPixelCentroids(region_obj)
#    
#    top_row, last_row = Raster_row_boundaries(region)
#    left_col, right_col = Raster_col_boundaries(region)
#    region = region[top_row:last_row+1,left_col:right_col+1]
#    
#    PRISM_centroids = np.zeros((len(region.flatten()),3))
#    
#    n=0
#    for i in range(top_row,last_row+1):
#        for j in range(left_col,right_col+1):
#            PRISM_centroids[n,0] = n 
#            PRISM_centroids[n,1] = centroids_X[0,j]
#            PRISM_centroids[n,2] = centroids_Y[i,0] 
#            
#            n = n + 1
#        
#    np.savetxt(main_path + '/' + 'PRISM_Raster_Centroids.csv', PRISM_centroids, delimiter=",")
#    
#    #%%
#    nrows = 0
#    for year in range(begin_year,end_year+1):
#        path = main_path + '/' + str(year)
#        fnames = os.listdir(path)
#        nrows = nrows + len(fnames)
#        
#    ncols = len(region.flatten()) + 3 # Year, Mon, Day
#    
#    #%%
#    PRISM_Data = np.ones((nrows,ncols))*-999
#    nrow = 0
#    for year in range(begin_year,end_year+1):
#    
#        path = main_path + '/' + str(year)
#        fnames = os.listdir(path)
#        
#        for name in fnames:
#            #if len(name) < 45: # added for daily temperature values
#            if len(name) < 41: # added for daily min temp.
#                dash_index = [i for i in range(len(name)) if name[i] == '_']
#                datestr = name[dash_index[-2]+1:dash_index[-1]]
#                PRISM_Data[nrow,0] = int(datestr[0:4])
#                PRISM_Data[nrow,1] = int(datestr[4:6])
#                PRISM_Data[nrow,2] = int(datestr[6:])
#                
#                pathzip = path + '/' + name
#                UnzipData(pathzip,pathunzip)
#                
#                temp_fnames = os.listdir(pathunzip)
#                main_bil_file = [tname for tname in temp_fnames if '.bil' in tname[-4:]]
#                print main_bil_file[0]
#                data = ReadBilFile(pathunzip + '/' + main_bil_file[0])
#                data = data[top_row:last_row+1,left_col:right_col+1]
#                data = np.transpose(data.flatten())
#                PRISM_Data[nrow,3:3+len(data)] = data
#                nrow = nrow + 1
#    
#    np.savetxt(main_path + '/' + csv_file_name, PRISM_Data, delimiter=",")