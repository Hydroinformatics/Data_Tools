# -*- coding: utf-8 -*-

from tools.utils import data_mgt_utils as dmgt

class Parser:
    
    def __init__(self):
        self.data_mgt_obj = dmgt.DataManager()
        self.data_mgt_obj._data_url = 'ftp://ftp.nass.usda.gov/download/res'
        self.zip_dir = 'CDL_Zip_Data'
        
        self.spec_years = None
        self.flist = None
    
    def Download(self, destdir, spec_years = None, del_destdir = 0):
        if spec_years:
            self.flist = []
            for year in spec_years:
                self.flist.append(str(year) + '_30m_cdls.zip')
            
        self.data_mgt_obj.Download(destdir, self.zip_dir , self.flist, del_destdir)
    
    def Unzip_Data(self , datadir = None, uzip_path = None):
        self.data_mgt_obj.Unzip_Data(datadir, uzip_path)
        
