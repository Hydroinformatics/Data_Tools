# -*- coding: utf-8 -*-

import os, sys, ftplib, logging, urlparse, zipfile, shutil

#%%
# FTP object used to fetch data from an FTP remote server
# Inputs:
# 1- FTP url - full path to the location of the dir in the remote server
# 2- Destination dir - dir where data will be downloaded/saved to
class FTP_Parser:
    
    def __init__(self, url = None, destdir = os.getcwd()):     
        self.ftp = ftplib.FTP()
        
        if url is None:
            sys.exit("Must provide an URL FTP path to the desire location. None was provided")
        
        self.url = urlparse.urlparse(url)
        
        self.dest_dir = destdir.replace('\\','/')
        if self.dest_dir[-1] != '/':
            self.dest_dir = self.dest_dir + '/'
            
        self.dl_data_list = []

   
    def list_dir(self):
        dir_listing = []
        self.ftp.dir(lambda x: dir_listing.append(x))
        return [(line[0].upper() == 'D', line.rsplit()[-1]) for line in dir_listing]
    
    def chdir(self):
        dir_levels = self.url.path.split('/')
        for ilevel in dir_levels:
            if ilevel == '':
                temp_dirs = self.list_dir()
            else:
                if (True, ilevel) in temp_dirs:
                    self.ftp.cwd(ilevel)
                    temp_dirs = self.list_dir()
    
    def Close(self):
        self.ftp.quit()
        
    # Creates a connection with the remote server
    #Inputs:
    # -- 1. Port,  -- 2. Username, -- 3. Password           
    def Connect(self, port = 0, user='anonymous', password = ''):
        
        server_ip = self.url.netloc
        
        try:
            self.ftp.connect(server_ip, port)
            #user and password vars must be in the config file
            self.ftp.login(user=user,  passwd=password)
            logging.info("Connecting to ftp server...")
            
        except:
            logging.error("Failed to setup connection")
            raise

    # Fetches and downloads the required data to the Destination dir 
    def Fetch(self, flist):
        try:
            self.chdir()
        except:
            logging.error('Dir path:' + self.url.path.split + ' not found.')
            raise
            
        self.server_files = self.list_dir()
        
        if flist:
            self.server_files = [(sbool, sfile) for sbool, sfile in self.server_files if sfile in flist]
            
        for sbool, sfile in self.server_files:
            
            local_file = os.path.join(self.dest_dir, sfile)
            
            try:
                print('Downloading: ' + sfile)
                self.ftp.retrbinary('RETR %s' %sfile, open(local_file, 'w').write)
                self.dl_data_list.append(sfile)
            except:
                logging.error("Failed to download " + sfile)
                raise
  
#%%
class DataManager():
    
    def __init__(self):
        self.dest_dir = None
        self.zip_dir = None
        self.data_path = None
        self.zipf_list = None
        self._data_url = None
        
        #self.pathuzip = None

    def Download(self, destdir, zip_dir = 'Zip_Data', flist = None, del_destdir = 0):
        
        self.dest_dir = destdir
        
        zip_dir = destdir + '/' + zip_dir
        
        if os.path.isdir(zip_dir) and del_destdir == 1:
            shutil.rmtree(zip_dir)
        
        if not os.path.isdir(zip_dir):
            os.makedirs(zip_dir)
    
        ftp = FTP_Parser(self._data_url, zip_dir)
        ftp.Connect()
        ftp.Fetch(flist)
        ftp.Close()
        
        self.zip_dir = zip_dir
        self.zipf_list = ftp.dl_data_list
        #self.pathuzip = destdir
        
    def Unzip_Data(self, pathuzip = None, datapath = None):
        
        if datapath:
            self.zip_dir = datapath
            self.zipf_list = [fname for fname in os.listdir(datapath) if fname[-3:] == 'zip']
            
        if pathuzip:
            self.pathuzip = pathuzip
            
        #if not pathuzip:
        
        for zipf in self.zipf_list:
            try:
                print self.pathuzip + '/' + zipf
                with zipfile.ZipFile(self.zip_dir + '/' + zipf, "r") as z:
                    z.extractall(self.pathuzip + '/' + zipf.strip('.zip'))
                
            except:
                logging.error("Failed to unzip file " + zipf)
                raise