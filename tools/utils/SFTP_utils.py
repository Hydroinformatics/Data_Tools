# -*- coding: utf-8 -*-

import os, sys, pysftp, logging, urlparse

#%%

class SFTP_Parser:
    
    def __init__(self, url = None, user='anonymous', password = '', destdir = os.getcwd()):     
        
        if url is None:
            sys.exit("Must provide an URL FTP path to the desire location. None was provided")
        
        self.url = urlparse.urlparse(url)
        self.dest_dir = destdir.replace('\\','/')
        if self.dest_dir[-1] != '/':
            self.dest_dir = self.dest_dir + '/'
        
        try:
            cnopts = pysftp.CnOpts()
            cnopts.hostkeys = None 
            self.sftp = pysftp.Connection(self.url.netloc, username=user, password=password, cnopts=cnopts)
            logging.info("Connecting to ftp server...")
        except:
            logging.error("Failed to setup connection")
            raise
    
    def list_dir(self):
        dir_listing = []
        self.sftp.listdir_attr(lambda x: dir_listing.append(x))
        return [(line[0].upper() == 'D', line.rsplit()[-1]) for line in dir_listing]
    
#    def ftp_chdir(self):
#        dir_levels = self.url.path.split('/')
#        for ilevel in dir_levels:
#            if ilevel == '':
#                temp_dirs = self.ftp_dir()
#            else:
#                if (True, ilevel) in temp_dirs:
#                    self.ftp.cwd(ilevel)
#                    temp_dirs = self.ftp_dir()
#
#    def Download(self):
#        try:
#            self.ftp_chdir()
#        except:
#            logging.error('File path:' + self.url.path.split + ' not found.')
#            raise
#            
#        self.server_files = self.ftp_dir()
#        for sbool, sfile in self.server_files:
#            local_file = os.path.join(self.dest_dir, sfile)
#            print('Downloading: ' + local_file)
#            self.ftp.retrbinary('RETR %s' %sfile, open(local_file, 'w').write)
#    
#    def ftp_quit(self):
#        self.ftp.quit()
        
#%%
destdir=r'C:\Users\sammy\Documents\GitHub\Data_Tools\QGIS_Tools\Download_test'
url = 'ftp://ftp.nass.usda.gov/download/res'

ftp = SFTP_Parser(url, destdir)
aa= ftp.sftp_dir()