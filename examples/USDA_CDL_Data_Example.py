# -*- coding: utf-8 -*-

import os, sys
os.chdir('..')
sys.path.append(os.getcwd() + '\tools')

from tools.USDA_Data import USDA_CDL_Parser as usda_parser

destdir= r'C:\Users\sammy\Documents\GitHub\Data_Tools\examples\Download_test'

# Fetch data from Remote Server
parser = usda_parser.Parser()
#parser.Download(destdir, [2019], 1)
parser.Unzip_Data(destdir + '\CDL_Zip_Data')
