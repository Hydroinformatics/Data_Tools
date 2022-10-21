# -*- coding: utf-8 -*-

import os, sys
os.chdir('..')
sys.path.append(os.getcwd() + '\tools')

from tools.PRISM_Data import PRISM_Parser as prism_parser

destdir= r'...\Data_Tools\examples\Download_test'

# Fetch data from Remote Server
parser = prism_parser.Parser()
parser.Download(destdir)
parser.Unzip_Data()
