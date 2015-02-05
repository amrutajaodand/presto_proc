import numpy as np
import os,sys


folder_list = ['20140821']
tags_list= ['guppi_56890_SwiftJ1749.4-2807_0002_0001.fits','guppi_56890_SwiftJ1749.4-2807_0002_0002.fits','guppi_56890_SwiftJ1749.4-2807_0002_0003.fits']

for folder in folder_list:
	print 'folder is %s'%folder
	for tag in tags_list:
		print 'moving over to the file  %s_%s'%(folder,tag)
		os.system('python presto_proc_03122014.py --folder_tag %s --file_tag %s > log%s.log'%(folder,tag,tag))
