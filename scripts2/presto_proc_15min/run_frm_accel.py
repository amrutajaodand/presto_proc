import numpy as np
import os,sys


folder_list = ['20140721','20140809','20140811']
tags_list= ['guppi_56859_SwiftJ1749.4-2807_0002_0001.fits','guppi_56878_SAXJ1808.4-3658_0001_0001.fits','guppi_56880_IGRJ17511-3057_0001_0001.fits']


for i in range(0,len(folder_list)):
	filename,fileext = os.path.splitext(tags_list[i])
	print 'moving over to the file  %s_%s'%(folder_list[i],tags_list[i])
	os.system('nohup python presto_proc_05122014_frm_accel.py --folder_tag %s --file_tag %s > log%s.log &'%(folder_list[i],tags_list[i],filename))

	#print 'nohup python presto_proc_05122014_frm_accel.py --folder_tag %s --file_tag %s > log%s.log'%(folder_list[i],tags_list[i],filename)
