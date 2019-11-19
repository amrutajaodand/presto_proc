""" script to run the presto proc over three sources for whom rfifind works. 
Modified on : 16/01/2015
Modifications(newest to oldest):
The script was modified to run over multiple sources  
The script prestp_proc_05122014 during this run was modified to have .fft passed to accelsearch and 
include effect of downsampling on prepsubband
Editor: Amruta Jaodand
"""
import numpy as np
import os,sys

N_cores = 8 # Number of cores used for processing 
folder_list = ['20140811']#['20140721']#,'20140809','20140811']
tags_list=['guppi_56880_IGRJ17511-3057_0001_0001.fits']#['guppi_56859_SwiftJ1749.4-2807_0002_0001.fits']#,'guppi_56878_SAXJ1808.4-3658_0001_0001.fits','guppi_56880_IGRJ17511-3057_0001_0001.fits']


for i in range(0,len(folder_list)):
	filename,fileext = os.path.splitext(tags_list[i])
	print 'moving over to the file  %s_%s'%(folder_list[i],tags_list[i])
	os.system('nohup python presto_proc_multicore.py --folder_tag %s --file_tag %s --N_cores %d > log_16012015_multi_%s.log &'%(folder_list[i],tags_list[i],N_cores,filename))
	
	print 'nohup python presto_proc_multicore.py --folder_tag %s --file_tag %s --N_cores %d > log_16012015_multi_%s.log &'%(folder_list[i],tags_list[i],N_cores,filename)

