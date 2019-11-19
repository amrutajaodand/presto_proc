import numpy as np
import os,sys


folder_list = ['20140721','20140721','20140809','20140809','20140811','20140811']
tags_list= ['guppi_56859_SwiftJ1749.4-2807_0002_0001.fits','guppi_56859_XTEJ1814-338_0003_0001.fits','guppi_56878_IGRJ00291+5934_0003_0001.fits','guppi_56878_SAXJ1808.4-3658_0001_0001.fits','guppi_56880_IGRJ17498-2921_0002_0001.fits','guppi_56880_IGRJ17511-3057_0001_0001.fits']
for i in range(0,len(folder_list)):
	print 'moving over to the file  %s_%s'%(folder_list[i],tags_list[i])
	os.system('nohup python presto_proc_05122014.py --folder_tag %s --file_tag %s > log%s.log &'%(folder_list[i],tags_list[i],tags_list[i]))

	#print 'nohup python presto_proc_03122014.py --folder_tag %s --file_tag %s > log%s.log'%(folder_list[i],tags_list[i],tags_list[i])
