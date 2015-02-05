import numpy as np
import os,sys


folder_list = ['20140809']
tags_list= ['a','bb','cc','ddd']
for folder in folder_list:
	for tag in tags_list:
		print 'generating ass %s_%s'%(folder,tag)
		os.system('nohup python ../script.py --folder_tag %s --file_tag %s > log%s.log &'%(folder,tag,tag))
