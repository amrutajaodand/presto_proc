import numpy as np
import os,sys, time



os.system(


"""
tags_list= ['a','bb','cc','ddd','b','s','t','f','f','l','d']
n = 2	# number of cores
k = 0 	# vary k between 0 and int(tags_list/n)

for i in range(k*n , (k+1)*n):
	if i < len(tags_list):
		tag = tags_list[i]
		print 'generating ass %s'%tag
		os.system('nohup python script.py --file_tag %s > log%s.log &'%(tag,tag))
"""
