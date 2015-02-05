import numpy as np
import os, sys, getopt
from datetime import datetime
import re


##  input parameters 
### use the arguments given by the user #################################################################################
args=sys.argv[1:]   # list of arguments
try:
   opts, args = getopt.getopt(args,"h",["folder_tag=","file_tag="])#, "m=","incl_angle="])

except getopt.GetoptError, exc:			# display right format to parse arguments
   print exc.msg
   sys.exit(2)  # exit if the format of arguments is not right
print opts
#print args 
for opt, arg in opts:	
	if opt in ('--folder_tag'):
		folder_tag = str(arg)
	if opt in ('--file_tag'):			
		file_tag = str(arg)
  

printfile = 'test_%s_%s.dat'%(folder_tag,file_tag)
fw= open(printfile,'w+')
fw.write('testing  %s_%s'%(folder_tag,file_tag))
