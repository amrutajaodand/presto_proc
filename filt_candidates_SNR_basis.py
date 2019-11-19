"""This script filters candidates from each folder based on the SNR of the filtered candidates""" 
#################################################################################
#importing modules
import numpy as np
import os, sys, getopt, re, shutil, math, time
from datetime import datetime
import multiprocessing as mp
import ddplan_parallelizer as ddp_par
from glob import glob
import sys 
#################################################################################
##  input parameters 
### use the arguments given by the user ##########################################
args=sys.argv[1:]   # list of arguments
try:
   opts, args = getopt.getopt(args,"h",["folder_tag=","file_tag=","N_cores=","f_num="])#parameters to be passed on to the system 
# where the file with file_tag is and N_cores for number of cores to be used for the processing.

except getopt.GetoptError, exc:                 # display right format to parse arguments
   print exc.msg
   sys.exit(2)  # exit if the format of arguments is not right
print opts
#print args 
for opt, arg in opts:
        if opt in ('--folder_tag'):
                folder_tag = str(arg)
        if opt in ('--file_tag'):
                file_tag = str(arg)
        if opt in ('--N_cores'):
                N_cores = int(arg)
        if opt in ('--f_num'):
                f_num = str(arg)
#################################################################################
# import the files to be processed with PRESTO
path = "/projects/0/lotaas2/data/jaodand/GBT14B-036_AMXPs/%s"%folder_tag #remember to change the folder number ot a tag passed by the run script
infilename = "%s"%file_tag
filename,fileext = os.path.splitext(infilename)

##Initialise the folders to store the results of the processing 
#first initialise the key to name and call files 
key =filename #used for naming all the files removes last number from the observation 

#create folder for a given observation date
outpath= "/projects/0/lotaas2/data/jaodand/results/multicore"

#create folder specific to the particular data on a given date(results_dir)
results_dir = os.path.join(outpath,path[50:],key)
print results_dir

os.chdir(results_dir)
os.system("awk '$3 > 10.0' *candidates*.txt > candidates_SNR_filtered.txt")

##getting the directories list based on above foiltered candiates
dir_list =  []  
h_cands_SNR_filt = open("candidates_SNR_filtered.txt",'rw')

for line in h_cands_SNR_filt:
	sep = ('_A')
	filename1 =  (line.split()[0]).split(sep,1)[0]
	filename1 = filename1.split("out_",1)[-1]
	filename2 = ''.join([filename,"_",filename1])
	print filename2
	dir_list.append(filename2)
print dir_list


#Now form subfolders inside the result folders to store only the essential results with threshold cuts there 
for i,j in enumerate(glob("*ephmer*")):
	os.chdir(j)

	print os.getcwd()
	if not os.path.exists("%s_cands_SNR_filtered"%j):
		os.mkdir("%s_cands_SNR_filtered"%j)
	for l,m in enumerate(dir_list):	
		os.system("ls *%s*.tar.gz | xargs -i cp {}  %s_cands_SNR_filtered"%(m,j))

	os.chdir("../")

