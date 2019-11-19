"""This script filters candidates from each folder based on the SNR of the filtered candidates""" 
#################################################################################
#importing modules
import numpy as np
import os, sys, getopt, re, shutil, math, time
from datetime import datetime
import multiprocessing as mp
import ddplan_parallelizer as ddp_par
from glob import glob
import sys,tarfile 
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

##getting the directories list based on above foiltered candiates
dir_list = glob("prepfold_ephmer_*")
print dir_list  

datafile = open('plot_data_%s.txt'%key, 'w')
datafile.write("#Plt_title\t#DM\t#DELTA_TASC\t#RED_CHI_SQ\n")

for i,j in enumerate(dir_list):
	print os.getcwd()
	os.chdir(j)

	##making a directory to store candidates greater than 1.5 red_chi_sqr 
	if not os.path.exists("Candidates_red_chi_gt_1.75"):		
		os.mkdir("Candidates_red_chi_gt_1.75")
	
	##opening a tarred subfolder to read bestprof files 
	for tarf in glob("*.tar.gz"):
		tar = tarfile.open(tarf)

		###getting names of the bestprof files in the tar subfolder 
		for fl in tar.getnames():
			print "here we are",fl
			##working to get info from the .bestprof files 
			if '.bestprof' in fl:
				f_bprof = fl.split('/')[-1]	## f_bprof is the name of the file 
				f_ps = f_bprof[:-9]+'.ps' ##fps is the .ps file corresponding to f_bprof  
				print "f_bprof is",fl,f_bprof
				print "f_ps is",f_ps
				####getting a reduced-chi_sq, dm and tasc from a .bestprof file for plotting
				file_bestprof_data = (tar.extractfile(fl)).read()
				Red_chi_sq_line = [line for line in file_bestprof_data.split('\n') if "Reduced" in line]
				red_chi_sq = float(Red_chi_sq_line[0].split()[-1])
				dm = float(re.search('DM(.+?)_par', f_bprof).group(1))
				if '2015' in f_bprof: 
					delta_tasc = 0.1*float(re.search('par_2015_(.*?)_PSR', f_bprof).group(1)) 
				else:
					delta_tasc = 0.1*float(re.search('par_(.*?)_PSR', f_bprof).group(1)) 
				print red_chi_sq,dm,delta_tasc

				datafile.write("%s\t%0.2f\t%0.2f\t%f\n"%(re.search('(.+?)_DM', tarf).group(1),dm,delta_tasc,red_chi_sq))

				if red_chi_sq >=1.75:
					print "Candidate with reduced chi sq > than 1.5 %s"%(fl)

					print "tar -xvf %s -C Candidates_red_chi_gt_1.75 %s"%(tarf,tarf[:-7]+"/"+f_bprof) 				
					os.system("tar -xvf %s -C Candidates_red_chi_gt_1.75 %s"%(tarf,tarf[:-7]+"/"+f_bprof)) 				

					print "tar -xvf %s -C Candidates_red_chi_gt_1.75 %s"%(tarf,tarf[:-7]+"/"+f_ps) 				
					os.system("tar -xvf %s -C Candidates_red_chi_gt_1.75 %s"%(tarf,tarf[:-7]+"/"+f_ps)) 				

	os.chdir("../")
datafile.close()
