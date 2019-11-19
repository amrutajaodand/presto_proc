"""This script can be used to obtain a ddplan for parallelizing the script 
Author : Amruta 
Created: 09/01/2015

"""
#################################################################################
#importing modules
import numpy as np
import os, sys, getopt, re, shutil, math
from datetime import datetime
#################################################################################
#################################################################################
def ddplan_par(ncore,results_dir,ddplan_file_tag):
	print ncore, results_dir,ddplan_file_tag
	h = open(os.path.join(results_dir,ddplan_file_tag),'rw')
	num_lines = sum(1 for line in h) #number of subbands

	i = 0 
	lowdm = 0 

	h = open(os.path.join(results_dir,ddplan_file_tag),'rw')
	f = open(os.path.join(results_dir,"%s_parallel"%ddplan_file_tag),'w+')

	for i, line in enumerate(h,1):        #i counter for no. of subbands,  
		if i ==1:
			f.write(line)	      #Low DM    High DM     dDM  DownSamp  dsubDM   #DMs  DMs/call  calls  WorkFract

		if i in range(2,num_lines-2):   # start from 2 as the line number 1 is for header end at num_lines-2 as last three are 3 empty lines at end. 
			temp = line.split()	# Getting a given line from he ddplan 
			dividend = int(temp[7])/ncore # number of call being distributed per core temp[7] - number of , here there will be some calls left as remainder which we will accommodate in second looping. 
			nDM_core = dividend * int(temp[6]) #number of DMs/core are being obtained no.of calls/core*temp[6](DMs/cal)
			print dividend,nDM_core          
		
			for j in range(1,ncore+1):
				highdm = float(temp[0])+(j*float(temp[2])*nDM_core)# Now, setting a limit for the highDM on given core, int(temp[2]) is the dDM.  
				lowdm = float(temp[0])+((j-1)*float(temp[2])*nDM_core)#Setting the lowDM for a given core. 
				f.write("%f    %f     %f  %d  %f   %d  %d  %d  %f  %d\n"%(lowdm,highdm,float(temp[2]),int(temp[3]),float(temp[4]),nDM_core,int(temp[6]),dividend,float(temp[8]),j))
				print lowdm, highdm, j 

			#Strating the process for the leftover calls form above and distributing them evenly through all the cores. 
			highdm_core_end = highdm 		
			remainder = int(temp[7])-int(dividend*ncore) #number of calls remaining 

	#		remDM_core = (remainder*int(temp[6]))/ncore #numbers of DMs reamianing per core obtained by remaining calls*DMs/cals*(1/ncore)#take care that this is divisible by ncore!!! 
			print highdm_core_end, remainder
			# now, since number of calls < ncore therefore we need to allot one call per core by force thus dsubDM become the difference between lowDM and highDM
			if remainder ==0: 
				print "entered the loop here nothing to do here"

			else:
				p = (float(temp[1])-highdm_core_end)/float(remainder)  #Care should be taken to check if the (float(temp[1]) i.e highDM in previous ddplan - highdm_core_end) is divisible by dsubDM. Thus, I adopt a reverse method to see if the number obtained through this subtraction divided by remainder calls, yeilds the dsubDM	
				print round(p,16) 
				if round(p,4) == round(float(temp[4]),4):	

					print "entered the loop here"
					for j in range(1,remainder+1):#ncore+1):
						highdm = highdm_core_end+(j*float(temp[2])*int(temp[6]))#(j*float(temp[2])*remDM_core)
						lowdm = highdm_core_end+((j-1)*float(temp[2])*int(temp[6])) #((j-1)*float(temp[2])*remDM_core)
						print lowdm, highdm,j
						numcall = 1
						f.write("%f    %f     %f  %d  %f   %d  %d  %d  %f  %d\n"%(lowdm,highdm,float(temp[2]),int(temp[3]),float(temp[4]),int(temp[6]),int(temp[6]),numcall,float(temp[8]),j))#remDM_core,int(temp[6]),remainder,float(temp[8])))

	f.close()
	return
