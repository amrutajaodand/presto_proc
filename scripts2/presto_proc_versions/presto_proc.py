"""This script can be used to runs all presto commands over a data
A.Jaodand, 23/10/2014
Last modification: 
"""
#################################################################################
#Deciding the input parameters 

highDM = 500.0 # Highest DM value
nsub = 32	#number of subbands 
r = 0.5		#acceptable resolution in time
N = 30 		# Integer value which will serve as upper and lower limit for Tasc range 
dt = 0.1	#Step size for varying Tasc 
#################################################################################
import numpy as np
import os
import re 
import math
import shutil 
#################################################################################
# import the files to be processed with PRESTO
path = "/lustre/pulsar/scratch/ajaodand/GBT14B-036_AMXPs/20140721"
infilename = "guppi_56859_M28A_0001_0001.fits"
infilepath = os.path.join(path,infilename)
print infilepath	

# folder where the script is located 
src_dir = "/lustre/pulsar/scratch/ajaodand/GBT14B-036_AMXPs/scripts"

##Initialise the files to store the results of the processing 
#create folder for a given date
outpath= "/lustre/pulsar/scratch/ajaodand/GBT14B-036_AMXPs/results"
outdir = os.path.join(outpath,"20140721")
if not os.path.exists("%s"%outdir):
	os.mkdir("%s"%outdir)

#create folder specific to the particular data on a given date(results_dir)
filename,fileext = os.path.splitext(infilename)
print filename
results_dir = os.path.join(outpath,"20140721",filename)
if not os.path.exists("%s"%results_dir):
	os.mkdir("%s"%results_dir)

#create folder to store final complete raw data folding results
if not os.path.exists(os.path.join(results_dir,"prepfold_raw_data")):
	os.mkdir(os.path.join(results_dir,"prepfold_raw_data"))

#file containing the latest ephimeridis
dir_ephmer = "/lustre/pulsar/scratch/ajaodand/GBT14B-036_AMXPs/20140809/ephimeridis"
f_ephmer="/lustre/pulsar/scratch/ajaodand/GBT14B-036_AMXPs/20140809/SAXJ1808.4-3658.par"
#################################################################################
##define required functions
def roundup(x):
#	print x
	return (int(math.ceil(float(x)/(10**(len(str(x))-2)))))*(10**(len(str(x))-2))

def input_params(f,y,g):     
    for line in f:
             if y in line:
                     temp = line.split();
                    # print temp
    	     	     return temp[g]  
#################################################################################
####Start the processing with PRESTO 

##readfile and store metdata
os.system("readfile %s > %s"%(infilepath,os.path.join(results_dir,"metdata_%s"%filename)))

key = "GBT"  #used for naming all the files
key2 = "GBT_prep" #used for naming the .dat files from the prepsubband 
print key2 
f = open("/lustre/pulsar/scratch/ajaodand/GBT14B-036_AMXPs/results/20140721/guppi_56859_M28A_0001_0001/metdata_guppi_56859_M28A_0001_0001",'rw')
numout = roundup(input_params(f,'Spectra per file',4))
"""
##find rfi and move the o/p files to results_dir
os.system("rfifind -time 2.0 -o %s %s"%(key,infilepath))
os.system("mv %s/*%s* %s"%(src_dir,key,results_dir))

##looking for persistent low level rfi/prepdata
#First obtaining the numout from metdata file
f = open("/lustre/pulsar/scratch/ajaodand/GBT14B-036_AMXPs/results/20140721/guppi_56859_M28A_0001_0001/metdata_guppi_56859_M28A_0001_0001",'rw')
numout = roundup(input_params(f,'Spectra per file',4))

#now running the prepadata command
#print "the loc directory is %s/%s filename"%(results_dir, key)
os.system("prepdata -nobary -o %s_topo_DM0.00 -dm 0.0 -mask %s/%s_rfifind.mask -numout %d %s > prepdata_out"%(key,results_dir,key,numout,infilepath))

#exploring the .dat file and obtaining the fft
#os.system("exploredat %s_topo_DM0.00.dat"%(key))
os.system("realfft %s_topo_DM0.00.dat"%(key))

#Running the accelsearch 
os.system("accelsearch -numharm 4 -zmax 0 %s_topo_DM0.00.dat > accelsearch_out"%(key))

###Don't forget to add the .birds steps here!!!


#DDplan
##Obtaining the input arguments for DDplan.py 

f = open(os.path.join(results_dir,"metdata_%s"%filename),'rw')
nchan = int(input_params(f,'Number of channels',4))	#number of channels

f = open(os.path.join(results_dir,"metdata_%s"%filename),'rw')
dbw = float(input_params(f,'Total Bandwidth',4))	   	#Channel width 

f = open(os.path.join(results_dir,"metdata_%s"%filename),'rw')
dt = float(input_params(f,'Sample time',4))
dt = dt*10**(-6)					#Sample time 

f = open(os.path.join(results_dir,"metdata_%s"%filename),'rw')
f0 =  float(input_params(f,'Central freq',4))		#Central frequency 

##Running DDplan.py
os.system("DDplan.py -d %d -n %d -b %f -t %f -f %f -s %d -r %f -o %s_ddplan> %s/ddplan_out_%s"%(highDM,nchan,dbw,dt,f0,nsub,r,results_dir,results_dir,filename))
g = open(os.path.join(results_dir,"ddplan_out_%s"%filename),'rw')
contents1, sentinel, contents2 = g.read().partition("Low DM    High DM     dDM  DownSamp  dsubDM   #DMs  DMs/call  calls  WorkFract\n") 
#splitting the o/p of DDplan to obtain the information about different runs of prepsubbabnd

##saving DDplan step to be used for prepsubband to a file
h = open(os.path.join(results_dir,"ddplan_%s"%filename),'w+') 
h.write(sentinel) #writing the header of the file -- Low DM    High DM     dDM  DownSamp  dsubDM   #DMs  DMs/call  calls  WorkFract
h.write(contents2) # ddplan steps for subbands 
h.read()
h.close()
#Prepsubband
##Executing  prepsubband for every subband for corresponding ddm and number of calls

h = open(os.path.join(results_dir,"ddplan_%s"%filename),'rw')
num_lines = sum(1 for line in h) #number of subbands

i = 0

h = open(os.path.join(results_dir,"ddplan_%s"%filename),'rw')
for i, line in enumerate(h,1):        #i counter for no. of subbands,  
	if i in range(2,num_lines-2):   # start from 2 as the line number 1 is for header end at num_lines-3 as last three are 4 empty lines at end. 
#		print i 
		temp = line.split()
#		print temp 
		for j in range(0,int(temp[7])):
			lodm = float(temp[0])+(j*float(temp[2])*float(temp[6]))
	#		print "prepsubband -nsub %d -lodm %f -dmstep %f -numdms %d -numout %d -downsamp %d -mask %s/%s_rfifind.mask -o %s %s"%(nsub,lodm,float(temp[2]),float(temp[6]),numout,int(temp[3]),results_dir,key,key2,infilepath)
			os.system("prepsubband -nsub %d -lodm %f -dmstep %f -numdms %d -numout %d -downsamp %d -mask %s/%s_rfifind.mask -o %s %s"%(nsub,lodm,float(temp[2]),float(temp[6]),numout/int(temp[3]),int(temp[3]),results_dir,key,key2,infilepath))

	print "did not enter loop"

#mv all the recently generated .dat files to the results_dir 
os.system("mv *%s* %s"%(key,results_dir))
os.system("mv *%s* %s"%(key2,results_dir))
os.system("mv *%s* %s"%('out',results_dir))


#Prepare to search the data
##Use xargs to fft and zap the data 
os.system("ls %s/*.dat | xargs -n 1 realfft > %s/realfft_out"%(results_dir,results_dir))
#!!!Add zapping commands later 

#Searching for periodic signals 
os.system("ls %s/*.dat | xargs -n 1 accelsearch -zmax 0 > %s/accelsearch_data"%(results_dir,results_dir))

#Running single pulse search 
os.system("single_pulse_search.py %s/*.dat > %s_single_pulse_search_out"%(results_dir,key))

#sifting periodic candidates 
shutil.copy('ACCEL_sift.py',results_dir)
os.chdir("%s"%results_dir)
os.system("python ACCEL_sift.py > cands.txt")
#os.chdir("/lustre/pulsar/scratch/ajaodand/GBT14B-036_AMXPs/scripts")


os.chdir("%s"%results_dir) #changing the home directory for future operations to reduce the path length 

# we need to save all the accelcands from the cands file generated above for next prepfold steps
g = open('cands.txt','rw') 
contents1, sentinel, contents2 = g.read().partition("                           file:candnum                               DM     SNR    sigma   numharm    ipow     cpow      P(ms)          r         z     numhits \n")

h = open("cands_filt.txt",'w+') 
h.write(sentinel)
h.write(contents2)
h.close()
f = open("cands_filt.txt",'r')      #the part of the cands file containing the candidates are stored 
m = open("candidates.txt",'wa+')      #from the cands_filt file only the candidates which will be used for prepfold are stored here 
m.write(sentinel)
for line in f:
	if key2 in line:
		 m.write(line)
m.close()

os.chdir("%s"%results_dir)
array = []
r = open('candidates.txt','r+')
for line in r.readlines():
	p = line.split(',')
	array.append(p[0].split())

#Folding the pusar candidate 
i = 0 
for i in range(1,len(array)):
	temp = array[i][0].split(':')   #splitting the zeroth element of the array containing candidates info to obtain filename and candidate number 
	os.system("prepfold -noxwin -accelcand %d -accelfile %s.cand %s_DM%3.2f.dat > prepfold_out"%(int(temp[1]),temp[0],key2,float(array[i][1])))

i=0
for i in range(1,len(array)):
	temp = array[i][0].split(':')   #splitting the zeroth element of the array containing candidates info to obtain filename and candidate number 
	bf = open("%s_DM%3.2f_ACCEL_Cand_%d.pfd.bestprof"%(key2,float(array[i][1]),int(temp[1])),'rw')	
	p_bary = float(input_params(bf,'P_bary',4))*10**(-3)  #period to fold raw data with given candidate's period from .bestprof file 
	bf = open("%s_DM%3.2f_ACCEL_Cand_%d.pfd.bestprof"%(key2,float(array[i][1]),int(temp[1])),'rw')	
	red_chi_sq = float(input_params(bf,'Reduced chi-sqr',4))
	if red_chi_sq >= 3.00:
	#now as the number of channels are 2048 everywher taking nsub = 64 and number of bins = 64 for all candidates doing the prepfold for whole data for all given candidates. 
#		print "reduced chi-square is %1.3f"%red_chi_sq #"reached here",float(array[i][1])
		print "file is %s_DM%3.2f.dat"%(key2,float(array[i][1]))#"reached here",float(array[i][1])
		os.system("prepfold -noxwin -n 64 -nsub 64 -p %0.17f -dm %3.2f %s_DM%3.2f.dat > prepfold_rawdata_out"%(p_bary,float(array[i][1]),key2,float(array[i][1]))) 


os.chdir("%s"%results_dir)
os.system("ls *.ps | xargs -n 1 ps2pdf")
os.system("mv *ms* prepfold_raw_data")#%(key2,float(array[i][1]),int(temp[1]))) 
os.system(" mv *.pfd.pdf* prepfold_results/")
os.system("ls prepfold_raw_data/*.ps | xargs -n 1 ps2pdf")

#Running single pulse search 
os.system("single_pulse_search.py %s/*%s*.*dat > %s_single_pulse_search_out"%(results_dir,key,key))

"""

##Script to generate ephimeridis files
f = open(f_ephmer,"rw")
mjd_current = float(input_params(f,'TASC',1))
# Need to find  way to reduce this script but limited by Pyhton's way of writing files

for p in range(-N*10,N*10,int(dt*10)):  #multiplication by 10^(-(power of 10 in the dt)) as the for loop does not accept the decimal values for stepsizes
     ddt = dt*1.1574074074074073e-05*p #stepsize * sec(in MJD units) * step used
     mjd = mjd_current +ddt
     fin = open(f_ephmer,"rw")
     fout = open("%s/SAXJ1808.4-3658_TASC_%d.par"%(dir_ephmer,p),"wa+")
     for line in fin:
        fout.write(line.replace('TASC            %f'%mjd_current,'TASC            %f'%mjd))
     fout.close()
     fin.close()
     fin = open("%s/SAXJ1808.4-3658_TASC_%d.par"%(dir_ephmer,p),"rw")
     fout = open("%s/SAXJ1808.4-3658_PEPOCH_%d.par"%(dir_ephmer,p),"wa+")
     for line in fin:
        fout.write(line.replace('PEPOCH          %f'%mjd_current,'PEPOCH          %f'%mjd))
     fout.close()
     fin.close()
     fin = open("%s/SAXJ1808.4-3658_PEPOCH_%d.par"%(dir_ephmer,p),"rw")
     fout = open("%s/SAXJ1808.4-3658_%d.par"%(dir_ephmer,p),"wa+")
     for line in fin:
        fout.write(line.replace('POSEPOCH        %f'%mjd_current,'POSEPOCH        %f'%mjd))
     fout.close()
     fin.close()

os.chdir(dir_ephmer)
os.system("rm *TASC*")
os.system("rm *PEPOCH*")


"""
for p in range(-N*10,N*10,int(dt*10)):  #multiplication by 10^(-(power of 10 in the dt)) as the for loop does not accept the decimal values for stepsizes 
	ddt = dt*1.1574074074074073e-05*p #stepsize * sec(in MJD units) * step used 
	mjd = mjd_current +ddt
	fin = open("M28A.par")
	fout = open("M28A_%d.par"%p,"wt")
	for line in fin:
		fout.write(line.replace('Tasc','Pasc'))
	fin.close()
	fout.close()
"""