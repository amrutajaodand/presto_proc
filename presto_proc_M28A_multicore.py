"""This script can be used to the presto proc multicore script on M28A: testing
A.Jaodand, 20/1/2015
Last modification:15 Jan 2015
"""
#################################################################################
#Deciding the input parameters 

highDM = 1000.0 # Highest DM value
nsub = 32	#number of subbands 
r = 0.1		#acceptable resolution in time
N = 30          # Integer value which will serve as upper and lower limit for Tasc range 
dt = 0.1        #Step size for varying Tasc 
Nf = 1 		#number of .fits files being used here
Ncore_max = 24   #Upper limit on core usage 
#################################################################################
#importing modules
import numpy as np
import os, sys, getopt, re, shutil, math, time
from datetime import datetime
import multiprocessing as mp
import ddplan_parallelizer as ddp_par
import glob
#################################################################################
##  input parameters 
### use the arguments given by the user ##########################################
args=sys.argv[1:]   # list of arguments
try:
   opts, args = getopt.getopt(args,"h",["folder_tag=","file_tag=","N_cores="])#, "m=","incl_angle="]) #folder_tag name of the folder
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
#################################################################################
# import the files to be processed with PRESTO
path = "/projects/0/lotaas/data/jaodand/GBT14B-036_AMXPs/%s"%folder_tag #remember to change the folder number ot a tag passed by the run script

infilename = "%s"%file_tag
filename,fileext = os.path.splitext(infilename)

os.system("ls %s"%os.path.join(path,"%s.fits"%(filename)))#[:-1],'[1-3]')))
infilepath = os.path.join(path,"%s.fits"%(filename))#[:-1],'[1-3]'))
print infilepath	

# folder where the script is located 
src_dir = "/projects/0/lotaas/data/jaodand/scripts/presto_proc_15min"

##Initialise the files to store the results of the processing 
#create folder for a given date
outpath= "/projects/0/lotaas/data/jaodand/results/multicore"
if not os.path.exists("%s"%outpath):
        os.mkdir("%s"%outpath)

outdir = os.path.join(outpath,"%s"%folder_tag)   
if not os.path.exists("%s"%outdir):
	os.mkdir("%s"%outdir)

key =filename[:-5] #used for naming all the files removes last number from the observation filename like 000[1-12]
key2 = "%s_out"%filename[:-5] #used for naming the .dat files from the prepsubband - important key as we need this for file movements later  
print key2
 

#create folder specific to the particular data on a given date(results_dir)
results_dir = os.path.join(outpath,"%s"%folder_tag,key) 
if not os.path.exists("%s"%results_dir):
	os.mkdir("%s"%results_dir)

#create folder to store final raw data, timeseries and ephimerids based folding results

if not os.path.exists(os.path.join(results_dir,"prepfold_raw_data_%s"%key)):
	os.mkdir(os.path.join(results_dir,"prepfold_raw_data_%s"%key))
if not os.path.exists(os.path.join(results_dir,"prepfold_timeseries_%s"%key)):
	os.mkdir(os.path.join(results_dir,"prepfold_timeseries_%s"%key))
if not os.path.exists(os.path.join(results_dir,"prepfold_ephmer_%s"%key)):
	os.mkdir(os.path.join(results_dir,"prepfold_ephmer_%s"%key))
	

#file containing the latest ephimeridis and initiating folder to store ephimeridis generated
#dir_ephmer = "ephimeridis_%s"%(filename)
#if not os.path.exists(os.path.join(path,"%s"%dir_ephmer)):
#        os.mkdir(os.path.join(path,"%s"%dir_ephmer))
#f_ephmer="%s/%s.par"%(path,infilename[12:-15]) #infilename [12:-15] gives the name of the source from the infilename
#################################################################################
##define required functions
def roundup(x):
#	print x
	return (int(math.ceil(float(x)/(10**(len(str(x))-2)))))*(10**(len(str(x))-2))

def input_params(f,y,g):     
    for line in f:
             if y in line:
#	             print line
                     temp = line.split();
                    # print temp
    	     	     return temp[g]  
#################################################################################
####Start the processing with PRESTO 

##readfile and store metdata
os.system("readfile %s > %s"%(os.path.join(path,infilename),os.path.join(results_dir,"metdata_%s"%filename)))

f = open(os.path.join(results_dir,"metdata_%s"%filename),'rw')
numout = Nf*int(input_params(f,'Spectra per file',4))
print "numout is %d"%numout
"""
##find rfi and move the o/p files to results_dir
print "rfifind -time 2.0 -o %s %s"%(key,infilepath)
os.system("rfifind -time 2.0 -o %s %s"%(key,infilepath))

print "mv %s/*%s* %s"%(src_dir,key,results_dir)
os.system("mv %s/*%s* %s"%(src_dir,key,results_dir))

##looking for persistent low level rfi/prepdata
#First obtaining the numout from metdata file
f = open(os.path.join(results_dir,"metdata_%s"%filename),'rw')
numout =Nf*int(input_params(f,'Spectra per file',4))
print numout 

#now running the prepadata command
print "the loc directory is %s"%(results_dir)
print "prepdata -nobary -o %s_topo_DM0.00 -dm 0.0 -mask %s/%s_rfifind.mask -numout %d %s > prepdata_out_%s"%(key,results_dir,key,numout,infilepath,key)
os.system("prepdata -nobary -o %s_topo_DM0.00 -dm 0.0 -mask %s/%s_rfifind.mask -numout %d %s > prepdata_out_%s"%(key,results_dir,key,numout,infilepath,key))

#exploring the .dat file and obtaining the fft
#os.system("exploredat %s_topo_DM0.00.dat"%(key))
print "realfft %s_topo_DM0.00.dat"%(key)
os.system("realfft %s_topo_DM0.00.dat"%(key))

#Running the accelsearch 
print "accelsearch -numharm 4 -zmax 100 %s_topo_DM0.00.fft > accelsearch_out_%s"%(key,key)
os.system("accelsearch -numharm 4 -zmax 100 %s_topo_DM0.00.fft > accelsearch_out_%s"%(key,key))

###Don't forget to add the .birds steps here!!!

#DDplan
##Obtaining the input arguments for DDplan.py 

print "open(os.path.join(results_dir,metdata_filename)),'rw')"
f = open(os.path.join(results_dir,"metdata_%s"%filename),'rw')
nchan = int(input_params(f,'Number of channels',4))	#number of channels
print "number of channels is %d"%nchan

f = open(os.path.join(results_dir,"metdata_%s"%filename),'rw')
dbw = float(input_params(f,'Total Bandwidth',4))	   	#Channel width 
print "bandwidth is %f"%dbw 

f = open(os.path.join(results_dir,"metdata_%s"%filename),'rw')
dt = float(input_params(f,'Sample time',4))
dt = dt*10**(-6)					#Sample time 
print "sample time is %f"%dt

f = open(os.path.join(results_dir,"metdata_%s"%filename),'rw')
f0 =  float(input_params(f,'Central freq',4))		#Central frequency 
print "central frequency is %f"%f0

##Running DDplan.py

print "DDplan.py -d %d -n %d -b %f -t %f -f %f -s %d -r %f -o %s_ddplan> %s/ddplan_out_%s"%(highDM,nchan,dbw,dt,f0,nsub,r,results_dir,results_dir,key)
os.system("DDplan.py -d %d -n %d -b %f -t %f -f %f -s %d -r %f -o %s_ddplan> %s/ddplan_out_%s"%(highDM,nchan,dbw,dt,f0,nsub,r,results_dir,results_dir,key))
g = open(os.path.join(results_dir,"ddplan_out_%s"%key),'rw')#out is to store the output from ddplan 
contents1, sentinel, contents2 = g.read().partition("Low DM    High DM     dDM  DownSamp  dsubDM   #DMs  DMs/call  calls  WorkFract\n") #splitting the o/p of DDplan to obtain the information about different runs of prepsubbabnd

##saving DDplan step to be used for prepsubband to a file
h = open(os.path.join(results_dir,"ddplan_%s"%key),'w+') #without out in the name this file contains the ddplan to be followed in the prepsubband 
h.write(sentinel) #writing the header of the file -- Low DM    High DM     dDM  DownSamp  dsubDM   #DMs  DMs/call  calls  WorkFract
h.write(contents2) # ddplan steps for subbands 
h.read()
h.close()

#creating a distributed ddplan for a given number of cores 
ddp_par.ddplan_par(N_cores,results_dir,"ddplan_%s"%key)

#Prepsubband

##Executing  prepsubband for every subband for corresponding ddm and number of calls

h = open(os.path.join(results_dir,"ddplan_%s_parallel"%key),'rw')
num_lines = sum(1 for line in h) #number of subbands

i = 0
datapoints = []
h = open(os.path.join(results_dir,"ddplan_%s_parallel"%key),'rw')
for i, line in enumerate(h,1):        #i counter for no. of subbands,  
        if i in range(2,num_lines+1):   # start from 2 as the line number 1 is for header end at num_lines-3 as last three are 4 empty lines at end. 
                temp = line.split()
		print temp
                datapoints.append(temp)

print datapoints
print len(datapoints)
for nc in range(1,N_cores+1):	
	if not os.path.exists(os.path.join(src_dir,"%s_core_%d"%(key2,nc))):
		os.mkdir(os.path.join(src_dir,"%s_core_%d"%(key2,nc)))
        
#############################################################
#function to parallelize the prepsubband 
def prepsubband_parallel_proc(n):
        temp = datapoints[n]
        print temp[7]
	for j in range(0,int(temp[7])):
                lodm = float(temp[0])+(j*float(temp[2])*float(temp[6]))
                print lodm
		logfile = "log_%s_%d.log"%(key,int(temp[9]))#logfile to store prepsubband log for a given core. 
		temp_folder = "%s_core_%d"%(key2,int(temp[9]))#folder used  to store o/p for a given core.
                print "prepsubband -nsub %d -lodm %f -dmstep %f -numdms %d -numout %d -downsamp %d -mask %s/%s_rfifind.mask -o %s %s"%(nsub,lodm,float(temp[2]),float(temp[6]),numout/int(temp[3]),int(temp[3]),results_dir,key,key2,infilepath)
		os.chdir("%s"%temp_folder)
                os.system("prepsubband -nsub %d -lodm %f -dmstep %f -numdms %d -numout %d -downsamp %d -mask %s/%s_rfifind.mask -o %s %s"%(nsub,lodm,float(temp[2]),float(temp[6]),numout/int(temp[3]),int(temp[3]),results_dir,key,key2,infilepath))
		#os.system("prepsubband -nsub %d -lodm %f -dmstep %f -numdms %d -numout %d -downsamp %d -mask %s/%s_rfifind.mask -o %s %s >> %s"%(nsub,lodm,float(temp[3]),float(temp[6]),numout/int(temp[3]),int(temp[3]),results_dir,key,key2,infilepath,logfile)) #numout/downsampling in arguments !! 
		os.chdir("%s"%src_dir)
#############################################################
#Passing arguments and parallelized prepsubband
start_time = time.time()
#
if N_cores > Ncore_max: # safety check
        print 'Error: Asking for too many cores, reduce N_cores and try again'
        exit(-1)
N_datapoints = len(datapoints)
index = range(N_datapoints)                                            # input for parallel_process
print index
print "index was given above"
pool = mp.Pool(processes=N_cores)                                       # specifies number of parallel processes
pool.map(prepsubband_parallel_proc,index)                               # index sequence passed to pool for reference to datapoint at that index 
pool.close()
end_time = time.time()                                                  # time at which maximization ends
print 'time taken =%f s'%(end_time-start_time)
#############################################################

os.system("mv *%s_core_*/*%s* %s"%(key2,key2,results_dir))
os.system("rm -r *%s_core_*"%(key2))

#mv all the recently generated .dat files to the results_dir 
print "mv *%s* %s"%(key,results_dir)
os.system("mv *%s* %s"%(key,results_dir))
os.system("mv *%s* %s"%(key2,results_dir))

#Prepare to search the data
##Use xargs to fft and zap the data 
print "ls %s/*.dat | xargs -n 1 realfft > %s/realfft_out_%s"%(results_dir,results_dir,key)
os.system("ls %s/*.dat | xargs -n 1 realfft > %s/realfft_out_%s"%(results_dir,results_dir,key))

#!!!Add zapping commands later 

#Searching for periodic signals 
##First making an array with names of all files over which accelsearch will be run (purpose:multiprocessing)
os.system("ls %s/*out*.fft > %s/accelsearch_input_%s"%(results_dir,results_dir,key))#saving list of all the .fft files.
j = 0
accel_datapoints = []
h_accel = open("%s/accelsearch_input_%s"%(results_dir,key),'rw')
for line in h_accel:        #i counter for no. of subbands,  
        temp = line.split()
        accel_datapoints.append(temp[0])
print len(accel_datapoints)
#############################################################
##function to parallelise the accelsearch 
def accel_parallel_proc(n):
        temp = accel_datapoints[n]

        print "accelsearch -zmax 100 %s > %s/accelsearch_log_%d"%(temp,results_dir,n)
        os.system("accelsearch -zmax 100 %s > %s/accelsearch_log_%d"%(temp,results_dir,n))

#        time.sleep(1)           #FIXME this is just for testing, be sure to delete this afterwards
#############################################################
start_time = time.time()
#
if N_cores > Ncore_max: # safety check
        print 'Error: Asking for too many cores, reduce N_cores and try again'
        exit(-1)
N_accel_datapoints = len(accel_datapoints)
index = range(N_accel_datapoints)                                            # input for parallel_process
print index
pool = mp.Pool(processes=N_cores)                                       # specifies number of parallel processes
pool.map(accel_parallel_proc,index)                                     # index sequence passed to pool for reference to datapoint at that index 
pool.close()
end_time = time.time()                                                  # time at which maximization ends
print 'time taken =%f s'%(end_time-start_time)

#concatenate all the accelsearch log files 
files = sorted(glob.glob('%s/*accelsearch_log*'%results_dir))
with open('accelsearch_data_out.txt','w') as result:
        for file_ in files:
                for line in open( file_, 'r'):
                        result.write(line)
os.system("rm -rf %s/*accelsearch_log*"%results_dir)
#############################################################

#sifting periodic candidates 
shutil.copy('ACCEL_sift.py',results_dir)
os.chdir("%s"%results_dir)
print "python ACCEL_sift.py"
results = os.popen("python ACCEL_sift.py").readlines()
res = ''.join(results)
h = open("cands_%s.txt"%key,'w+')
h.write(res)
h.close()

# we need to save all the accelcands from the cands file generated above for next prepfold steps
g = open('cands_%s.txt'%key,'rw') 
contents1, sentinel, contents2 = g.read().partition("                           file:candnum                               DM     SNR    sigma   numharm    ipow     cpow      P(ms)          r         z     numhits \n")

h = open("cands_filt_%s.txt"%key,'w+') 
h.write(sentinel)
h.write(contents2)
h.close()
f = open("cands_filt_%s.txt"%key,'r')      #the part of the cands file containing the candidates are stored 
m = open("candidates_%s.txt"%key,'wa+')      #from the cands_filt file only the candidates which will be used for prepfold are stored here 
m.write(sentinel)

for line in f:
	if key2 in line:
		 m.write(line)
m.close()
"""
os.chdir("%s"%results_dir)
array = []
r = open('candidates_%s.txt'%key,'r+')
for line in r.readlines():
	p = line.split(',')
	array.append(p[0].split())

print len(array)
"""
#Folding the pusar candidate 
for i in range(1,len(array)):
	temp = array[i][0].split(':')   #splitting the zeroth element of the array containing candidates info to obtain filename and candidate number 
	print "prepfold -noxwin -accelcand %d -accelfile %s.cand %s_DM%3.2f.dat > prepfold_out"%(int(temp[1]),temp[0],key2,float(array[i][1]))
	os.system("prepfold -noxwin -accelcand %d -accelfile %s.cand %s_DM%3.2f.dat > prepfold_out"%(int(temp[1]),temp[0],key2,float(array[i][1])))

os.system("ls *.ps | xargs -n 1 ps2pdf")
os.system(" mv *.pfd.pdf* prepfold_timeseries_%s/"%key)
os.system("rm *.pfd")
os.system("rm *out**.ps")
"""
############################################################################################################################################################
# the sections for single_pulse_search, prepfold a/c ephimeridis nd prepfold of the .fits files for candidates through accelsearch are below not being used currently.


i=0
for i in range(1,len(array)):
	temp = array[i][0].split(':')   #splitting the zeroth element of the array containing candidates info to obtain filename and candidate number 
	bf = open("%s_DM%3.2f_ACCEL_Cand_%d.pfd.bestprof"%(key2,float(array[i][1]),int(temp[1])),'rw')	
	p_bary = float(input_params(bf,'P_bary',4))*10**(-3)  #period to fold raw data with given candidate's period from .bestprof file 
	bf = open("%s_DM%3.2f_ACCEL_Cand_%d.pfd.bestprof"%(key2,float(array[i][1]),int(temp[1])),'rw')	
	red_chi_sq = float(input_params(bf,'Reduced chi-sqr',4))
	if red_chi_sq >= 3.00:
	#now as the number of channels are 2048 everywher taking nsub = 64 and number of bins = 64 for all candidates doing the prepfold for whole data for all given candidates. 
		print "reduced chi-square is %1.3f"%red_chi_sq #"reached here",float(array[i][1])
		os.chdir("%s"%path)
		os.system("prepfold -noxwin -n 64 -nsub 64 -p %f -dm %3.2f %s > prepfold_raw_data_out"%(p_bary,float(array[i][1]),infilename))
		os.chdir("%s"%results_dir)
#
os.chdir("%s"%results_dir)
os.system("ls *.ps | xargs -n 1 ps2pdf")
os.system("mv %s/*ms* prepfold_raw_data_%s"%(path,key))#,float(array[i][1]),int(temp[1]))) 
os.system("ls prepfold_raw_data_%s/*.ps | xargs -n 1 ps2pdf"%key) 
"""
os.chdir("%s"%results_dir)
##Script to generate ephimeridis files
p = 0
f = open(f_ephmer,"rw")
mjd_current = float(input_params(f,'TASC',1))

# Need to find  way to reduce this script but limited by Pyhton's way of writing files

############producing and saving all the ephimeridis files. 
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
os.system("rm %s/*TASC*"%dir_ephmer)
os.system("rm %s/*PEPOCH*"%dir_ephmer)

N = 0.1
##################using the ephimeridis files generated above to fold every .dat file.
i = 0 
ing the ephimeridis files generated above to fold every .dat file.
i = 0 
for i in range(1,2):#len(array)):       #counter over the DM values as given by the accel_sift 
        if not os.path.exists(os.path.join(results_dir,"prepfold_ephmer_%s/prepfold_ephmer_%s_DM%3.2f"%(filename,filename,float(array[i][1])))):
                os.mkdir(os.path.join(results_dir,"prepfold_ephmer_%s/prepfold_ephmer_%s_DM%3.2f"%(filename,filename,float(array[i][1])))) #folder to store every .dat (DM)'s results for folding it using all ephimeridis files  

        os.chdir(os.path.join(results_dir,"prepfold_ephmer_%s/prepfold_ephmer_%s_DM%3.2f"%(filename,filename,float(array[i][1])))) #changing path to above initialised folder 

        #Now folding .dat file with every ephimeridis file
        for p in range(int(-N*10),int(N*10),int(dt*10)):  #multiplication by 10^(-(power of 10 in the dt)), for loop does not accept the decimal values for stepsizes
                os.system("prepfold -noxwin -par %s/ephimeridis/SAXJ1808.4-3658_%d.par -fine -mask ../../%s_rfifind.mask ../../%s_DM%3.2f.dat"%(path,p,key,key2,float(array[i][1])))
        #       print "prepfold -noxwin -par %s/ephimeridis/SAXJ1808.4-3658_%d.par -fine -mask ../../%s_rfifind.mask ../../%s_DM%3.2f.dat > prepfold_ephmer_out"%(path,p,key,key2,float(array[i][1]))
        os.chdir("%s"%results_dir)
                
#               os.system("ls | sed 'p;s/PSR/%d_PSR/'| xargs -n2 mv"%p)         
#               os.system("mv *%s_DM%3.2f_%d_PSR* prepfold_ephmer_%s/prepfold_ephmer_%s_DM%3.2f"%(key2,float(array[i][1]),p,key,key,float(array[i][1])))
#       os.chdir(os.path.join(results_dir,"prepfold_ephmer_%s"%filename))       #changing to the directory for keeping the zipped format of folder with all the .ps files 
#       os.system("tar -zcvf prepfold_ephmer_%s_DM%3.2f.tar.gz prepfold_ephmer_%s/prepfold_ephmer_%s_DM%3.2f"%(filename,float(array[i][1]),filename,filename,float(array[i][1])))   

os.system("ls prepfold_ephmer_%s/*.ps | xargs -n 1 ps2pdf"%key) 


#########singlepulse_search.py is being run here
os.chdir("%s"%results_dir)
#Running single pulse search 
os.system("/bin/ls %s*.*dat | xargs -n 1 single_pulse_search.py -p> %s_single_pulse_search_out"%(key2,key2))
#os.system("cat *singlepulse > %s_singlepulse_out"%key)
os.system("single_pulse_search.py *.singlepulse")
"""
