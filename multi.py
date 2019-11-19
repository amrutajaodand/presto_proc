import os, sys, time
import numpy as np
import multiprocessing as mp


""" Function to parallely find Fitting Factor using several starting points"""
def parallel_process(n):

	# The 'nth' dataing value 	
	m = datapoints[n]				# data value of m		
	

	max_match_array = np.zeros(5)				# initialize array to store return values 
	max_match_array[0]= m			
	max_match_array[1] = m**2
	max_match_array[2] = m**3
	max_match_array[3] = m**4
	max_match_array[4] = m**5
	
	time.sleep(1)		#FIXME this is just for testing, be sure to delete this afterwards


	return max_match_array
start_time = time.time()  						# time at which maximization starts

datapoints = np.array([1,2,3,4,5,6])	# can also be a multidimensional array
N_cores = 6	# number of cores to be used at once (to test the script vary this from 1,2,3 etc)
if N_cores > 8:	# safety check
	print 'Error: Asking for too many cores, reduce N_cores and try again'
	exit(-1)
N_datapoints = len(datapoints) 
index = range(N_datapoints)                                            # input for parallel_process
pool = mp.Pool(processes=N_cores)                                       # specifies number of parallel processes
match_results = pool.map(parallel_process,index)        # stores the output of all parallel processes
pool.close()         
end_time = time.time()							# time at which maximization ends

print match_results                                                                   # close all worker processes
print 'time taken =%f s'%(end_time-start_time)
