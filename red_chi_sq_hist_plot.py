"This scripts is written to plot the histograms for each 15 min integration absed on the reduced chi square" 
################################################################################################
from glob import glob 
from matplotlib import pyplot as plt 
import numpy as np 
import os, re
################################################################################################
filelist = glob("/projects/0/lotaas2/data/jaodand/results/multicore/*201408*/*/*data*.txt")

plt.ioff()

	

for f in filelist:
     key = re.search('data_(.+?).txt', f).group(1)
     plt.figure()
     labels = np.loadtxt(f,delimiter='\t',usecols=[0],dtype = np.str)
     values = np.loadtxt(f,delimiter='\t',usecols=[1,2,3])
 
     ##plotting the 2013 ephmeris cands with both p and pdot search 
     idx_1 = np.where(labels == 'prepfold_ephmer_%s'%key)
     val_1 = values[idx_1]	
     counts_1 , bin_edges_1 = np.histogram(val_1[:,-1], bins = 200)
     plt.semilogy(bin_edges_1[:-1],counts_1,color='cyan',label='prepfold_ephmer')

     ##plotting the 2013 ephmeris cands with both no p and but pdot search 
     idx_2 = np.where(labels == 'prepfold_ephmer_nopsearch_%s'%key)
     val_2 = values[idx_2]
     counts_2 , bin_edges_2 = np.histogram(val_2[:,-1], bins = 200)
     plt.semilogy(bin_edges_2[:-1],counts_2,color='blue',label='prepfold_ephmer_nopsearch')


     ##plotting the 2015 ephmeris cands with both p and pdot search 
     idx_3 = np.where(labels == 'prepfold_ephmer_2015_%s'%key)
     val_3 = values[idx_3]
     counts_3 , bin_edges_3 = np.histogram(val_3[:,-1], bins = 200)
     plt.semilogy(bin_edges_3[:-1],counts_3,color='red',label='prepfold_ephmer_2015')

     ##plotting the 2015 ephmeris cands with both no p and but pdot search 
     idx_4 = np.where(labels == 'prepfold_ephmer_2015_nopsearch_%s'%key)
     val_4 = values[idx_4]
     counts_4 , bin_edges_4 = np.histogram(val_4[:,-1], bins = 200)
     plt.semilogy(bin_edges_4[:-1],counts_4,color='green',label='prepfold_ephmer_2015_nopsearch')

     counts, bin_edges = np.histogram(values[:,-1], bins = 200)
     plt.semilogy(bin_edges[:-1],counts,color='black')
     plt.legend(loc=1,fontsize=8.0)
     plt.xlabel("Red_chi_sq")
     print (f.split("/"))[-1]
     plt.savefig("%s.pdf"%((f.split("/"))[-1][:-4]))

     #plt.hist(values[:-1])
