#!/bin/bash -x
#SBATCH -p normal
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=A.D.Jaodand@uva.nl
##########################################################
#
# LOTAAS Single Pulse search job launcher
#
# written by Daniele Michilli
#modiefied :A.D.Jaodand
#
# Not written for direct use, use LSPs_launcher.sh instead
#
##########################################################
#python presto_proc_M28A_multicore.py --folder_tag 20140721 --file_tag 'guppi_56859_M28A_0001_0001.fits' --N_cores 16 > test.out 
python presto_proc_M28A_multicore.py --folder_tag 20140809 --file_tag 'guppi_56878_M28A_0002_0001.fits' --N_cores 16 > test.out 
exit 0
