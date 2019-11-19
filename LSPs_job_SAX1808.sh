#!/bin/bash -x
#SBATCH -p normal
#SBATCH -t 80:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=A.D.Jaodand@uva.nl
##########################################################
#
# LOTAAS Single Pulse search job launcher
#
# written by Daniele Michilli
#modified :A.D.Jaodand
#
# Not written for direct use, use LSPs_launcher.sh instead
#
##########################################################
#python presto_proc_M28A_multicore.py --folder_tag 20140721 --file_tag 'guppi_56859_M28A_0001_0001.fits' --N_cores 16 > test.out 

#python presto_proc_multicore.py --folder_tag 20140809 --file_tag 'guppi_56878_SAXJ1808.4-3658_0001_0001.fits' --N_cores 24 --f_num [1-3]> guppi_56878_SAX1808_1_3.out
python presto_proc_multicore_2015.py --folder_tag 20140809 --file_tag 'guppi_56878_SAXJ1808.4-3658_0001_0001.fits' --N_cores 24 --f_num [1-3]> guppi_56878_SAX1808_1_3_2015.out

#python presto_proc_multicore.py --folder_tag 20140809 --file_tag 'guppi_56878_SAXJ1808.4-3658_0001_0004.fits' --N_cores 24 --f_num [4-6]> guppi_56878_SAX1808_4_6.out
#python presto_proc_multicore_2015.py --folder_tag 20140809 --file_tag 'guppi_56878_SAXJ1808.4-3658_0001_0004.fits' --N_cores 24 --f_num [4-6]> guppi_56878_SAX1808_10_12_2015.out

#python presto_proc_multicore.py --folder_tag 20140809 --file_tag 'guppi_56878_SAXJ1808.4-3658_0001_0007.fits' --N_cores 24 --f_num [7-9]> guppi_56878_SAX1808_7_9.out
#python presto_proc_multicore_2015.py --folder_tag 20140809 --file_tag 'guppi_56878_SAXJ1808.4-3658_0001_0007.fits' --N_cores 24 --f_num [7-9]> guppi_56878_SAX1808_7_9_2015.out

#python presto_proc_multicore.py --folder_tag 20140809 --file_tag 'guppi_56878_SAXJ1808.4-3658_0001_0010.fits' --N_cores 24 --f_num [0-2]> guppi_56878_SAX1808_10_12.out
#python presto_proc_multicore_2015.py --folder_tag 20140809 --file_tag 'guppi_56878_SAXJ1808.4-3658_0001_0010.fits' --N_cores 24 --f_num [0-2]> guppi_56878_SAX1808_10_12_2015.out
#python presto_proc_nopsearch_multicore.py --folder_tag 20140809 --file_tag 'guppi_56878_SAXJ1808.4-3658_0001_0010.fits' --N_cores 24 --f_num [0-2]> guppi_56878_SAX1808_10_12_nopsearch.out
#python presto_proc_nopsearch_multicore_2015.py --folder_tag 20140809 --file_tag 'guppi_56878_SAXJ1808.4-3658_0001_0010.fits' --N_cores 24 --f_num [0-2]> guppi_56878_SAX1808_10_12_nopsearch_2015.out

#python presto_proc_multicore.py --folder_tag 20140822 --file_tag 'guppi_56891_SAXJ1808.4-3658_0006_0004.fits' --N_cores 24 --f_num [4-6]> guppi_56891_SAX1808_4_6.out
#python presto_proc_multicore_2015.py --folder_tag 20140822 --file_tag 'guppi_56891_SAXJ1808.4-3658_0006_0004.fits' --N_cores 24 --f_num [4-6]> guppi_56891_SAX1808_4_6_2015.out
#python presto_proc_nopsearch_multicore.py --folder_tag 20140822 --file_tag 'guppi_56891_SAXJ1808.4-3658_0006_0004.fits' --N_cores 24 --f_num [4-6]> guppi_56891_SAX1808_4_6_nopsearch.out 
#python presto_proc_nopsearch_multicore_2015.py --folder_tag 20140822 --file_tag 'guppi_56891_SAXJ1808.4-3658_0006_0004.fits' --N_cores 24 --f_num [4-6]> guppi_56891_SAX1808_4_6_nopsearch_2015.out 

#python presto_proc_multicore.py --folder_tag 20140822 --file_tag 'guppi_56891_SAXJ1808.4-3658_0006_0001.fits' --N_cores 24 --f_num [1-3]> guppi_56891_SAX1808_1_3.out
#python presto_proc_multicore_2015.py --folder_tag 20140822 --file_tag 'guppi_56891_SAXJ1808.4-3658_0006_0001.fits' --N_cores 24 --f_num [1-3]> guppi_56891_SAX1808_1_3_2015.out
#python presto_proc_nopsearch_multicore.py --folder_tag 20140822 --file_tag 'guppi_56891_SAXJ1808.4-3658_0006_0001.fits' --N_cores 24 --f_num [1-3]> guppi_56891_SAX1808_1_3_nopsearch.out 
#python presto_proc_nopsearch_multicore_2015.py --folder_tag 20140822 --file_tag 'guppi_56891_SAXJ1808.4-3658_0006_0001.fits' --N_cores 24 --f_num [1-3]> guppi_56891_SAX1808_1_3_nopsearch_2015.out 



exit 0
