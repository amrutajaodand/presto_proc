"""This script calculates the phase modulation with orbit due to Tasc
Author : Amruta 
Created: 30/01/2015

"""
#################################################################################
#importing modules
import numpy as np
import os, sys, getopt, re, shutil, math
from datetime import datetime
#################################################################################
#################################################################################
def phs_mod_lim(M1,M2,iota,P_orb,P_spin):
	#units of input parameters 
	#P_orb(hr), P_spin(sec), M1 & M2 (M_sun),iota (rad)
	M_sun = 4.9255*10**(-6)

	M = (M1+M2)*M_sun #total mass
	f1 = (M2*M_sun*np.sin(iota))**3/(M)**2 #Mass funstion 

	Phs_mod_amp = 7.375*f1**(1./3)*P_orb**(2./3)/P_spin 
	Phs_mod_ang = 2*np.pi*((t0-Tasc)./P_orb)+np.pi/2
	return(Phs_mod_amp,Phs_mod_ang)
