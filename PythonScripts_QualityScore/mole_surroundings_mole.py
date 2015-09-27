#! /home/nandi/anaconda/bin/python

import getopt
import gzip
import io
import MySQLdb
import numpy as np
import os
import pickle
import shutil
import struct
import sys
import copy
import matplotlib.pyplot as plt
import matplotlib.image as mpimg



from string import translate, maketrans, punctuation 
from os import path

import lib_frames_index as fi  			# model1: get access to frames in database
import lib_molecule_file_control as mfc 	# model2: operation about molecule files
import lib_get_mole_surrounding as gms 

###################################################################################
## FUNCTION: get the intensity of surrounding pixels of one molecule region aligned to certain interval
##           in all frames it stepping across

## $python get_mole_surroundings_mole.py IntervalNum GroupNum MoleculeNum


################################################################################
## global variables

Input_MoleParentFolder = "/aspen/nandi/MF_cap348/maps_inca34"
Output_moleParentFolder = "/exports/home/cwu269/cwu/Project_QualityScore/Data/surrounding/" 
Output_plotParentFolder = "/exports/home/cwu269/cwu/Project_QualityScore/Plots/surrounding/" 

ncol = 1392
nrow = 1024





#####################################################################################
## main: get the surrounding intensity of one molecule
def main(argv):

    	try:                                
        	opts, args = getopt.getopt(argv, "hs:e:", ["help", "grammar="])
    	except getopt.GetoptError:          
        	usage()                         
        	sys.exit(2)   
   
	GroupNum = int(args[0])
	MoleculeNum = int(args[1])

	gms.surrounding_OneMolecule_allframe(Input_MoleParentFolder, Output_moleParentFolder, Output_plotParentFolder, GroupNum, MoleculeNum, 2, True)
	


#####################################################################################
## calling the main function

if __name__ == "__main__":
	main(sys.argv[1:])










	
		












		
