#! /home/nandi/anaconda/bin/python

import getopt
import io
import numpy as np
import os
import pickle
import shutil
import struct
import sys
import copy


from string import translate, maketrans, punctuation 
from os import path

import lib_frames_index as fi  			# model1: get access to frames in database
import lib_molecule_file_control as mfc 	# model2: operation about molecule files
import lib_get_mole_surrounding as gms		# model3: get the surrounding pixels(donut or sandwich) around molecule backbone
import lib_get_dust as gd			# model4: detect the nearby noisy pixels and the backbone pixels effected by garbage

################################################################################
## global variables

Input_MoleParentFolder = "/exports/aspen/cwu/Project_QualityScore/Data/molecule_file/1pixel_region"
Input_SurroundingParentFolder = "/exports/home/cwu269/cwu/Project_QualityScore/Data/MF" 


################################################################################
def main(argv):
	try:                                
        	opts, args = getopt.getopt(argv, "hs:e:", ["help", "grammar="])
    	except getopt.GetoptError:          
        	usage()                         
        	sys.exit(2)   

	IntervalNum = int(args[0])
	GroupNum = int(args[1])
	MoleculeNum = int(args[2])
#	layer = int(args[3])
	layer = 1

	IntervalName = "refFrag_" + str(IntervalNum)
	Input_MoleParentFolder_interval = os.path.join(Input_MoleParentFolder, IntervalName)
	Input_SurroundingParentFolder_interval = os.path.join(Input_SurroundingParentFolder, IntervalName)

	
	threshold = 5800 ## just for test

	gd.Effected_OneMolecule(Input_MoleParentFolder_interval, Input_SurroundingParentFolder_interval, GroupNum, MoleculeNum, layer, threshold)



#################################################################################
if __name__ == "__main__":
	main(sys.argv[1:])




