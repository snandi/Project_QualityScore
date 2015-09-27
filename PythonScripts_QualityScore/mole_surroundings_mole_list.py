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
import lib_get_mole_surrounding as gms 		# model3: get the surrounding pixels(donut or sandwich) around molecule backbone


###################################################################################
## FUNCTION: get the intensity of surrounding pixels of one molecule region aligned to certain interval
##           in all frames it stepping across

## $python get_mole_surroundings_mole.py IntervalNum GroupNum MoleculeNum


################################################################################
## global variables

Input_moleListFolder = "/exports/aspen/cwu/Project_QualityScore/Data/molecule_file"
Input_MoleParentFolder = "/exports/aspen/cwu/Project_QualityScore/Data/molecule_file/5pixel_region"
Output_moleParentFolder = "/exports/home/cwu269/cwu/Project_QualityScore/Data/MF" 
Output_plotParentFolder = "/exports/home/cwu269/cwu/Project_QualityScore/Plots/MF" 


ncol = 1392
nrow = 1024

#####################################################################################
## main: get the surrounding intensity of the molecules in a list
def main(argv):

    	try:                                
        	opts, args = getopt.getopt(argv, "hs:e:", ["help", "grammar="])
    	except getopt.GetoptError:          
        	usage()                         
        	sys.exit(2)   
   
	## arguments
	IntervalNum = argv[0]
	IntervalName = "refFrag_" + str(IntervalNum)
#	ChrNum = int(argv[1])
	ChrNum = 1

	## set the input and output folders for one interval
	Input_MoleParentFolder_interval = os.path.join(Input_MoleParentFolder, IntervalName)
	Output_moleParentFolder_interval = os.path.join(Output_moleParentFolder, IntervalName)
	Output_plotParentFolder_interval = os.path.join(Output_plotParentFolder, IntervalName)
	
	## read in the molecule list
	ListFileName = 'Chr' + str(ChrNum) + '_interval' + str(IntervalNum) + '.txt'
	FileReadPath = os.path.join(Input_moleListFolder, ListFileName)
	fileid = open(FileReadPath, 'r')
	mole_list = fileid.read()
	mole_list = mole_list.split(' ')

	## get the data file and image of surrounding pixels for each molecule
	for mole in mole_list:
		GroupNum = int(mole[0:7])
		MoleculeNum = int(mole[15:19])
		
		print GroupNum, MoleculeNum

		
		GroupFolder = mfc.get_grouppath(GroupNum, Input_MoleParentFolder_interval)
		if os.path.exists(GroupFolder):
			fileReadPath = mfc.get_molecule_filepath(GroupFolder, MoleculeNum)
			if os.path.exists(fileReadPath):
				gms.surrounding_OneMolecule_allframe(Input_MoleParentFolder_interval, Output_moleParentFolder_interval, Output_plotParentFolder_interval, GroupNum, MoleculeNum, 2, False)
			else: 
				print mole, "not exist!"
		else:
			print mole, "not exist!"

	print "end !"	

    	
	
	


#####################################################################################
## calling the main function

if __name__ == "__main__":
	main(sys.argv[1:])










	
		












		
