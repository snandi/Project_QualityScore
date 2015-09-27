#! /home/nandi/anaconda/bin/python

import getopt
import io
import numpy as np
import os
import pickle
import shutil
import struct
import sys
import math

from os import path

import lib_frames_index as fi  			# model1: get access to frames in database
import lib_molecule_file_control as mfc 	# model2: operation about molecule files
import lib_get_mole_surrounding as gms 		# model3: get the surrounding pixels(donut or sandwich) around molecule backbone
import lib_simplify_line as sl			# model5: functions related to the simplify line algorithm

###################################################################################
## FUNCTION: get the straightness measure of the molecule fragments with simplify algorithm


##############################################
## global variables

Input_moleListFolder = "/exports/aspen/cwu/Project_QualityScore/Data/molecule_file"
Input_MoleculeParentFolder = "/exports/aspen/cwu/Project_QualityScore/Data/molecule_file/5pixel_region"
Output_MoleculeParentFolder = "/exports/aspen/cwu/Project_QualityScore/Data/MF"

threshold = 3

###################################################################################
## main
def main(argv):

    	try:                                
        	opts, args = getopt.getopt(argv, "hs:e:", ["help", "grammar="])
    	except getopt.GetoptError:          
        	usage()                         
        	sys.exit(2) 

	IntervalNum = int(argv[0])
	IntervalName = "refFrag_" + str(IntervalNum)
#	ChrNum = int(argv[1])
	ChrNum = 1

	## set the input and output folders for one interval
	Input_MoleculeParentFolder_interval = os.path.join(Input_MoleculeParentFolder, IntervalName)
	Output_MoleculeParentFolder_interval = os.path.join(Output_MoleculeParentFolder, IntervalName)

	## read in the molecule list
	ListFileName = 'Chr' + str(ChrNum) + '_interval' + str(IntervalNum) + '.txt'
	FileReadPath = os.path.join(Input_moleListFolder, ListFileName)
	fileid = open(FileReadPath, 'r')
	mole_list = fileid.read()
	mole_list = mole_list.split(' ')

	## get the file of straight score
	os.chdir(Output_MoleculeParentFolder_interval)
	molecule_filename_output = "refFrag" + str(IntervalNum) + "_StraightScore.txt"
	fileid = open(molecule_filename_output, "w")
#	fileid.write(" ".join(["distance threshold = ", str(threshold), "\n\n"]))
	fileid.write(" ".join(["MoleculeID","FrameNum","theta_net","theta_max","score","\n"]))

	for mole in mole_list:
		if len(mole) == 0: continue
		GroupNum = int(mole[0:7])
		MoleculeNum = int(mole[15:19])

		GroupFolder = mfc.get_grouppath(GroupNum, Input_MoleculeParentFolder_interval)
		if os.path.exists(GroupFolder):
			fileReadPath = mfc.get_molecule_filepath(GroupFolder, MoleculeNum)
			if os.path.exists(fileReadPath):
			
				straightness_list = sl.simplify_lines_AllFrame(Input_MoleculeParentFolder_interval, GroupNum, MoleculeNum, threshold)
				for straightness in straightness_list:
					[FrameNum, theta_net, theta_max, test_score_sum] = straightness
					print "group", GroupNum, "molecule", MoleculeNum, "frame", FrameNum
					fileid.write(" ".join([mole[0:19],str(FrameNum), str(theta_net), str(theta_max), str(test_score_sum), "\n"]))
			else:
				print "not exist!"
		else:
			print "not exist"
	fileid.close()

	
#####################################################################################
## calling main	

if __name__ == "__main__":
	main(sys.argv[1:])











