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
## FUNCTION: get the straightness measure of a molecule fragment with simplify algorithm


##############################################
## global variables
Input_MoleParentFolder = "/exports/aspen/cwu/Project_QualityScore/Data/molecule_file/5pixel_region"

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
	GroupNum = int(argv[1])
	MoleculeNum = int(argv[2])

	IntervalName = "refFrag_" + str(IntervalNum)
	Input_MoleParentFolder_interval = os.path.join(Input_MoleParentFolder, IntervalName)

	straightness_list = sl.simplify_lines_AllFrame(Input_MoleParentFolder_interval, GroupNum, MoleculeNum, threshold)
		
	for straightness in straightness_list:
		[FrameNum, theta_net, theta_max, test_score_sum] = straightness

		print "group", GroupNum, "molecule", MoleculeNum, "frame", FrameNum
		print "net of angles:", theta_net
		print "max of angles:", theta_max
		print "score based on net & max:", test_score_sum

#####################################################################################
## calling main	

if __name__ == "__main__":
	main(sys.argv[1:])











