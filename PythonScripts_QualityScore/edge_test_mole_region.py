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
import math

from string import translate, maketrans, punctuation 
from os import path

import lib_frames_index as fi  		# model1: get access to frames in database
import lib_molecule_file_control as mfc 	# model2: operation about molecule files
import lib_get_mole_surrounding as gms


###################################################################################################
## global variables

Input_MoleParentFolder = "/exports/aspen/cwu/Project_QualityScore/Data/molecule_file/5pixel_region"
#Input_ListFolder = "/exports/aspen/cwu/Project_QualityScore/Data/molecule_file"

ncol = 1392
nrow = 1024

edge_range = 50

###################################################################################################
## func1
def select_outframe(GroupNum, FrameNum):
	grayData_array = fi.get_grayData_frominf(GroupNum, FrameNum)
	
	Zeros_Position = np.where(grayData_array == 0)	## [0]: row
							## [1]: column
	
	Zeros_row = list(Zeros_Position[0])
	Zeros_column = list(Zeros_Position[1])
	
	nCol = list(np.repeat(ncol, len(Zeros_column)))

	Zeros_cellnum = map(mfc.get_index, Zeros_row, Zeros_column, nCol)
	
	return Zeros_cellnum

###################################################################################################
## func2
def eage_test_OnePixel(OutEdge_CellList, origin_cellnum):

	upper_bound, upper_condition = gms.target_pixel(origin_cellnum, edge_range, 0)
	lower_bound, lower_condition = gms.target_pixel(origin_cellnum, -edge_range, 0)

	if(upper_condition and lower_condition):
		return max(OutEdge_CellList.count(upper_bound), OutEdge_CellList.count(lower_bound))
	else:
		return 1

###################################################################################################
## func3

def edge_test_Oneframe(fileContents, GroupNum, FrameNum, MoleculeNum):

	OutEdge_CellList = select_outframe(GroupNum, FrameNum)
	
	Contents_setframe = mfc.get_content_setframe(fileContents, FrameNum)
	CellNumList = mfc.get_CellNumbers(Contents_setframe)
	
	OutOfEdge_Positions = []

	for cellnum in CellNumList:
		if(eage_test_OnePixel(OutEdge_CellList, cellnum)):
			OutOfEdge_Positions.append(cellnum)
		
	return OutOfEdge_Positions

###################################################################################################
## func4

def edge_test(Input_MoleParentFolder_interval, GroupNum, MoleculeNum):
	
	fileContents = mfc.read_molecule_files_frominf(Input_MoleParentFolder_interval, GroupNum, MoleculeNum)
	connect_count, frames, cells = mfc.get_frame_count_mole(fileContents)

	edge = []
	for FrameNum in frames:
		print "group", GroupNum, "molecule", MoleculeNum, "frame", FrameNum
		OutOfEdge_Positions = edge_test_Oneframe(fileContents, GroupNum, FrameNum, MoleculeNum)
		if len(OutOfEdge_Positions):
			print "at the edge!  ", len(OutOfEdge_Positions), "pixels within ", str(edge_range)+"-pixel area of the edge"
			edge.append(1)
		else:
			print "not at the edge..."
			edge.append(0)
		
	return edge

###################################################################################################
## calling main	

def main(argv):
	
	try:                                
        	opts, args = getopt.getopt(argv, "hs:e:", ["help", "grammar="])
    	except getopt.GetoptError:          
        	usage()                         
        	sys.exit(2)   

	IntervalNum = int(args[0])
	GroupNum = int(args[1])
	MoleculeNum = int(args[2])

	IntervalName = "refFrag_" + str(IntervalNum)
	Input_MoleParentFolder_interval = os.path.join(Input_MoleParentFolder, IntervalName)

	edge = edge_test(Input_MoleParentFolder_interval, GroupNum, MoleculeNum)

	print edge ## list for edge test, 1: within edge area; 0: out of edge area
	
###################################################################################################
## calling main	

if __name__ == "__main__":
	main(sys.argv[1:])

	


