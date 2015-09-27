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

import lib_frames_index as fi  		# model1: get access to frames in database
import lib_molecule_file_control as mfc 	# model2: operation about molecule files

###################################################################################
## FUNCTION: get the 1-pixel backbones from 5-pixel backbones for one interval

## $python from5to1.py IntervalNum

###################################################################################
## global variables

Input_ListFolder = "/exports/aspen/cwu/Project_QualityScore/Data/molecule_file"
Input_moleParentFolder = "/exports/home/cwu269/cwu/Project_QualityScore/Data/molecule_file/5pixel_region"
Output_moleParentFolder = "/exports/home/cwu269/cwu/Project_QualityScore/Data/molecule_file/1pixel_region"
  
onepixel = 1

###################################################################################
## func1: get the 1 pixel backbone molecule file from a 5 pixels backcone molecule file,
##        given groupNum, moleculeNum, Input_moleParentFolder and Output_moleParentFolder 

def npixels_file_to_1pixel_file(Input_moleParentFolder, Output_moleParentFolder, GroupNum, MoleculeNum):
	GroupName = "group1-" + str(GroupNum) + "-inca34-outputs/"
	
	MoleculeName = "molecule" + str(MoleculeNum) + ".txt"
	molecule_file_5pixel = os.path.join(Input_moleParentFolder, GroupName, MoleculeName)
	molecule_file_1pixel = os.path.join(Output_moleParentFolder, GroupName, MoleculeName)
	
	if os.path.exists(molecule_file_5pixel):
		
		fileContents_5pixel = mfc.read_molecule_files(molecule_file_5pixel)
        	numPixels = mfc.get_numPixels(fileContents_5pixel)
        	numLines = len(fileContents_5pixel)
        	SkipMolecule = mfc.checkforMinusOne(fileContents_5pixel, numPixels)
        	if SkipMolecule:		## If there is a -1 in the third column, skip that molecule
         	   exit(1)
	
       	 	FrameNumbers = mfc.get_FrameNumbers(fileContents_5pixel)
        	CellNumbers = mfc.get_CellNumbers(fileContents_5pixel)

		frames = list(set(FrameNumbers))
		grayDataList = []
		for FrameNum in frames:
			grayDataList.append(fi.get_grayData_frominf(GroupNum, FrameNum))
		onepixel = 1
	
	
		Output_GroupFolder = os.path.join(Output_moleParentFolder, GroupName)
		if(not(os.path.exists(Output_GroupFolder))):
			os.makedirs(Output_GroupFolder)
		
		fileWrite = open(molecule_file_1pixel,"w+")
        	fileWrite.seek(0)
        	fileWrite.write(fileContents_5pixel[0])
        		
        	for LineNum in range(1, numPixels+1):
      			Line_Old = fileContents_5pixel[LineNum]
        		FrameNumIndex = frames.index(FrameNumbers[LineNum-1])
        		Line_New = str(mfc.change_Intensity(Line_Old, onepixel, grayDataList[FrameNumIndex]))
        		fileWrite.write(Line_New + "\n")
	
	        for LineNum in range(numPixels+1, numLines):
	        	fileWrite.write(fileContents_5pixel[LineNum])
	         	
	        fileWrite.close()
	else:
		print  molecule_file_5pixel, "no input!"	
	
###################################################################################
## Main                                                                          ##
###################################################################################
def main(argv):

    	try:                                
     	   	opts, args = getopt.getopt(argv, "hg:d", ["help", "grammar="])
    	except getopt.GetoptError:          
		usage()                         
      	  	sys.exit(2)      
    
	IntervalNum = int(argv[0])
	IntervalName = "refFrag_" + str(IntervalNum)
#	ChrNum = int(argv[1])
	ChrNum = 1	

	ListFileName = 'Chr' + str(ChrNum) + '_interval' + str(IntervalNum) + '.txt'
	FileReadPath = Input_ListFolder + "/" + ListFileName
	fileid = open(FileReadPath, 'r')
	mole_list = fileid.read()
	mole_list = mole_list.split(' ')
	fileid.close()

	for mole in mole_list:
		GroupNum = int(mole[0:7])
		MoleculeNum = int(mole[15:19])
		
		print GroupNum, MoleculeNum
		
		Input_moleParentFolder_interval = os.path.join(Input_moleParentFolder, IntervalName)
		Output_moleParentFolder_interval = os.path.join(Output_moleParentFolder, IntervalName)

		npixels_file_to_1pixel_file(Input_moleParentFolder_interval, Output_moleParentFolder_interval, GroupNum, MoleculeNum)	

	print "end !"	



###############################################################################
	
if __name__ == "__main__":
 	main(sys.argv[1:])

