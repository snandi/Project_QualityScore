#! /home/nandi/anaconda/bin/python

import getopt
import io
import re
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

import lib_frames_index as fi  			# model1: get access to frames in database
import lib_molecule_file_control as mfc 	# model2: operation about molecule files



###################################################################################
## FUNCTION: get the molecule files containing the region aligned to one interval

## $python get_molefiles_region.py IntervalNum


################################################################################
## global variables

convert = 209.0

Input_moleParentFolder = "/aspen/nandi/MF_cap348/maps_inca34"
Output_moleParentFolder = "/exports/home/cwu269/cwu/Project_QualityScore/Data/molecule_file/5pixel_region"
Input_inforParentFolder='/exports/aspen/cwu/Project_QualityScore/Data/molecule_file/'
InforFileName = 'MF_cap348_inca34_cf209_minSize50_minFrag5_Aligned.alignmentChunks'


#################################################################################
## convert base pair coordinates to pixel indices
def convert_bp2pixel(bp_coordinate, LengthRatio):
	bp_per_pixel = convert / LengthRatio
	if(bp_coordinate == 0):
		pixel_coordinate = 0
	else:
		pixel_coordinate = int(math.floor((bp_coordinate-1) / bp_per_pixel))
	
	return pixel_coordinate
	
#################################################################################
## generate the molecule file containing certain region(aligned)
def get_region_file(Input_moleParentFolder, Output_moleParentFolder, GroupNum, MoleculeNum,start_pixelind, end_pixelind, orientation):
	GroupName = "group1-" + str(GroupNum) + "-inca34-outputs/"
	
	MoleculeName = "molecule" + str(MoleculeNum) + ".txt"

	molecule_file_whole = os.path.join(Input_moleParentFolder, GroupName, MoleculeName)
	molecule_file_region = os.path.join(Output_moleParentFolder, GroupName, MoleculeName)
	
	fileContents_whole = mfc.read_molecule_files(molecule_file_whole)
	
        numPixels_whole = mfc.get_numPixels(fileContents_whole)
        numLines_whole = len(fileContents_whole)

        SkipMolecule = mfc.checkforMinusOne(fileContents_whole, numPixels_whole)
        if SkipMolecule:		## If there is a -1 in the third column, skip that molecule
            exit(1)

	if(start_pixelind+1>numPixels_whole or end_pixelind+1>numPixels_whole):
#		print "out of range!"
#		print  start_pixelind, end_pixelind, numPixels_whole
		return numPixels_whole - max(start_pixelind+1,end_pixelind+1)

	Pixels_whole = fileContents_whole[1:(numPixels_whole + 1)]
	


	Output_GroupFolder = os.path.join(Output_moleParentFolder, GroupName)
	if(not(os.path.exists(Output_GroupFolder))):
		os.makedirs(Output_GroupFolder)
	if(not(os.path.exists(Output_GroupFolder))):
		print "create problem"
	
	os.chdir(Output_GroupFolder)
#	print os.getcwd()

	fileWrite = open(molecule_file_region,"w+")
        fileWrite.seek(0)
	Line0 = fileContents_whole[0].split(" ")
	Line0[0] = str(abs(start_pixelind - end_pixelind) + 1)
        fileWrite.write(" ".join(Line0))
        
	if orientation>0 :
		step = 1 
		for LineNum in range(start_pixelind+1, end_pixelind+2, step):
      			Line = fileContents_whole[LineNum]
        		fileWrite.write(Line)
	else:
		step = -1
		for LineNum in range(start_pixelind+1, end_pixelind, step):
      			Line = fileContents_whole[LineNum]
        		fileWrite.write(Line)

        for LineNum in range(numPixels_whole+1, numLines_whole):
        	fileWrite.write(fileContents_whole[LineNum])
         	
        fileWrite.close()
	
	return 1


#################################################################################
## main: get the molecule files containing the region aligned to one interval,
## 	 given interval number(and Chromosome number)	 
def main(argv):
	
	try:                                
     	   	opts, args = getopt.getopt(argv, "hg:d", ["help", "grammar="])
    	except getopt.GetoptError:          
		usage()                         
      	  	sys.exit(2)      
	
	IntervalNum = int(argv[0])
#	ChrNum = int(argv[1])
	ChrNum = 1

	## get the molecule list
	Molecule_ListFileName = 'Chr' + str(ChrNum) + '_interval' + str(IntervalNum) + '.txt'
	Molecule_ListFile = Input_inforParentFolder + Molecule_ListFileName
	mp = open(Molecule_ListFile,"r")
	mole_list = mp.read()
	mole_list = mole_list.split(' ')
	mp.close()
	
	print len(mole_list)

	## get information about alignments	
	InforReadPath = Input_inforParentFolder + InforFileName
	fp = open(InforReadPath, "r")
	content = fp.read()
	fp.close()

	## processing
	k = 0
	for mole in mole_list:

		moleID = str(re.findall(r"(.*)", mole)[0])
		GroupNum = int(moleID[0:7])
		MoleculeNum = int(moleID[16:19])
		
		regex = 'chr' + str(ChrNum) + '\s' + str(IntervalNum) + '\s' + str(IntervalNum+1) + '\s' + moleID + '\s[0-9]+\s[0-9]+\s[0-9]+\s[0-9]+\s([0-9]*)\s([0-9]*)\s(.+?)\s(.*)'
		s = re.findall(regex, content)

		for region in s:
			region = list(region)
			start_bp = int(region[0])
			end_bp = int(region[1])
			orientation = int(region[2])
			LengthRatio = float(region[3])

			start_pixelind = convert_bp2pixel(start_bp, LengthRatio)
			end_pixelind = convert_bp2pixel(end_bp, LengthRatio)
			a = get_region_file(Input_moleParentFolder, Output_moleParentFolder, GroupNum, MoleculeNum,start_pixelind, end_pixelind, orientation)
			if(a < 0): 
				print moleID, "   discrepancy:  " ,-a
				k += 1


	print k, "out of range"
	

#################################################################################
## calling main function

if __name__=="__main__":
	main(sys.argv[1:])




