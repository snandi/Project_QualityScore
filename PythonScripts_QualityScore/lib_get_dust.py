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
import matplotlib.pyplot as plt

from string import translate, maketrans, punctuation 
from os import path
from matplotlib.legend_handler import HandlerLine2D

import lib_frames_index as fi  			# model1: get access to frames in database
import lib_molecule_file_control as mfc 	# model2: operation about molecule files
import lib_get_mole_surrounding as gms 		# model3: get the surrounding pixels(donut or sandwich) around molecule backbone

###################################################################################
## model4: detect the nearby noisy pixels and the backbone pixels effected by garbage

################################################################################
## global variables

version = "1.0"

ncol = 1392
nrow = 1024

layerMAX = 3
effect_range = 5


################################################################################
## func0: multiply for list
def multiply(x,y):
	return x * y


################################################################################
## func1_1: get cellnumber of surrounding pixels for certain molecule/fragment
def get_donut_coordinates(SurroundingfileContents, layer):
	CellNumList = []
	for p in range(1, len(SurroundingfileContents)):
		line = SurroundingfileContents[p].split(" ")
		cellnum = int(line[layerMAX + layer - 1])
		CellNumList.append(cellnum)

	if CellNumList.count(0) != 0:
		CellNumList.remove(0)

	return CellNumList


################################################################################
## func1_2: get intensity of surrounding pixels for certain molecule/fragment
def get_donut_intensity(SurroundingfileContents, layer):
	IntensityList = []
	for p in range(1, len(SurroundingfileContents)):
		line = SurroundingfileContents[p].split(" ")
		intensity = int(line[layer - 1])
		IntensityList.append(intensity)

	if IntensityList.count(0) != 0:	
		IntensityList.remove(0)

	return IntensityList

################################################################################
## func2_1: get cellnumber of garbage pixels for certain molecule/fragment

def get_garbage_cooridnate(SurroundingfileContents, layer, threshold):
	
	GarbageList_coordinate = []
	CellNumList = get_donut_coordinates(SurroundingfileContents, layer)
	IntensityList = get_donut_intensity(SurroundingfileContents, layer)

	for i in range(0, len(CellNumList)):
		if IntensityList[i] > threshold:
			GarbageList_coordinate.append(CellNumList[i])

	return GarbageList_coordinate

################################################################################
## func2_2: get intensity of garbage pixels for certain molecule/fragment
def get_garbage_intensity(SurroundingfileContents, layer, threshold):
	
	IntensityList = get_donut_intensity(SurroundingfileContents, layer)

	GarbageList_intensity = filter(lambda x:x > threshold, IntensityList)
	
	return GarbageList_intensity

################################################################################
## func3_1: get cellnumber of pixels within the effect range of certain pixel
def effected_pixels_foronepoint(origin_cellnum):
	EffectedList = []
	for hori_dis in range(-effect_range, effect_range + 1):
		[arm_pixel, arm_condition] = gms.target_pixel(origin_cellnum, 0, hori_dis)
		for vent_dis in range(-effect_range, effect_range + 1):
			[effected_pixel, condition] = gms.target_pixel(arm_pixel, vent_dis, 0)
			if condition:
				EffectedList.append(effected_pixel)
	
	return EffectedList
	

################################################################################
## func3_2: get cellnumber of pixels within the effect range of all garbage pixels
def effected_pixels(SurroundingfileContents, layer, threshold):
	
	GarbageList_coordinate = get_garbage_cooridnate(SurroundingfileContents, layer, threshold)

	EffectedListAll = []
	for dust_pixel in GarbageList_coordinate:
		effected_list_local = effected_pixels_foronepoint(dust_pixel)
		EffectedListAll.extend(effected_list_local)
	
	EffectedListAll = list(set(EffectedListAll))

	return EffectedListAll
		
################################################################################
## func4_1: determine one pixel is within certain region or not
def determine_within_foronepoint(test_cellnum, EffectedList):
	for effected_pixel in EffectedList:
		if test_cellnum == effected_pixel:
			return 1
	
	return 0


################################################################################
## func4_2: determine which pixels in backbone are effected by surrounding garbage
def determine_within_forbackbone(BackbonefileContents, SurroundingfileContents, layer, threshold):
	CellNumList_backbone = mfc.get_CellNumbers(BackbonefileContents)
	EffectedList = effected_pixels(SurroundingfileContents, layer, threshold)
	
	ConditionList = []
	for test_cellnum in CellNumList_backbone:
		condition = determine_within_foronepoint(test_cellnum, EffectedList)
		ConditionList.append(condition)
	
	return ConditionList

################################################################################
## func_application1: get the effected condition of one molecule/fragment

def Effected_OneMolecule(Input_MoleParentFolder_interval, Input_SurroundingParentFolder_interval, GroupNum, MoleculeNum, layer, threshold):

#def main(argv):

#	try:                                
#        	opts, args = getopt.getopt(argv, "hs:e:", ["help", "grammar="])
#    	except getopt.GetoptError:          
#        	usage()                         
#        	sys.exit(2)   

#	GroupNum = int(args[0])	
#	MoleculeNum = int(args[1])
#	layer = int(args[2])
	
	GroupNum = int(GroupNum)	
	MoleculeNum = int(MoleculeNum)
	layer = int(layer)

	print GroupNum, MoleculeNum

	GroupName = "group1-" + str(GroupNum) + "-inca34-outputs/"

	
	
	GroupFolder = mfc.get_grouppath(GroupNum, Input_MoleParentFolder_interval)
	if os.path.exists(GroupFolder):
		molecule_fileContents = mfc.read_molecule_files_frominf(Input_MoleParentFolder_interval, GroupNum, MoleculeNum)
		connect_count, frames, cells = mfc.get_frame_count_mole(molecule_fileContents)
		for frame in frames:
			BackbonefileContents = mfc.get_content_setframe(molecule_fileContents, frame)
			
			FrameNum = frame
			surrounding_fileName = "molecule" + str(MoleculeNum) + "_frame" + str(FrameNum) + ".txt"
			SurroundingReadPath = os.path.join(Input_SurroundingParentFolder_interval, GroupName, surrounding_fileName)
			SurrID = open(SurroundingReadPath, 'r')
			SurroundingfileContents = SurrID.readlines()

			ConditionList = determine_within_forbackbone(BackbonefileContents, SurroundingfileContents, layer, threshold)
			

			## plotting
			fig = plt.figure(frame)
			BackboneIntensity = mfc.get_IntenAll(BackbonefileContents)
		
			attention_x = list(map(multiply, range(0, len(BackboneIntensity)), ConditionList))
			attention_x = filter(lambda x:x is not 0, attention_x)
			if ConditionList[0] != 0:
				attention_x.insert(0, 0)

			attention_y = map(multiply, BackboneIntensity, ConditionList)
			if attention_y.count(0) != 0:
				attention_y = filter(lambda y:y is not 0, attention_y)
		
			line, = plt.plot(BackboneIntensity, label = "intensity porfile")
			dot, = plt.plot(attention_x, attention_y, "ro", color = "red")

			Title = "Group" + str(GroupNum) + "_Molecule" + str(MoleculeNum) + " Backbone1 intensity" + "(default threshold" + str(threshold) + ")"
			plt.title(Title)			

			plt.xlabel("position(orientation depend on alignment)/pixels")
			plt.ylabel("intensity")

#			plt.legend([line], "intensity porfile", loc = 4)
			plt.legend([line, dot],  ["intensity porfile", "effected pixels(from layer"+str(layer) + ")"])
#			print ConditionList

		plt.show()
	else:
		print "group not exist!"




################################################################################	
#if __name__ == "__main__":
#	main(sys.argv[1:])


