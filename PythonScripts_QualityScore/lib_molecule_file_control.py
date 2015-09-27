#! /home/nandi/anaconda/bin/python

import getopt
import gzip
import io
import MySQLdb

from numpy import * 
import os
import pickle
import shutil
import struct
import sys

from string import translate, maketrans, punctuation 
from os import path


###################################################################################
## model2: operation about molecule files

version = "1.0"

###################################################################################
## func0: get the fullpath of one group containing molecule files
def get_grouppath(GroupNum, ParentFolder):
	GroupNum = str(GroupNum)
    	GroupFolder = "group1-" + GroupNum + "-inca34-outputs"
    	GroupFolder = os.path.join(ParentFolder, GroupFolder)
    	return GroupFolder

###################################################################################
## func1_0: get fullpath&name of a molecule file, given group number and molecule number
def get_molefilepath(GroupNum, MoleNum, ParentFolder):
    	GroupFolder = "group1-" + GroupNum + "-inca34-outputs"
   	GroupFolder = os.path.join(ParentFolder, GroupFolder)
    	MoleFilename = "molecule" + MoleNum + ".txt"
    	MoleFilePath = os.path.join(GroupFolder, MoleFilename)
    	return MoleFilePath

###################################################################################
## func1_1: get all the molecule filenames in one group as a list, given the group number and input folder
def get_molecule_filenames(GroupFolder):
    	molecule_filenames = [f for f in sorted(os.listdir(GroupFolder)) if 'molecule' in f]
    	return molecule_filenames


###################################################################################
## func1_2: get fullpath&name of a molecule file, given a group folder and molecule number
def get_molecule_filepath(GroupFolder, MoleNum = 0, molecule_filename = ""):
	MoleNum = str(MoleNum)
	molecule_filenames = get_molecule_filenames(GroupFolder)
	if(molecule_filename == ""):
		molecule_filename = "molecule" + MoleNum + ".txt"
	molecule_filepath = os.path.join(GroupFolder, molecule_filename)
	for f in molecule_filenames:
		if(f==molecule_filename):
			return molecule_filepath
	print "something wrong! maybe the file doesn't exist...."
	return molecule_filepath


###################################################################################
## func1_3: read the content in a molecule file, inputting fullpath&name of the file, and readlines as mode

def read_molecule_files(filepathRead):
    if os.path.exists(filepathRead) and os.path.isfile(filepathRead):
        fileRead = open(filepathRead,"r")
       	fileContents = fileRead.readlines()
        fileRead.close()
        return fileContents 


###################################################################################
## summary1: read-in molecule file from name, return fileContents
def read_molecule_files_fromname(ParentFolder, GroupNum, molecule_filename):
	GroupFolder = get_grouppath(GroupNum, ParentFolder)
	filepathRead = get_molecule_filepath(GroupFolder, 0, molecule_filename)
	fileContents = read_molecule_files(filepathRead)
	return fileContents


###################################################################################
## summary2: read-in molecule file from information, return fileContents

def read_molecule_files_frominf(ParentFolder, GroupNum, MoleculeNum):
	GroupFolder = get_grouppath(GroupNum, ParentFolder)
	filepathRead = get_molecule_filepath(GroupFolder, MoleculeNum)
	fileContents = read_molecule_files(filepathRead)
	return fileContents



###################################################################################
## func2_0: Check if there is a -1 in the third column of the molecule.txt files in Prabu's folder.

def checkforMinusOne(fileContents, numPixels):
    for LineNum in range(0, numPixels):
	IntensityValue = int(fileContents[LineNum].split(" ")[2])
        if IntensityValue < 0:
            SkipMolecule = 1
            break
        else:
            SkipMolecule = 0
    return SkipMolecule 

###################################################################################
## func2_1: get the number of pixels in this molecule file, inputting content list of the file

def get_numPixels(fileContents):
    Line1 = fileContents[0].split(" ")
    numPixels = int(Line1[0])
    return numPixels 

###################################################################################
## func2_2: get the frame number of certain line in a molecule file, inputting this line as str

def get_FrameNum(Line):
    Line2 = Line.split(" ")
    FrameNum = int(Line2[0])
    return FrameNum 

###################################################################################
##func2_3:get all frame numbers for a molecule file as a numeric list

def get_FrameNumbers(fileContents):
    numPixels = get_numPixels(fileContents)
    FrameNumList = []
    for LineNum in range(1, numPixels+1):
	Line2 = fileContents[LineNum].split(" ")
        FrameNumList.append(int(Line2[0]))

    return FrameNumList

###################################################################################
## func2_4: get the cell number of certain line in a molecule file, inputting this line as str
def get_CellNum(Line):
    Line2 = Line.split(" ")
    CellNum = int(Line2[1])
    return CellNum 


###################################################################################
##func2_5:get all cell numbers for a molecule file as a numeric list

def get_CellNumbers(fileContents):

    	CellNumList = []
    	Line0 = fileContents[0].split(" ")
    	numpixels = int(Line0[0])

    	for l in range(1, numpixels + 1):
    		Line2 = fileContents[l].split(" ")
    		CellNumList.append(int(Line2[1]))

    	return CellNumList 


####################################################################################
## func2_6: get the intensity value of each line in a molecule file, inputting this line as str

def get_Inten(Line):
	Line2 = Line.split(" ")
	Intensity = int(Line2[2])
	return Intensity

####################################################################################
## func2_7: get the intensity value of all lines in a molecule file as a numeric list

def get_IntenAll(fileContents):
	numPixels = get_numPixels(fileContents)
	IntensityList = []
	for LineNum in range(1, numPixels+1):
		IntensityList.append(get_Inten(fileContents[LineNum]))

	return IntensityList

#####################################################################################
## func2_8: determine whether this molecule steps across frame, and return the number of frame(s and positions of frame connection)

def get_frame_count_mole(fileContents):
    connect_count = 0
    cellnumlist = get_CellNumbers(fileContents)
    framenumlist = get_FrameNumbers(fileContents)
    frames = [framenumlist[0]]
    cells = []

    for ind in range(0, (len(framenumlist)-2)):
	if(framenumlist[ind] != framenumlist[ind+1]):
		connect_count += 1 
		frames.append(framenumlist[ind+1])
		cells.append([cellnumlist[ind],cellnumlist[ind+1]])

    return connect_count, frames, cells

#####################################################################################
## func2_9_1: determine whether a molecule is in the frame or not

def find_mole_setframe(fileContents, frame):
	numPixels = get_numPixels(fileContents)
	inframe = 0
	for i in range(1, numPixels + 1):
		Line_franmenum = get_FrameNum(fileContents[i])
		if(Line_franmenum == frame):
			inframe = 1
	return inframe

#####################################################################################
## func2_9_2: get the content of a molecule in one frame

def get_content_setframe(fileContents, frame):
	numPixels = get_numPixels(fileContents)
	Contents_setframe = [fileContents[0]]
	newPixels = 0
	for i in range(1, numPixels + 1):
		Line_franmenum = get_FrameNum(fileContents[i])
		if(Line_franmenum == frame):
			newPixels += 1
			Contents_setframe.append(fileContents[i])
	Line0 = fileContents[0].split(" ")
	Line0[0] = str(newPixels)
	Contents_setframe[0] = " ".join(Line0)
	
	return Contents_setframe
		

#####################################################################################
## func2_10: get molecule names in one frame

def get_names_setframe(GroupFolder, frame):
	molecule_filenames = get_molecule_filenames(GroupFolder)
	molecule_filenames_setframe = []
	for filename in molecule_filenames:
		filepathRead = get_molecule_filepath(GroupFolder, molecule_filename = filename)
		fileContents = read_molecule_files(filepathRead)
		inframe = find_mole_setframe(fileContents, frame)
		if(inframe):
			molecule_filenames_setframe.append(filename)
	if(molecule_filenames_setframe == []):
		print "no molecule in this group involving this frame!"
	return molecule_filenames_setframe
			

#####################################################################################
## func3_1: get the coordinate of one cell

def get_coordinate(cellnum,ncol):
	coordinate = zeros(2,dtype=int)
	coordinate[0] = cellnum /ncol 	###    # of row 
	coordinate[1] = cellnum % ncol	###    # of col
	coordinate = list(coordinate)
	return coordinate

#####################################################################################
## func3_2: get the coordinate of all points in a molecule

def get_map(ncol,CellNumList):
	mole_map = [] 
	for cell in CellNumList:
		mole_map.append(get_coordinate(ncol, cell))
	return mole_map


#####################################################################################
## func3_3: get the cellnum of one coordinate

def get_index(row, column, ncol):
	cellnum = row * ncol + column
	return cellnum

#####################################################################################
## func4_1: get the top and bottom position of a molecule in certain frame

def boundary_y(ncol, fileContents, Frame):
	# get molecule information ready
	CellNums = array(get_CellNumbers(fileContents))
	FrameNums = array(get_FrameNumbers(fileContents))
	
	try:
		# select points in this frame
		indices = where(FrameNums == Frame)
		CellNums_inframe = CellNums[indices]
	
		# get the position of these points
		Mole_Map = array(get_map(ncol, CellNums_inframe))

		# calculate the boundary
		top_y = Mole_Map.min(axis = 0)[1]
		bottom_y = Mole_Map.max(axis = 0)[1]	
  
	#####  additional function 		
	#	Mole_Map_y = Mole_Map[:,1]
	#
	#	top_indices = where(Mole_Map_y == top_y)
	#	bottom_indices = where(Mole_Map_y == bottom_y)
	#	
	#	top_position = Mole_Map[top_indices]
	#	bottom_position = Mole_Map[bottom_indices]

		return top_y, bottom_y  #, top_position, bottom_position 
	except:
		print "not in this frame!"


#####################################################################################
## func4_2: get the head and tail position of a molecule in certain frame

def boundary_x(ncol, CellNumList):

	head_cells = []
	tail_cells = []
	CellPosition = []

	for cellnum in CellNumList:
		CellPosition.append(get_coordinate(cellnum, ncol))
	Cell_x = array(CellPosition)[:,1]

	cell_x_min = Cell_x.min()
	cell_x_max = Cell_x.max()

	indices_min = where(Cell_x == cell_x_min)
	indices_max = where(Cell_x == cell_x_max)

	for ind in indices_min:
		head_cells.append(CellNumList[ind])
	for ind in indices_max:
		tail_cells.append(CellNumList[ind])

	return head_cells, tail_cells

#####################################################################################
## func4_3: determine a cell is in a region or not

def find_cell_inregion_x(cellnum, select_region, ncol):
	cell_position = get_coordinate(cellnum, ncol)
	inregion = 1
	if(cell_position[1] < select_region[0] or cell_position[1] > select_region[1]):
		inregion = 0
	return inregion


###################################################################################
## func5_1: calculate the sum intensity of ventral and dorsal pixels of original pixel 

def getNewIntensity(onepixel, CellNum, grayData):
   
    nRow = grayData.shape[0]
    nCol = grayData.shape[1]
    
    CellNumber = CellNum
    CellRow = (CellNumber/nCol)
    CellCol = CellNumber % nCol
    Intensity5 = grayData[CellRow,CellCol] + grayData[CellRow-1,CellCol] + grayData[CellRow-2,CellCol] + grayData[CellRow+1,CellCol] + grayData[CellRow+2,CellCol] 
    Intensity3 = grayData[CellRow,CellCol] + grayData[CellRow-1,CellCol] + grayData[CellRow+1,CellCol]
    Intensity1 = grayData[CellRow,CellCol]
    if onepixel == 1:
        NewIntensity = Intensity1
    elif onepixel == 3:
        NewIntensity = Intensity3
    elif onepixel == 5:
        NewIntensity = Intensity5
    else:
        NewIntensity = Intensity5
        
    return NewIntensity

###################################################################################
## func5_2: change the intensity information in one line of molecule file

def change_Intensity(Line, onepixel, grayData):
    LineStr = Line.split(" ")
    LineInt = map(int, LineStr)
    FrameNumber = LineInt[0]
    CellNumber = LineInt[1]
    LineInt[2] = getNewIntensity(onepixel, CellNumber, grayData)
    LineStr_New = map(str, LineInt)
    Line_New = " ".join(LineStr_New)
    return Line_New

######################################################################################
## func_end: return name of this model

def modelName():
	return __name__






















