#! /home/nandi/anaconda/bin/python

import getopt
import gzip
import io
import re
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

###################################################################################
## model3: functions for getting the surrounding pixels(donut or sandwich) around molecule backbone
## NOTE!! problem: because most molecule files are truncation,
##		   surrounding pixels might overlap with nearby molecule backbone...

version = "1.0"

## constants: shape of frame

ncol = 1392
nrow = 1024

colormap = 'gist_ncar'
image_max = 30000
mask = 24000


###############################################################################
##func1: calculate the cell number of target surrounding pixel
def target_pixel(origin_cellnum, vertical_dis = 0, horizontal_dis = 0):
	origin_position = mfc.get_coordinate(origin_cellnum, ncol)
	inframe = 1
	target_position = [origin_position[0] + vertical_dis, origin_position[1] + horizontal_dis]
	target_cellnum = target_position[0]*ncol + target_position[1]
	
	if (target_position[0] < 0) or (target_position[0] >= nrow) or (target_position[1] < 0) or (target_position[1] >= ncol):
		inframe = 0
	return [target_cellnum, inframe]

################################################################################
## func2: get the intensity of target pixel from frame data

def get_pixel_intensity(target_cellnum, grayData_array):
	target_position = list(mfc.get_coordinate(target_cellnum, ncol))
	target_intensity = grayData_array[target_position[0], target_position[1]]
	return target_intensity

################################################################################
## func3: get head or tail surroundings
def head_tail_box(origin_cellnum, head = True, tail = False, vertical_range = 3, horizontal_range = 3):
	target_cells = []
	if(head):
		box_central = target_pixel(origin_cellnum, 0, -horizontal_range)
		target_cells.append(box_central)
		for i in range(1, vertical_range + 1):
			target_cell = target_pixel(origin_cellnum, i, -horizontal_range)
#			target_cell.append(origin_cellnum)
			target_cells.append(target_cell)

			target_cell = target_pixel(origin_cellnum, -i, -horizontal_range)
#			target_cell.append(origin_cellnum)
			target_cells.append(target_cell)

		for j in range(1, horizontal_range + 1):
			target_cell = target_pixel(origin_cellnum, vertical_range, -j)
#			target_cell.append(origin_cellnum)
			target_cells.append(target_cell)

			target_cell = target_pixel(origin_cellnum, -vertical_range, -j)
#			target_cell.append(origin_cellnum)
			target_cells.append(target_cell)

	if(tail):
		box_central = target_pixel(origin_cellnum, 0, horizontal_range)
		target_cells.append(box_central)
		for i in range(1, vertical_range + 1):
			target_cell = target_pixel(origin_cellnum, i, horizontal_range)
#			target_cell.append(origin_cellnum)
			target_cells.append(target_cell)

			target_cell = target_pixel(origin_cellnum, -i, horizontal_range)
#			target_cell.append(origin_cellnum)
			target_cells.append(target_cell)

		for j in range(1, horizontal_range + 1):

			target_cell = target_pixel(origin_cellnum, vertical_range, j)
#			target_cell.append(origin_cellnum)
			target_cells.append(target_cell)

			target_cell = target_pixel(origin_cellnum, -vertical_range, j)
#			target_cell.append(origin_cellnum)
			target_cells.append(target_cell)

	return target_cells
	

#################################################################################
## func4: get ventral and dorsal surrundings
def ventral_dorsal_string(CellNumList, vertical_range = 3):
	target_cells = []
	for origin_cellnum in CellNumList:
		target_cell = target_pixel(origin_cellnum, vertical_range, 0)
#		target_cell.append(origin_cellnum)
		target_cells.append(target_cell)

		target_cell = target_pixel(origin_cellnum, -vertical_range, 0)
#		target_cell.append(origin_cellnum)
		target_cells.append(target_cell)
	return target_cells

#################################################################################
## func5: get all the target surrounding pixels of molecule that are in this frame 
def surrounding_circle(CellNumList,backbone_range = 2, surround_range = 1, close = False):

	vertical_range = backbone_range + surround_range
	horizontal_range = vertical_range
	
	surrounding_condition = ventral_dorsal_string(CellNumList, vertical_range)

	if(close):
		head_cells, tail_cells = mfc.boundary_x(ncol, CellNumList)
		for head in head_cells:
			head_box = head_tail_box(head, True, False, vertical_range, horizontal_range)
			for head_pixel in head_box:
				surrounding_condition.append(head_pixel)
		for tail in tail_cells:
			tail_box = head_tail_box(tail, False, True, vertical_range, horizontal_range)
			for tail_pixel in tail_box:
				surrounding_condition.append(tail_pixel)

	surrounding_condition = set(map(tuple, surrounding_condition))
	surrounding_condition = list(map(list, surrounding_condition))
	
	surrounding_condition = np.array(surrounding_condition)
	
	surrounding_inframe = surrounding_condition[:,1]
	surrounding_circle = surrounding_condition.compress(surrounding_inframe, axis = 0)[:,0]

	### origin_backbone = surrounding_condition.compress(surrounding_inframe, axis = 0)[:,[0,2]]
	### ready to output the origin position(cellnums)?

	return list(surrounding_circle)  #list(surrounding_origin)

####################################################################################
## func6: get the position of surroundings

def surrounding_band(CellNumList, backbone_range = 2, surround_range = 1, close = False):
	surrounding_band = []
	for i in range(1, surround_range + 1):
		surrounding_circles = surrounding_circle(CellNumList, backbone_range, i, close)
		for sc in surrounding_cirles:
			surrounding_band.append(sc)
	
	surrounding_band = list(set(surrounding_band))

	### surrounding_origin = set(map(tuple,surrounding_band))
	### surrounding_origin = np.array(list(map(list, surrounding_band)))
	#### surrounding_band = surrounding_origin[:,0]

	######################### ready to store surrounding_origin(cellnumbs)?

	return surrounding_band

####################################################################################
## func7: get the intensity of surroundings

def surrounding_intensity_circle(CellNumList, grayData_array, grayData_surround, backbone_range = 2, surround_range = 1, close = False):
	
	surrounding_intensity_circle = []
	surrounding_circles = surrounding_circle(CellNumList, backbone_range, surround_range, close)
	
	for pixel in surrounding_circles:
		target_position = list(mfc.get_coordinate(pixel, ncol))
		surrounding_intensity_circle.append(get_pixel_intensity(pixel, grayData_array))
		grayData_surround[target_position[0], target_position[1]] = mask
		
	return surrounding_intensity_circle #, surrounding_origin

#####################################################################################
## func8_1: input GroupNum, FrameNum and MoleculeNum, 
##        get the surrounding intensity(data file and image) of one molecule in one frame
##        close -> True: donut surrounding/ False: sandwich surrounding

def surrounding_OneMolecule(Input_MoleParentFolder, Output_moleParentFolder, Output_plotParentFolder, GroupNum, FrameNum, MoleculeNum, select_region = [-1,-1], backbone_range = 2, close = False):

	GroupNum = int(GroupNum)
	FrameNum = int(FrameNum) 
	MoleculeNum = int(MoleculeNum)
	
	### data input
	grayData_array = fi.get_grayData_frominf(GroupNum, FrameNum)
	display_array1 = copy.deepcopy(grayData_array)
	display_array2 = copy.deepcopy(grayData_array)
	display_array3 = copy.deepcopy(grayData_array)

	fileContents = mfc.read_molecule_files_frominf(Input_MoleParentFolder, GroupNum, MoleculeNum)
	inframe = mfc.find_mole_setframe(fileContents, FrameNum)
	if(not(inframe)):
		print "the molecule is not in this frame! "
		exit(1)
	
	molecule_contents_setframe = mfc.get_content_setframe(fileContents, FrameNum)

	CellNumList_setframe = mfc.get_CellNumbers(molecule_contents_setframe)
	
	heads, tails = mfc.boundary_x(ncol, CellNumList_setframe)
	x_min = mfc.get_coordinate(heads[0],ncol)[1]
	x_max = mfc.get_coordinate(tails[0],ncol)[1]

	if(select_region[0] == -1): select_region[0] = x_min
	if(select_region[1] == -1): select_region[1] = x_max 
	

	CellNumList_setframe_region = []
	for cellnum in CellNumList_setframe:
		inregion = mfc.find_cell_inregion_x(cellnum, select_region, ncol)
		if(inregion):
			CellNumList_setframe_region.append(cellnum)
	
	if(len(CellNumList_setframe) != 0):
		### surrounding intensity acquisition
		surrounding_intensity_circle_1 = surrounding_intensity_circle(CellNumList_setframe_region, grayData_array, display_array1, 2, 1, close)
		surrounding_intensity_circle_2 = surrounding_intensity_circle(CellNumList_setframe_region, grayData_array, display_array2, 2, 2, close)
		surrounding_intensity_circle_3 = surrounding_intensity_circle(CellNumList_setframe_region, grayData_array, display_array3, 2, 3, close)
		
		surrounding_coordinate_1 = surrounding_circle(CellNumList_setframe_region, 2, 1, close)
		surrounding_coordinate_2 = surrounding_circle(CellNumList_setframe_region, 2, 2, close)
		surrounding_coordinate_3 = surrounding_circle(CellNumList_setframe_region, 2, 3, close)

		length = len(surrounding_intensity_circle_3)
		surrounding = np.zeros((length, 6), dtype = int)
		surrounding[0:(len(surrounding_intensity_circle_1)),0] = surrounding_intensity_circle_1
		surrounding[0:(len(surrounding_intensity_circle_2)),1] = surrounding_intensity_circle_2
		surrounding[0:(len(surrounding_intensity_circle_3)),2] = surrounding_intensity_circle_3
		surrounding[0:(len(surrounding_coordinate_1)),3] = surrounding_coordinate_1
		surrounding[0:(len(surrounding_coordinate_2)),4] = surrounding_coordinate_2
		surrounding[0:(len(surrounding_coordinate_3)),5] = surrounding_coordinate_3
	
		### data output
		GroupName = "group1-" + str(GroupNum) + "-inca34-outputs/"
		FrameName = "frame" + str(FrameNum)
		output_groupfolder = os.path.join(Output_moleParentFolder, GroupName)

		## for storing data
		if(not(os.path.exists(output_groupfolder))):
			os.makedirs(output_groupfolder)
		os.chdir(output_groupfolder)
	
		filename = "molecule" + str(MoleculeNum) + "_frame" + str(FrameNum) + ".txt"
		fileid = open(filename, "w")
		fileid.write(" ".join(["intensity1", "intensity2", "intensity3", "coordinate1", "coordinate2", "coordinate3", "\n"]))
		for line in range(0, length):
			line_str = map(str, surrounding[line, :])
			line_str.append("\n")
			fileid.write(" ".join(line_str))
		fileid.close()
		
		print "group", GroupNum, "_frame", FrameNum, "_molecule", MoleculeNum, ": surrounding file created!"
	
		## plotting the frame
		mole_start_position = mfc.get_coordinate(CellNumList_setframe[0],ncol)
		text_position = [mole_start_position[0]-15, mole_start_position[1]]
		text = "mole_" + str(MoleculeNum)
		
		plt.figure(FrameNum, figsize=(16, 12), dpi = 100)
    		mergeplot = plt.imshow(display_array1)
    		mergeplot.set_cmap(colormap)
		plt.text(text_position[1], text_position[0], text, color = "yellow", size = 'medium', weight = 400)
    		plt.colorbar()
    		mergeplot.set_clim([np.min(display_array3), image_max])
	
		framename_pdf = "molecule" + str(MoleculeNum) + "_frame" + str(FrameNum) + ".pdf"
		framename_png = "molecule" + str(MoleculeNum) + "_frame" + str(FrameNum) + ".png"
		plt.savefig(framename_pdf)
		plt.savefig(framename_png)
#		plt.show()
		plt.close()


		## save plots to another folder
		output_groupfolder_plot = os.path.join(Output_plotParentFolder, GroupName)
		if(not(os.path.exists(output_groupfolder_plot))):
			os.makedirs(output_groupfolder_plot)
		framename_pdf_plot = framename_pdf
		framename_png_plot = framename_png
		shutil.copyfile(os.path.join(output_groupfolder,framename_pdf),os.path.join(output_groupfolder_plot,framename_pdf_plot))
		shutil.copyfile(os.path.join(output_groupfolder,framename_png),os.path.join(output_groupfolder_plot,framename_png_plot))
		os.remove(os.path.join(output_groupfolder,framename_pdf))
		os.remove(os.path.join(output_groupfolder,framename_png))


	else:
		print("selected region wrong! no pixel of the molecule in this region..")

#####################################################################################
## func8_1: input GroupNum, FrameNum and MoleculeNum, 
##        get the surrounding intensity(data file and image) of one molecule in one frame
##        close -> True: donut surrounding/ False: sandwich surrounding
##        used specificly for "surrounding_OneFrame"

def surrounding_OneMolecule_forframe(display_array, Input_MoleParentFolder, Output_moleParentFolder, Output_plotParentFolder, GroupNum, FrameNum, MoleculeNum = 0, backbone_range = 2, close = False, molecule_filename = " "):

	GroupNum = int(GroupNum)
	FrameNum = int(FrameNum) 
	MoleculeNum = int(MoleculeNum)

	### data input
	grayData_array = fi.get_grayData_frominf(GroupNum, FrameNum)
	
	if molecule_filename == " ":
		molecule_filename = "molecule" + MoleculeNum + ".txt"
	else:
		MoleculeNum = int(re.findall(r"molecule([0-9]*).txt", molecule_filename)[0])
	
	fileContents = mfc.read_molecule_files_fromname(Input_MoleParentFolder, GroupNum, molecule_filename)
	inframe = mfc.find_mole_setframe(fileContents, FrameNum)
	if(not(inframe)):
		print "the molecule is not in this frame! "
		exit(1)
	
	molecule_contents_setframe = mfc.get_content_setframe(fileContents, FrameNum)

	CellNumList_setframe = mfc.get_CellNumbers(molecule_contents_setframe)

	### surrounding intensity acquisition
	surrounding_intensity_circle_1 = surrounding_intensity_circle(CellNumList_setframe, grayData_array, display_array, 2, 1, close)
	surrounding_intensity_circle_2 = surrounding_intensity_circle(CellNumList_setframe, grayData_array, display_array, 2, 2, close)
	surrounding_intensity_circle_3 = surrounding_intensity_circle(CellNumList_setframe, grayData_array, display_array, 2, 3, close)

	surrounding_coordinate_1 = surrounding_circle(CellNumList_setframe, 2, 1, close)
	surrounding_coordinate_2 = surrounding_circle(CellNumList_setframe, 2, 2, close)
	surrounding_coordinate_3 = surrounding_circle(CellNumList_setframe, 2, 3, close)

	length = len(surrounding_coordinate_3)
	surrounding = np.zeros((length, 6), dtype = int)
	surrounding[0:(len(surrounding_intensity_circle_1)),0] = surrounding_intensity_circle_1
	surrounding[0:(len(surrounding_intensity_circle_2)),1] = surrounding_intensity_circle_2
	surrounding[0:(len(surrounding_intensity_circle_3)),2] = surrounding_intensity_circle_3
	surrounding[0:(len(surrounding_coordinate_1)),3] = surrounding_coordinate_1
	surrounding[0:(len(surrounding_coordinate_2)),4] = surrounding_coordinate_2
	surrounding[0:(len(surrounding_coordinate_3)),5] = surrounding_coordinate_3

	### data output
	GroupName = "group1-" + str(GroupNum) + "-inca34-outputs/"
	FrameName = "frame" + str(FrameNum)

	output_groupfolder = os.path.join(Output_moleParentFolder, GroupName)
	output_framefolder = os.path.join(output_groupfolder, FrameName)

	## for storing data
	if(not(os.path.exists(output_framefolder))):
		os.makedirs(output_framefolder)
	os.chdir(output_framefolder)

	molecule_filename_output = "molecule" + str(MoleculeNum) + "_frame" + str(FrameNum)
	fileid = open(molecule_filename_output, "w")
	fileid.write(" ".join(["intensity1", "intensity2", "intensity3", "coordinate1", "coordinate2", "coordinate3","\n"]))
	for line in range(0, length):
		line_str = map(str, surrounding[line, :])
		line_str.append("\n")
		fileid.write(" ".join(line_str))
	fileid.close()
	
	print "group", GroupNum, "_frame", FrameNum, "_",molecule_filename, ": surrounding file created!"
	
	## creating folder for successive plotting
	output_groupfolder_plot = os.path.join(Output_plotParentFolder, GroupName)
	output_framefolder_plot = os.path.join(output_groupfolder_plot, FrameName)
	if(not(os.path.exists(output_framefolder_plot))):
		os.makedirs(output_framefolder_plot)

	return CellNumList_setframe[0]  ## return the start cell of this molecule


#####################################################################################
## func9: input GroupNum and MoleculeNum, 
##        get the surrounding intensity(data file and image) of one molecule in all frames it stepping across

def surrounding_OneMolecule_allframe(Input_MoleParentFolder, Output_moleParentFolder, Output_plotParentFolder, GroupNum, MoleculeNum, backbone_range = 2, close = False):

	GroupNum = int(GroupNum)
	MoleculeNum = int(MoleculeNum)
	
	fileContents = mfc.read_molecule_files_frominf(Input_MoleParentFolder, GroupNum, MoleculeNum)
	connect_count, frames, cells = mfc.get_frame_count_mole(fileContents)
	for FrameNum in frames:
		surrounding_OneMolecule(Input_MoleParentFolder, Output_moleParentFolder, Output_plotParentFolder, GroupNum, FrameNum, MoleculeNum, [-1,-1], backbone_range, close)


#####################################################################################
## func10: input GroupNum and FrameNum, 
##         get the surrounding intensity(data file and image) of all molecules in this frame

def surrounding_OneFrame(Input_MoleParentFolder, Output_moleParentFolder, Output_plotParentFolder, GroupNum, FrameNum, backbone_range = 2, close = False):
	GroupNum = int(GroupNum)
	FrameNum = int(FrameNum)
	GroupFolder = mfc.get_grouppath(GroupNum, Input_MoleParentFolder)

	
	grayData_array = fi.get_grayData_frominf(GroupNum, FrameNum)
	display_array = copy.deepcopy(grayData_array)

	molecule_filenames = mfc.get_names_setframe(GroupFolder, FrameNum)
	molecule_amount = len(molecule_filenames)
	start_cells = []
	
	if(molecule_amount):
		for molecule_filename in molecule_filenames:
			mole_start_cell = surrounding_OneMolecule_forframe(display_array, Input_MoleParentFolder, Output_moleParentFolder, Output_plotParentFolder, GroupNum, FrameNum, 0, backbone_range, close, molecule_filename)
			start_cells.append(mole_start_cell)
		print "end of frame!"

		plt.figure(FrameNum, figsize=(16, 12), dpi = 100)
    		mergeplot = plt.imshow(display_array)
	    	mergeplot.set_cmap(colormap)
    		plt.colorbar()
		for ind in range(0, molecule_amount):
			mole_start_position = mfc.get_coordinate(start_cells[ind], ncol)
			text_position = [mole_start_position[0]-15, mole_start_position[1]]
			text = molecule_filenames[ind].split(".")[0].split("e")[2]
			plt.text(text_position[1], text_position[0], text, color = "yellow", size = 'medium',weight = 400)
    		mergeplot.set_clim([np.min(display_array), image_max])
	else:
		plt.figure(1,figsize=(16, 12), dpi = 100)
		mergeplot = plt.imshow(display_array)
	    	mergeplot.set_cmap(colormap)
    		plt.colorbar()
		mergeplot.set_clim([np.min(display_array), image_max])

	GroupName = "group1-" + str(GroupNum) + "-inca34-outputs/"
	FrameName = "frame" + str(FrameNum)
	output_groupfolder = os.path.join(Output_moleParentFolder, GroupName)
	output_framefolder = os.path.join(output_groupfolder, FrameName)
	
	if(not os.path.exists(output_framefolder)):
		os.mkdir(output_framefolder)
	os.chdir(output_framefolder)

	framename_png = "frame" + str(FrameNum) + ".png"
	plt.savefig(framename_png)
	output_groupfolder_plot = os.path.join(Output_plotParentFolder, GroupName)
	output_framefolder_plot = os.path.join(output_groupfolder_plot, FrameName)
	if(not os.path.exists(output_framefolder_plot)):
		os.mkdir(output_framefolder_plot)
	
	framename_png_plot = framename_png
	shutil.copyfile(os.path.join(output_framefolder,framename_png),os.path.join(output_framefolder_plot,framename_png_plot))
	plt.close()
#	plt.show()

	return molecule_amount
        

#####################################################################################
## func11: input GroupNum,
##         get the surrounding intensity(data file and image) of all molecules in this group

def surrounding_OneGroup(Input_MoleParentFolder, Output_moleParentFolder, Output_plotParentFolder, GroupNum, backbone_range = 2, close = False):
	GroupNum = int(GroupNum)
	GroupFolder = mfc.get_grouppath(GroupNum, Input_MoleParentFolder)

	framenames = fi.get_framenames(GroupNum)
	for FrameNum in range(0, len(framenames)):
			print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ for frame ", FrameNum, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			surrounding_OneFrame(Input_MoleParentFolder, Output_moleParentFolder, Output_plotParentFolder, GroupNum, FrameNum, backbone_range, close)
	print "end of group!"




######################################################################################
## func_end: return name of this model

def modelName():
	return __name__






	
		












		
