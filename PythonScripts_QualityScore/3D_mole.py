#! /home/cwu269/cwu/program/bin


import getopt
import io
import numpy as np
import pickle
import shutil
import struct
import copy
import re
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import lib_frames_index as fi  			# model1: get access to frames in database
import lib_molecule_file_control as mfc 	# model2: operation about molecule files
import lib_get_mole_surrounding as gms

###################################################################################
## FUNCTION: plot 3D-surface figure for pixels within and around one molecule backbone

## python 3D_mole.py GroupNum MoleculeNum

## constants: shape of frame

ncol = 1392
nrow = 1024

###################################################################################################
## global variables

Input_MoleParentFolder = "/aspen/nandi/MF_cap348/maps_inca34"
Input_ListFolder = "/exports/aspen/cwu/Project_QualityScore/Data/molecule_file"
#Output_plotParentFolder = "/exports/home/cwu269/cwu/Project_QualityScore/Plots/MF"

x_range = 7

###################################################################################################
def diffusion_position(CellNum):
	diff_CellNumList = []
	diff_cell_condition_pre = gms.target_pixel(CellNum, 0, 0)
	for d in range(-1, -x_range-1, -1):
		diff_cell_condition = gms.target_pixel(CellNum, d, 0)
		if(diff_cell_condition[1]):
			diff_CellNumList.append(diff_cell_condition[0])
		else:
			diff_CellNumList.append(diff_cell_condition_pre[0])
		diff_cell_condition_pre = diff_cell_condition

		diff_cell_condition_pre = gms.target_pixel(CellNum, 0, 0)

	diff_cell_condition_pre = gms.target_pixel(CellNum, 0, 0)
	diff_CellNumList.reverse()
	for d in range(0, x_range+1, 1):
		diff_cell_condition = gms.target_pixel(CellNum, d, 0)
		if(diff_cell_condition[1]):
			diff_CellNumList.append(diff_cell_condition[0])
		else:
			diff_CellNumList.append(diff_cell_condition_pre[0])
		diff_cell_condition_pre = diff_cell_condition
	
	return diff_CellNumList 

###################################################################################################
def diffusion_intensity(CellNum, grayData_array):

	diff_intensity = []
	diff_CellNumList = diffusion_position(CellNum)

	for diff_cellnum in diff_CellNumList:
	
		diff_cell_intensity = gms.get_pixel_intensity(diff_cellnum, grayData_array)
		diff_intensity.append(diff_cell_intensity)

	return diff_intensity

		
###################################################################################################
def main(argv):

	try:                                
        	opts, args = getopt.getopt(argv, "hs:e:", ["help", "grammar="])
    	except getopt.GetoptError:          
        	usage()                         
        	sys.exit(2)   

	GroupNum = int(args[0])
	MoleculeNum = int(args[1])
	
	GroupName = "group1-" + str(GroupNum) + "-inca34-outputs/"
	MoleculeName = "molecule" + str(MoleculeNum) + ".txt"
	molecule_file = os.path.join(Input_MoleParentFolder, GroupName, MoleculeName)

	if os.path.exists(molecule_file):
		
		fileContents = mfc.read_molecule_files(molecule_file)
		numPixels = mfc.get_numPixels(fileContents)
		CellNumList = mfc.get_CellNumbers(fileContents)
		FrameNumList = mfc.get_FrameNumbers(fileContents)			

		frames = list(set(FrameNumList))
		grayDataList = []
		for FrameNum in frames:
			grayDataList.append(fi.get_grayData_frominf(GroupNum, FrameNum))

		z = []
		for i in range(0, numPixels):
			FrameNumIndex = frames.index(FrameNumList[i])
			row = diffusion_intensity(CellNumList[i], grayDataList[FrameNumIndex])
			z.append(row)

		z = np.array(z)
		x = np.arange(-x_range,(x_range + 1) ,1)
		y = np.arange(0, numPixels, 1)
		x, y = np.meshgrid(x,y)

#		PlotPath = os.path.join(Output_plotParentFolder, IntervalName, GroupName)
#		if(not(os.path.exists(PlotPath))):
#			os.makedirs(PlotPath)
#		os.chdir(PlotPath) 
		   	
		fig = plt.figure()
		ax = fig.gca(projection = '3d')
		
		surf = ax.plot_surface(x, y, z, rstride = 1, cstride = 1, cmap = cm.coolwarm, linewidth = 0, antialiased = False)
#		surf = ax.plot_wireframe(x, y, z, rstride = 1, cstride = 1, cmap = cm.coolwarm, linewidth = 1, antialiased = False)
#		ax.set_zlim()
		ax.zaxis.set_major_locator(LinearLocator(10))
#		ax.zaxis.set_majot_formatter(FormatStrFormatter('%.02f'))
		
		Title = "Group" + str(GroupNum) + "_Molecule" + str(MoleculeNum) + "_intensity_distribution"
		plt.title(Title)
		ax.set_xlabel('diffusion_range')
		ax.set_ylabel('length')
		ax.set_zlabel('intensity')
		fig.colorbar(surf, shrink = 0.5)
#		framename_pdf = "mole" + str(MoleculeNum) + "_backbone1_region" + ".pdf"
			
#		plt.savefig()
#		plt.close()
					
		plt.show()
	else:
		print mole, "not exist!"

############################################################################################
## calling main function

if __name__ == "__main__":
    main(sys.argv[1:])


