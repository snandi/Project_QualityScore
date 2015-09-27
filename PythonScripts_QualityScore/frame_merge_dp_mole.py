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
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from string import translate, maketrans, punctuation 
from os import path

import lib_frames_index as fi  			# model1: get access to frames in database
import lib_molecule_file_control as mfc 	# model2: operation about molecule files

###################################################################################
## FUNCTION:
## given the group number and molecule number
## merging frames this molecule stepping across, and labeling this molecule on the image


###################################################################################
## global variables
 
Input_MoleculeParentFolder = "/aspen/nandi/MF_cap348/maps_inca34"
Output_plotParentFolder = "/exports/home/cwu269/cwu/Project_QualityScore/Plots/surrounding"

nRow = 1024 #fron_framedata.shape[0]
nCol = 1392 #fron_framedata.shape[1]

###################################################################################
## Main                                                                          ##

def main(argv):

    	try:                                
        	opts, args = getopt.getopt(argv, "hg:d", ["help", "grammar="])
    	except getopt.GetoptError:          
        	usage()                         
        	sys.exit(2)      
    
    	GroupNum = str(args[0])
    	MoleculeNum = str(args[1])
    
    	## for molecule information
    	molecule_filepath = mfc.get_molefilepath(GroupNum, MoleculeNum, Input_MoleculeParentFolder)
    	mole_content = mfc.read_molecule_files(molecule_filepath)

    	Line1 = mole_content[1].split(" ")
    	LineEnd = mole_content[int(mole_content[0].split(" ")[0])].split(" ")
    	start_cell = int(Line1[1])
    	end_cell = int(LineEnd[1])

    	connect, connect_frames, connect_cells = mfc.get_frame_count_mole(mole_content)
   
    	connect_cells = np.array(connect_cells)

    	numPixels = mfc.get_numPixels(mole_content)
    	skipmolecule = mfc.checkforMinusOne(mole_content,numPixels)

    	dif = np.zeros((connect+2,2),dtype=int)
    	dif[0, :] = [0, 0]
    	dif[connect+1, :] = [0, 0]


    	## merging!
    	if(skipmolecule):
		print "broken! skip..."
    	else:
		
    		if(connect): 
			print "connection!"
			fron_framedata = fi.get_grayData_frominf(GroupNum, connect_frames[0])
	
			for i in range(1,connect+1):
			
				fron_connect_point = np.array(mfc.get_coordinate(connect_cells[i-1,0], nCol))
    				hind_connect_point = np.array(mfc.get_coordinate(connect_cells[i-1,1], nCol))
				dif[i,:] = hind_connect_point - fron_connect_point
					
			rows = nRow + abs(dif.sum(axis=0)[0])
			cols = nCol + abs(dif.sum(axis=0)[1])
	
			mergedata = np.zeros((rows, cols),dtype=int)
	
			if(dif[1,0] >= 0):
				for i in range(0,connect+1):
					framedata = fi.get_grayData_frominf(GroupNum, connect_frames[i])
					current_row = abs(sum(dif[(i+1):,0]))
					current_col = abs(sum(dif[0:(i+1),1]))
					mergedata[current_row:(current_row + nRow), current_col:(current_col + nCol)] = framedata
			if(dif[1,0] < 0):
				for i in range(0,connect+1):
					framedata = fi.get_grayData_frominf(GroupNum, connect_frames[i])
					current_row = abs(sum(dif[0:(i+1),0]))
					current_col = abs(sum(dif[0:(i+1),1]))
					mergedata[current_row:(current_row + nRow), current_col:(current_col + nCol)] = framedata
    		else:	
        		print "not connection!"
			mergedata = fi.get_grayData_frominf(GroupNum, connect_frames[0])		
			
		## plotting the merged image
		mole_start_position = mfc.get_coordinate(start_cell, nCol)
		mole_end_position = mfc.get_coordinate(end_cell, nCol)
		mole_text1_position = [mole_start_position[1], mole_start_position[0] + abs(sum(dif[:,0]))]
		mole_text2_position = [mole_end_position[1] + abs(sum(dif[:,1])), mole_end_position[0]]

		## showing
		plt.figure(connect+1, figsize = (16*(connect+1), 12), dpi = 100)
    		mergeplot = plt.imshow(mergedata)
    		mergeplot.set_cmap('gist_ncar')
		text1 = MoleculeNum + "_s"
		text2 = MoleculeNum + "_e"
		plt.text(mole_text1_position[0],mole_text1_position[1], text1, size = 'xx-large', color = "red", weight = 700)
		plt.text(mole_text2_position[0],mole_text2_position[1], text2, size = 'xx-large', color = "red", weight = 700)
    		plt.colorbar()
    		mergeplot.set_clim([np.min(mergedata), 30000])

		figname_png = "group" + str(GroupNum) + "_molecule" + str(MoleculeNum)
		os.chdir(Output_plotParentFolder)
		plt.savefig(figname_png)
		plt.show()



#################################################################################
## calling main function   

if __name__ == "__main__":
  	main(sys.argv[1:])

