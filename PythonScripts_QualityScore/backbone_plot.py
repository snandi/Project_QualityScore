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


###################################################################################
## FUNCTION: plot the intensity profiles from 1-pixel backbones for the molecule aligned to one interval

## $python backbone_plot.py IntervalNum



###################################################################################################
## global variables

Input_MoleParentFolder = "/exports/aspen/cwu/Project_QualityScore/Data/molecule_file/1pixel_region"
Input_ListFolder = "/exports/aspen/cwu/Project_QualityScore/Data/molecule_file"
Output_plotParentFolder = "/exports/home/cwu269/cwu/Project_QualityScore/Plots/MF"

###################################################################################################
def main(argv):

	try:                                
        	opts, args = getopt.getopt(argv, "hs:e:", ["help", "grammar="])
    	except getopt.GetoptError:          
        	usage()                         
        	sys.exit(2)   
   
	IntervalNum = args[0]
#	ChrNum = int(argv[1])
	ChrNum = 1
	
	ListFileName = 'Chr' + str(ChrNum) + '_interval' + str(IntervalNum) + '.txt'
	FileReadPath = Input_ListFolder + "/" + ListFileName
	fileid = open(FileReadPath, 'r')
	mole_list = fileid.read()
	mole_list = mole_list.split(' ')
	

	for mole in mole_list:
		GroupNum = int(mole[0:7])
		MoleculeNum = int(mole[15:19])
		
		IntervalName = "refFrag_" + str(IntervalNum)
	 	GroupName = "group1-" + str(GroupNum) + "-inca34-outputs/"
		MoleculeName = "molecule" + str(MoleculeNum) + ".txt"
		molecule_file = os.path.join(Input_MoleParentFolder, GroupName, MoleculeName)

		if os.path.exists(molecule_file):
		
    			f = open(molecule_file, "r")
    			data = f.readlines()
    			i = 0
    			for line in data:
    				linenum = re.split(' ', line)
    				if(len(linenum) == 4):
					if(i == 0): 
						arr = [int(linenum[2])]
					else:
						arr.append(int(linenum[2]))
					i = i + 1
		
			PlotPath = os.path.join(Output_plotParentFolder, IntervalName, GroupName)
			if(not(os.path.exists(PlotPath))):
				os.makedirs(PlotPath)
    			os.chdir(PlotPath) 
			   	
			plt.figure(1, figsize=(16,12), dpi=100)
			plt.plot(arr)
		
			Title = "Group" + str(GroupNum) + "_Molecule" + str(MoleculeNum) + " Backbone1 intensity"
			plt.title(Title)
			plt.xlabel('position(left to right)')
			plt.ylabel('intensity')
			
	
			framename_pdf = "mole" + str(MoleculeNum) + "_backbone1_region" + ".pdf"
			
			plt.savefig(framename_pdf)
			plt.close()
		
			
#			plt.show()
		else:
			print mole, "not exist!"


############################################################################################
## calling main function

if __name__ == "__main__":
    main(sys.argv[1:])
