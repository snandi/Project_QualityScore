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

#import lib_frames_index as fi  		# model1: get access to frames in database
#import lib_molecule_file_control as mfc 	# model2: operation about molecule files


#######################################################################################
## FUNCTION: get the molecule list file given interval number (and Chromosome number)

## $python get_molelist.py IntervalNumber 

#######################################################################################
## global variables

Input_inforParentFolder = '/exports/aspen/cwu/Project_QualityScore/Data/molecule_file/'
Output_ListParentFolder = '/exports/aspen/cwu/Project_QualityScore/Data/molecule_file/'
InforFileName = 'MF_cap348_inca34_cf209_minSize50_minFrag5_Aligned.alignmentChunks'

#######################################################################################
## main: get the molecule list file given interval number (and Chromosome number)

def main(argv):

	try:                                
        	opts, args = getopt.getopt(argv, "hg:d", ["help", "grammar="])
    	except getopt.GetoptError:          
        	usage()                         
        	sys.exit(2)      
    
    	IntervalNum = int(args[0])
#	ChrNum = int(args[1]) 
	ChrNum = 1

	ListFileName = 'Chr' + str(ChrNum) + '_interval' + str(IntervalNum) + '.txt'
	InforReadPath = Input_inforParentFolder + InforFileName
	fp = open(InforReadPath, "r")
	content = fp.read()
	fp.close()

	regex = 'chr' + str(ChrNum) + '\s' + str(IntervalNum) +  '\s' + str(IntervalNum + 1) + '\s' + '(.{19})' + '\s[0-9]+\s[0-9]+\s[0-9]+\s[0-9]+\s[0-9]*\s[0-9]*\s.+?\s.*'
	IntervalInfor = re.findall(regex, content)

	ListWritePath = Output_ListParentFolder + ListFileName
	lp = open(ListWritePath, "w")
	lp.seek(0)
	lp.write(" ".join(IntervalInfor))
	lp.close()


#######################################################################################
if __name__=="__main__":
	main(sys.argv[1:])
