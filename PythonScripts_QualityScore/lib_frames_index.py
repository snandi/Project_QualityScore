#! /home/nandi/anaconda/bin/python

import getopt
import gzip
import io
import MySQLdb
import numpy
import os
import pickle
import shutil
import struct
import sys

from string import translate, maketrans, punctuation 
from os import path

###################################################################################
## model1: get access to frames in database

version = '1.0'

## parameters of frame: set for convenience
col_num = 1392
row_num = 1024

###################################################################################
## func1: get the fullpath of one group in database, inputting group number and operation system type

def get_group_path(anotaggrpnum, ostype):
    
    # Open database connection
    #db = MySQLdb.connect(host="peach.lmcg.wisc.edu", user="test", passwd="9654", db="ommdb" )
    db = MySQLdb.connect(host="db01.lmcg.wisc.edu", user="o_user", passwd="o_pass", db="ommdb")
    
    # prepare a cursor object
    cursor = db.cursor()
    
    # Set up the SQL Query
    sql = "SELECT FullPath FROM Groups where GroupID = %s " % (anotaggrpnum)

    fullpath = None
    try:
        # execute the query
        cursor.execute(sql)
        # Fetch the result
        fullpath = cursor.fetchall()[0][0]
    except:
        print "Error: Unable to Run Query"
    
    # Disconnect from server
    db.close()
    
    if fullpath and "win" in ostype.lower():
        #if windows then we need to modify
        winfullpath = "".join(["//PEPPER/linuxfs", fullpath])
        winchangechars = "/"
        T = maketrans(winchangechars, '\\'*len(winchangechars))
        fullpath = translate(winfullpath, T)
        
    return fullpath

###################################################################################
## func2: the subfolder of one group containing image data 

def get_group_imagedir_path(groupommpath):
    return os.path.join(groupommpath, 'run0', 'flat')

###################################################################################
## func3: get a list storing all frame names in one group

def get_image_filenames(imagedir):
    ## Written by Nandi, same as the get_group_image_filenames function, but without the backfilenames
    #backfilenames and puncfilenames are list of strings
    backfilenames = []
    
    try:
        #get sorted filenames, filenames must contain keyword correct
        file_list = [f for f in sorted(os.listdir(imagedir)) if 'correct' in f]
        
        #split it into back and punc filenames
        #back filenames are those with indices [0:end:3]
        for i in range(0, len(file_list), 3):
            backfilenames.append(os.path.join(imagedir, file_list[i]))
            
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror), ':', imagedir
        exit(-1)
                    
    return backfilenames


###################################################################################
## func4: get the intensity gray data from one frame, inputting the fullpath&name of the frame file
def get_grayData_frompath(filename):

    if os.path.exists(filename) and os.path.isfile(filename):
    	try:
    	    #read the binary data into image_data
    	    omi_file = gzip.open(filename, 'rb')
    	    image_data = omi_file.read()
    	    omi_file.close()
    	    
            #read the header data
            header_data = [struct.unpack(">i", image_data[i:i+4])[0] for i in range(0, 19, 4)] 
                
            #allocate a required size array
            num_bytes = header_data[4]
            num_pixels = header_data[2]*header_data[3] 	
            gray_data = numpy.zeros(num_pixels, dtype=int)
    	    		
            #convert binary data into intensities in gray_data
            start_pos = 20
            stop_pos = start_pos + num_bytes
            for i in range(num_pixels):
                gray_data[i] = struct.unpack(">H", image_data[start_pos:stop_pos])[0]
                start_pos = stop_pos
                stop_pos = start_pos + num_bytes
    	    gray_data_reshaped = gray_data.reshape(header_data[2], header_data[3])  	  
                
            return gray_data_reshaped            
        except:
            print 'not able to read file'
    else:
        print filename, 'is invalid!'
        sys.exit(-1)

###################################################################################
## func5: get framenames in one group

def get_framenames(GroupNum):
    groupOmmPath = get_group_path(GroupNum, "lin")
    imageDir = get_group_imagedir_path(groupommpath=groupOmmPath)
    filenames = get_image_filenames(imagedir=imageDir)
    return filenames



###################################################################################
## summary for getting gray array

def get_grayData_frominf(GroupNum, FrameNum):
    groupOmmPath = get_group_path(GroupNum, "lin")
    imageDir = get_group_imagedir_path(groupommpath=groupOmmPath)
    filenames = get_image_filenames(imagedir=imageDir)
    #print len(filenames)
    grayData = get_grayData_frompath(filenames[int(FrameNum)])
    return grayData



###################################################################################
## func_end: return name of this model

def modelName():
	return __name__





