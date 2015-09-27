#! /home/nandi/anaconda/bin/python

import getopt
import io
import numpy as np
import os
import pickle
import shutil
import struct
import sys
import math

from os import path

import lib_frames_index as fi  			# model1: get access to frames in database
import lib_molecule_file_control as mfc 	# model2: operation about molecule files
import lib_get_mole_surrounding as gms 		# model3: get the surrounding pixels(donut or sandwich) around molecule backbone

###################################################################################
## model5: get the straightness measure of molecule fragment with simplify algorithm


##############################################
## global variables

#threshold = 3

ncol = 1392
nrow = 1024

###################################################################################
## func0_1: multiply, used for map

def product(a, b):
	return a * b

###################################################################################
## func0_1: get measure of angle(-pi to pi) from sin & cos

def measure_angle(sin, cos):

	angle = math.acos(cos)
	if sin < 0:
		angle = - angle
	
	return angle

###################################################################################
## func1: get the distance between a point and a line

def point_to_line_distance(A,B,C, x0, y0):

	d = float(abs(A*x0 + B*y0 + C) / math.sqrt(A*A + B*B))
	
	return d

###################################################################################
## func2_1: get consine of the angle between two vectors

def cosine_distance(vector1, vector2):

	vector1 = np.array(vector1)
	vector2 = np.array(vector2)

	inner_product = np.dot(vector1, vector2)
	
	length1 = math.sqrt(np.dot(vector1, vector1))
	length2 = math.sqrt(np.dot(vector2, vector2))
	
	cos = float(inner_product / (length1 * length2))
	
	return cos

###################################################################################
## func2_2: get consine of the angle between two vectors
def sine_distance(vector1, vector2):

	vector1 = np.array(vector1)
	vector2 = np.array(vector2)

	outer_product = np.cross(vector1, vector2)

	length1 = math.sqrt(np.dot(vector1, vector1))
	length2 = math.sqrt(np.dot(vector2, vector2))

	sin = float(outer_product / (length1 * length2))

	return sin

###################################################################################
## func3: get the max distance of a curve to a line

def max_distance(celllist):

	cell_head = celllist[0]
	cell_tail = celllist[len(celllist) - 1]
	[y1, x1] = mfc.get_coordinate(cell_head, ncol)
	[y2, x2] = mfc.get_coordinate(cell_tail, ncol)

	A = y1-y2
	B = x2-x1
	C = x1*y2 - y1*x2

	coordinates = map(mfc.get_coordinate, celllist, list(np.repeat(ncol,len(celllist))))
	coordinates = np.array(coordinates)
	X = list(coordinates[:,1])
	Y = list(coordinates[:,0])
	
	ds = map(point_to_line_distance, list(np.repeat(A,len(celllist))),list(np.repeat(B,len(celllist))), list(np.repeat(C,len(celllist))), X, Y)

	max_dis = max(ds)
	
	mid = ds.index(max_dis)

	return [max_dis, mid]	



###################################################################################
## func4: get the simplify lines of a curve(set recursion depth limit as 900)

def simplify_line(lst, low, high, list_knots, threshold, re_time):
	if low <= high and re_time[0] < 900:
		re_time[0] += 1
		[max_dis, mid] = max_distance(lst[low:(high+1)])
		mid += low
		if max_dis > threshold:
			list_knots.append(mid)
			simplify_line(lst, low, mid, list_knots, threshold, re_time)
			simplify_line(lst, mid, high, list_knots, threshold, re_time)
		else:
			return 1
	else:
		return 1
###################################################################################
## func5: check whether the simplification is finished

def simplify_line_test(lst_test, threshold):
	
	[max_dis, mid] = max_distance(lst_test)
	if max_dis > threshold:
		return 1  ## haven't finished
	else:
		return 0  ## finish!
	

		
###################################################################################
## func6: get cosine values of turning angles for a polyline

def angle_list(turning_knots_cellnums):
	
	if len(turning_knots_cellnums) == 2:
		cos_list = [1]
		sin_list = [0]
	else:
		cos_list = []
		sin_list = []
		for i in range(0, len(turning_knots_cellnums) - 2):
			[y1, x1] = mfc.get_coordinate(turning_knots_cellnums[i], ncol)
			[y2, x2] = mfc.get_coordinate(turning_knots_cellnums[i + 1], ncol)
			[y3, x3] = mfc.get_coordinate(turning_knots_cellnums[i + 2], ncol)
	
			vector1 = [x2-x1, y2-y1]
			vector2 = [x3-x2, y3-y2]
	
			cos_list.append(cosine_distance(vector1, vector2))
			sin_list.append(sine_distance(vector1, vector2))
	theta_list = map(measure_angle, sin_list, cos_list)
	return theta_list



###################################################################################
## func7: get all the simplify lines for one molecule

def simplify_lines_AllFrame(Input_MoleParentFolder_interval, GroupNum, MoleculeNum, threshold):

	fileContents = mfc.read_molecule_files_frominf(Input_MoleParentFolder_interval, GroupNum, MoleculeNum)
	connect_count, frames, cells = mfc.get_frame_count_mole(fileContents)

	straightness = []

	for FrameNum in frames:
		

		## readin molecule backbone information
		fileContents_setframe = mfc.get_content_setframe(fileContents, FrameNum)
		CellNumList_setframe = mfc.get_CellNumbers(fileContents_setframe)

		## create a list to store turning points of simplify polyline
		list_knots = [0, len(CellNumList_setframe)-1]

		## generate the simplify line
		simplify_test = simplify_line_test(CellNumList_setframe, threshold)
		while(simplify_test):

			## for first recursion
			if len(list_knots) == 2:
				re_time = [0]
				simplify_line(CellNumList_setframe, 0, len(CellNumList_setframe)-1, list_knots, threshold, re_time)
			## downstream recursion
			else:
				## for each segment
				for i in range(0, len(list_knots) - 2):
					## check whether this segment need more recursion
					if simplify_test_list[i]:
						re_time = [0]
						simplify_line(CellNumList_setframe, list_knots[i], list_knots[i+1], list_knots, threshold, re_time)
			list_knots = sorted(list(set(list_knots)))

			## check every segment, decide simplify polyline creation is finished or not
			simplify_test_list = []
			for i in range(0, len(list_knots) - 2):
				simplify_test_list.append(simplify_line_test(CellNumList_setframe[list_knots[i]:(list_knots[i+1] + 1)], threshold))
			simplify_test = max(simplify_test_list)
		
		## get turning positions
		turning_knots_cellnums = list(np.array(CellNumList_setframe)[list_knots])

		## get information related to turning angles
		theta_list = angle_list(turning_knots_cellnums)
		theta_max = max(map(abs, theta_list))
		theta_sum = sum(map(abs, theta_list))
		theta_net = sum(theta_list)
	
		## estimation
		test_score_sum = float(1 / math.exp((abs(theta_net) + theta_max) / 2))

		## result report
		straightness.append([FrameNum, theta_net, theta_max, test_score_sum])

	return straightness












