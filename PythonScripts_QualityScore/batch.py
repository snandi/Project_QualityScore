#! /home/nandi/anaconda/bin/python

import os

intervals = range(0, 39)

def main():

	for IntervalNum in intervals:
		os.system("python produce_molelist.py " + str(IntervalNum)) 
		os.system("python get_molefiles_region.py " + str(IntervalNum))
		os.system("python mole_surroundings_mole_list.py " + str(IntervalNum))
		os.system("python 5to1.py " + str(IntervalNum))
		os.system("python frame_merge_dp_mole_list.py " + str(IntervalNum))



if __name__ == "__main__":
	main()
