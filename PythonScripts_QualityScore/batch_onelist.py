
import os
import sys
import getopt

def main(argv):

    	try:                                
        	opts, args = getopt.getopt(argv, "hg:d", ["help", "grammar="])
    	except getopt.GetoptError:          
        	usage()                         
        	sys.exit(2)     

	IntervalNum = int(args[0])

	os.system("python produce_molelist.py " + str(IntervalNum)) 
	os.system("python get_molefiles_region.py " + str(IntervalNum))
	os.system("python mole_surroundings_mole_list.py " + str(IntervalNum))
	os.system("python 5to1.py " + str(IntervalNum))
	os.system("python frame_merge_dp_mole_list.py " + str(IntervalNum))



if __name__ == "__main__":
	main(sys.argv[1:])
