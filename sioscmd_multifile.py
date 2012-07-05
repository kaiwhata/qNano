#!/usr/bin/python

import re, os
import fnmatch
from time import time

def create_library(parent, pattern):
	datafile_library = []
	for path, dirs, files in os.walk(os.path.abspath(parent)):
		for filename in fnmatch.filter(files, pattern):
			datafile_library.append(os.path.join(path, filename))
	print datafile_library
	return datafile_library
	
def main():
	parent = 'C:\\Izon_Data\\Electrolyte_Concentration_june2012\\Electrolyte_Runs\\50mM\\'
	pattern = "*.irbd"
	datafile_library = create_library(parent, pattern)
	print "I have found "+str(len(datafile_library))+" files to process"
	print "This is likely to take around "+str(2*len(datafile_library))+" minutes"
	print "I would recommend going for tea or perhaps checking facebook"
	i = 0
	for trial in datafile_library:
		time_i = time()
		os.system("sioscmd -i "+str(trial)+" -csvBlockadeFile")
		time_f = time()
		timer = round(((time_f-time_i)/60.0), 2)
		print "This run took "+str(timer)+" minutes"
		print "I have about "+str(((len(datafile_library)-(i+1))*timer))+"mins left to run"
		i=i+1	

if __name__ == '__main__': 
	main() 


