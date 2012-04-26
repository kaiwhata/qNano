#!/usr/bin/python

import os
import xml.etree.ElementTree as xml
import fnmatch
import re, string, sys 


def create_library(parent, pattern):
	datafile_library = []
	for path, dirs, files in os.walk(os.path.abspath(parent)):
		for filename in fnmatch.filter(files, pattern):
			datafile_library.append(os.path.join(path, filename))
	print datafile_library
	return datafile_library

def correct_name(name):
	######find pH and particle charge from filename#####
	sep = re.split('/',name)
	
	trial_name = sep[-3]
	date = sep[-5]
	folder = sep[-4]
	voltage = sep[-2]
	print trial_name, date, folder, voltage
	return trial_name, date, folder, voltage

def main():
	parent = "/media/ADATA CH94/PhD Data/9.2011/"
	pattern = "*.xml"
	datafile_library = create_library(parent, pattern)
	
	output=open("/home/jeldridge/expts_summary.csv","a")  #save results to output file
	output.write('%s , %s , %s , %s , %s , %s , %s, %s\n'%('Trial Name','Date','Folder','Voltage','Experiment name','Pore','Pressure', 'Full Trial Name')) #print the header
	
	for trial in datafile_library:
		xml_file = trial
		trial_name = correct_name(trial)
		#xml_file = os.path.join("/home/jeldridge/", "sample.xml")
		try:
			tree = xml.parse(xml_file)
			rootElement = tree.getroot()
			experimentElement = tree.findtext("Experiment")
			poreElement = tree.findtext("Aperture")
			pressureElement = tree.findtext("Pressure")
			#poreElement = tree.findtext("ElectrolyteID")
			output.write('%s , %s, %s, %s, %s , %s , %s, %s\n'%(trial_name[0],trial_name[1],trial_name[2],trial_name[3],experimentElement, poreElement, pressureElement, trial )) #print data
		except Exception, inst:
			print "Unexpected error opening %s: %s" % (xml_file, inst) 
			return
			
	output.close
			
		
	

if __name__ == '__main__': 
	main() 

