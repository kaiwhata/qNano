#!/usr/bin/python

#locator.py - locates and returns a list of all 'blockade_1.csv' 
#type files in a directory, copies and renames them with their
#parent directory's name and sticks them in the parent directory.

import re, os
import fnmatch
import csv
import math

def calibrate(): #calculates the pressure for each event
	r = (46.5/1000)  #radius of cylinder in m
	p = (200.00/1000.00) #in litres of H2O
	area = math.pi*r*r
	initial_pressure = (p/area) #in mmH2O
	print "Initial applied pressure ="
	print initial_pressure


def load(): #locates all files with the pattern and returns a list of their locations
	parent = "/media/Expansion Drive/"
	pattern = "blockades*.csv"
	datafile_library = []

	for path, dirs, files in os.walk(os.path.abspath(parent)):
		for filename in fnmatch.filter(files, pattern):
			print os.path.join(path, filename) 
			datafile_library.append(os.path.join(path, filename))
	
	print datafile_library

def csv_append():
	for trial in datafile_library:
		file = trial
		print trial		
		ifile = open(file,"rb")
		###find name from directory##
		sep = re.split('/',trial)
		correct_name = sep[(len(sep)-3)]
		output_filename = parent+correct_name+".csv"
		#####################
		ofile = open(output_filename, "wb") #output_filename
		reader = csv.reader(ifile)
		writer = csv.writer(ofile, delimiter=',')# ,quotechar='"',quoting=csv.QUOTE_ALL)
		rownum = 0
		for row in reader:
			if rownum == 0: #writing header row in output.csv
				minlength = len(row)
				header =row
				header.append('cumulative_counts')
				#header.append('mass_flow')
				#header.append('total_pressure')
				#header.append('applied_pressure')
				print header
			elif len(row) < minlength: #ensures all included rows contain all data and discards any that dont
				print str(rownum)+" omitted due to incomplete data"
			else: #calculates rapplied pressure based no time and calibration flow rate
				#time_min = float(row[0])/60
				#mass_flow = (-0.058*time_min*time_min+30.142*time_min-0.6686)
				#total_pressure = (mass_flow/1000)/(math.pi*r*r)
				#applied_pressure = initial_pressure-total_pressure			
				row.append(rownum)
				#row.append(mass_flow)
				#row.append(total_pressure)
				#row.append(applied_pressure)			
				writer.writerow(row)
			rownum += 1
		ifile.close()
		ofile.close()

def main():
	calibrate() #may not be necessary for non-pressure change runs
	load()
	csv_append()

if __name__ == '__main__': 
	main() 
