#!/usr/bin/python

import re, string, sys 
import os
import numpy
import pylab
from pylab import *
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.stats import mode
import scipy 
import fnmatch



def create_library(parent, pattern):
	datafile_library = []
	for path, dirs, files in os.walk(os.path.abspath(parent)):
		for filename in fnmatch.filter(files, pattern):
			datafile_library.append(os.path.join(path, filename))
	print datafile_library
	return datafile_library
	
def splitting(data,number):
	#creating averaging variables
	average_number = number
	data_list = data.tolist()#convert arrays to indexed lists
	grouped_data_slices = zip(*(iter(data_list),) * average_number)#split app_pressure_list of tuples of length 'average_number'
	av_data_list = [] #define emptylists
	#calculates the average from each group of tuples and places it in the list 'av_app_pressure'
	for i in grouped_data_slices:
		values = i
		av_data_list.append(numpy.average(values))
	return av_data_list
	
def correct_name(name):
	######find pH and particle charge from filename#####
	sep = re.split('/',name)
	trial_name = sep[-1]
	full_name = re.split('_',sep[-1])
	particle_charge_string = full_name[0]
	particle_charge = int(particle_charge_string[1:])
	pH_string = full_name[1]
	pH = str(pH_string[2:])
	return trial_name, particle_charge, pH, 
	
def main(): 
	parent = "/home/jeldridge/current_data/Elf_variation_test/"
	pattern = "B121_pH8_0.5*.csv"
	datafile_library = create_library(parent, pattern)
	trial_library = []
	#creating averaging variables
	FWHM_average_number = 5
	blockade_average_number = 50
	
	for trial in datafile_library:
		file = trial
		print trial	
		csvdata = numpy.loadtxt(file, delimiter=',', usecols = (3,5,7,8,14,17)) #loading data #3=FWHM, 5=FWQM,8 = baseline,14 = cumulative counts, 15 = mass outflow,16 = height change,17= app pressure
		
		cumulative_counts = csvdata[:,4]
		app_pressure = csvdata[:,5]
		FWHM = csvdata[:,0]
		FWQM = csvdata[:,1]
		blockade_magnitude = csvdata[:,2]
		trace = csvdata[:,3]
			
		#calculates mean averages (only use as indication)
		average_blockade_magnitude = sum(blockade_magnitude)/len(blockade_magnitude)
		average_trace = sum(trace)/len(trace)
		#calculates deltaI/I for each event
		normalised_blockades = 	blockade_magnitude/trace
		
		##########################setting up averages for FWHM/FWQM plot#################
		av_norm_pressure = splitting(app_pressure,blockade_average_number)
		av_app_pressure = splitting(app_pressure,FWHM_average_number)
		av_blockade_magnitude = splitting(blockade_magnitude,blockade_average_number)
		av_FWHM = splitting(FWHM,FWHM_average_number)
		av_normalised_blockades = splitting(normalised_blockades,blockade_average_number)
				
		plot(app_pressure, trace)
		plt.xlabel('Average Applied Pressure (mmH20)')
		plt.ylabel('Baseline current (nA)')
			
		trial_name = correct_name(trial)[0]
		particle_charge = correct_name(trial)[1]
		pH = correct_name(trial)[2]	
		
		trial_library.append(trial_name)
				
		geometry_indicator = (average_blockade_magnitude/average_trace)
		geometry_text = '${\Delta I}/{I}$ = '+'%.3e' % (geometry_indicator) 
						
		plt.text(min(app_pressure), min(trace), geometry_text)
				
		#plt.ylim([0.94*max(trace),1.01*max(trace)]) ##use this for positive voltage plots
		#plt.ylim( [130, 165]) ##use this for negative voltage plots
		#plt.xlim([-40,50])
			
			#pylab.savefig(trial+".png") #save fig for each analysis as .png
	plt.legend(trial_library)
	plt.show()
		
		
	
	
if __name__ == '__main__': 
	main() 

