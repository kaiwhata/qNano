#############overlaid plots

###format of files for reference###
#StartTime(0),Duration(1),FWHM Start Time(2),Duration FWHM(3),FWQM Start Time(4),FWQM Duration(5),Time of Maximum(6),Magnitude(7),Baseline(8),Bias(9),DeltaX(10),SignalRMS(11),Forced Resets(12),Blockade Orientation(13),cumulative counts(14),mass outflow(15),height change(16),app pressure(17)

import re, string, sys 
import os
import numpy
import pylab
from pylab import *
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.stats import mode
from scipy import stats
import scipy 
import fnmatch


############## plotting and fitting s-curve data #############################

def extract_variables(path, file_name): #extracts variables from csv files and turns them into data lists (format based on command line sioscmd outputs)
	file = os.path.join(path, file_name)
	#print file_name		
	csvdata = numpy.loadtxt(file, delimiter=',', usecols = (3,5,8,14,17)) #loading data #3=FWHM, 5=FWQM,8 = baseline,14 = cumulative counts, 15 = mass outflow,16 = height change,17= app pressure
	app_pressure = csvdata[:,4]
	cumulative_counts = csvdata[:,3]
	FWHM = csvdata[:,0]
	FWQM = csvdata[:,1]
	normalised = cumulative_counts/(max(cumulative_counts))
	return app_pressure, cumulative_counts, FWHM, FWQM, normalised

def cubic_fit_data(dependant_v,independant_v): #does a lst squares cubic fit between dependant_v andindependant_v
	cubiccoef= numpy.polyfit(dependant_v,independant_v,3) # evaluate polynomial
	cubicfit= numpy.polyval(cubiccoef,dependant_v)
	x_inflectioncubic = (-cubiccoef[1])/(3*cubiccoef[0]) #inflection location from cubic
	return cubicfit, x_inflectioncubic

def gaussianfit(dependant_v,independant_v, x_plot_resolution): #gaussian fit between dependant and independant variables
	#default set x_plot_resolution as = dependant_v - resolution dictates how many points the resulting function is evaluated at
		
	# p[0] = amplitude
	# p[1] = mean
	# p[2] = sigma
	fitfunc = lambda p, x: p[0]*scipy.exp(-(x-p[1])**2/(2.0*p[2]**2))+p[3] #target function
	errfunc = lambda p, x, y: fitfunc(p,x)-y #residuals
	
	#find pressure at FWHM maxima 
	a = scipy.where(independant_v==max(independant_v))
	p_equ = dependant_v[a[0]]
	
	sigma = sqrt(abs(sum((independant_v-p_equ)**2*independant_v)/sum(independant_v))) # calculate std dev of distribution (this may not be correct)
	
	p0 = [max(independant_v),p_equ, sigma,min(independant_v)] #statistical guess at initial parameters
	
	p1, success = optimize.leastsq(errfunc, p0, args=(dependant_v, independant_v)) #optimize fitting parameters
	                                     
	corrfit = fitfunc(p1, x_plot_resolution)# compute the best fit function from the best fit parameters
	
	return corrfit, p1[1], p1[2]
		
def splitting_stats(data, average_number): #averages variable 'data' by mean, median and mode groups of 'average_number' datapoints
	#creating averaging variables
	data_list = data.tolist()#convert arrays to indexed lists
	grouped_data_slices = zip(*(iter(data_list),) * average_number)#split app_pressure_list of tuples of length 'average_number'
	av_data_list_mean = [] #define emptylists
	av_data_list_median = []
	av_data_list_mode = []
	#calculates the average from each group of tuples and appends it to the relevant list 
	for i in grouped_data_slices:
		values = i
		av_data_list_mean.append(numpy.average(values)) #mean
		av_data_list_median.append(numpy.median(values)) #media
		av_data_list_mode.append(scipy.stats.mode(values)[0][0]) #mode
	return av_data_list_mean, av_data_list_median, av_data_list_mode, 

def plot_and_gaussian(dict_of_raw_data, trial_name, stats_method, averaging_number, legend):
	#plots input trials averaged by 'stats_method' over 'averaging_number' data points and performs gaussian fits and plots to each
	if stats_method == 0:
		stats_name = 'Mean'
	elif stats_method == 1:
		stats_name = 'Median'
	elif stats_method == 2:
		stats_name = 'Mode'
	elif stats_method == 3:
		stats_name = 'Lognormal'
	else:
		stats_name = 'Unknown'
	
	#locating data in dictionary of raw data then splitting and averaging data	
	av_pressure = splitting_stats(dict_of_raw_data[trial_name+'_app_pressure'],averaging_number)[stats_method]
	av_FWHM = splitting_stats(dict_of_raw_data[trial_name+'_FWHM'] ,averaging_number)[stats_method]
	av_pressure_return = splitting_stats(dict_of_raw_data[trial_name+'return_app_pressure'],averaging_number)[stats_method]
	av_FWHM_return = splitting_stats(dict_of_raw_data[trial_name+'return_FWHM'],averaging_number)[stats_method]
	
	#fits to averaged data	
	pstfit = gaussianfit(av_pressure,av_FWHM,dict_of_raw_data[trial_name+'_app_pressure'])
	ngtfit = gaussianfit(av_pressure_return,av_FWHM_return,dict_of_raw_data[trial_name+'return_app_pressure'])
	
	#setting up plot
	p1 = plt.plot(av_pressure,av_FWHM,'bo')
	p2 = plt.plot(dict_of_raw_data[trial_name+'_app_pressure'],pstfit[0],marker='-', color='b')
	
	p3 = plt.plot(av_pressure_return,av_FWHM_return,'ro')
	p4 = plt.plot(dict_of_raw_data[trial_name+'return_app_pressure'],ngtfit[0],marker='-', color='r')
	
	#annotating plot
	pstroundedmu = round(pstfit[1], 2)
	pstplotinputmu = "$\mu_{"+str(stats_name)+" (n ="+str(averaging_number)+" )}$ = "+str(pstroundedmu)
	
	ngtroundedmu = round(ngtfit[1], 2)
	ngtplotinputmu = "$\mu_{"+str(stats_name)+" (n ="+str(averaging_number)+" )}$ = "+str(ngtroundedmu)
	
	plt.text(0.75*min(av_pressure), 0.6*max(av_FWHM), pstplotinputmu )
	plt.axvline(x=pstroundedmu, ymin=0.01, ymax=0.99, color='b',linestyle=':') #blue vertical line at FWHM max from fit of 0.5V
	
	plt.text(0.3*max(av_pressure_return), 0.6*max(av_FWHM_return), ngtplotinputmu )
	plt.axvline(x=ngtroundedmu, ymin=0.01, ymax=0.99, color='r',linestyle=':') #red vertical line at FWHM max from fit of 0.5V return
	
	plt.title(trial_name)
	plt.ylabel(str(stats_name)+" FWHM (n = "+str(averaging_number)+")")
	plt.xlim([-40,40])
	#legend with plot
	if legend == True:	
		plt.legend([p1,p3],[str(trial_name)+" 0.5V",str(trial_name)+" -0.5V"])
	else:
		print "No Legend"
	return

def plot_annotate_format(value):#rounds and converts to latex string
	rounded = round(value, 2)
	latexstring = "$P_{cubic}$ = "+str(rounded)
	return latexstring 

def cc_plot(dict_of_raw_data, trial_name):
	
	a1= plt.plot(dict_of_raw_data[trial_name+'_app_pressure'],dict_of_raw_data[trial_name+'_normalised'], 'bo')
	a2= plt.plot(dict_of_raw_data[trial_name+'return_app_pressure'],dict_of_raw_data[trial_name+'return_normalised'], 'ro')
	
	#fits
	pstfit = cubic_fit_data(dict_of_raw_data[trial_name+'_app_pressure'],dict_of_raw_data[trial_name+'_normalised'])
	ngtfit = cubic_fit_data(dict_of_raw_data[trial_name+'return_app_pressure'],dict_of_raw_data[trial_name+'return_normalised'])
	
	#trying to annotate plots with inflections etc
	plt.axvline(x=pstfit[1], ymin=0.6, ymax=0.99, color='b',linestyle=':')
	plt.text(0.8*pstfit[1]), 0.5, plot_annotate_format(pstfit[1]))
	
	plt.axvline(x=ngtfit[1], ymin=0.4, ymax=0.8, color='r',linestyle=':')
	plt.text(1.2*pstfit[1]), 0.5, plot_annotate_format(ngtfit[1]))
	
	plt.ylim([0,1.0])
	plt.legend([a1,a2],[trial_name+" 0.5V",trial_name+" -0.5V"])
			
	return

def main(): 
	###full trial set
	trials = {
	'B47':"B47_pH8_0.5V_19.4.2012_techtest_trial2.csv",
	'B47return' : "B47_pH8_-0.5V_19.4.2012_techtest_trial2return.csv",
	'B86':"B86_pH8_0.5V_19.4.2012_techtest.csv",
	'B86return' : "B86_pH8_-0.5V_19.4.2012_techtest_return.csv",
	'B121':"B121_pH8_0.5V_19.4.2012_techtest.csv",
	'B121return': "B121_pH8_-0.5V_19.4.2012_techtestreturn.csv"
	}
	###positive to negative pressure runs trial set
	trialspst05V = {
	'B47':"B47_pH8_0.5V_19.4.2012_techtest_trial2.csv",
	'B86':"B86_pH8_0.5V_19.4.2012_techtest.csv",
	'B121':"B121_pH8_0.5V_19.4.2012_techtest.csv",
	}
	###negative to positive pressure runs trial set
	trialsngt05V = {
	'B47return' : "B47_pH8_-0.5V_19.4.2012_techtest_trial2return.csv",
	'B86return' : "B86_pH8_-0.5V_19.4.2012_techtest_return.csv",
	'B121return': "B121_pH8_-0.5V_19.4.2012_techtestreturn.csv"
	}
	###set file location path##
	path = "/home/jeldridge/current_data/syringe_pump_returns/"
		
	#generate all variables
	trial_dict = {}
	for trial in trials:
		 #print trial
		 trial_dict[trial+'_app_pressure'] = extract_variables(path,trials[trial])[0]
		 trial_dict[trial+'_cumulative_counts'] = extract_variables(path,trials[trial])[1]
		 trial_dict[trial+'_FWHM']= extract_variables(path,trials[trial])[2]
		 trial_dict[trial+'_FWQM']= extract_variables(path,trials[trial])[3]
		 trial_dict[trial+'_normalised']= extract_variables(path,trials[trial])[4]
	
	#generate fits for CC
	trial_fits = {}
	for trial in trials:
		print trial
		trial_fits[trial+'_fit'] = cubic_fit_data(trial_dict[trial+'_app_pressure'], trial_dict[trial+'_normalised'])
			
	#plotting cumulative counts
	plt.figure(1)
	
	plt.title("Normalised counts, at pH 8")
	plt.ylabel('Number of Events')
	plt.xlabel('Applied Pressure (mm H2O)')
		
	plt.subplot(311)
	cc_plot(trial_dict, 'B121')
		
	plt.subplot(312)
	cc_plot(trial_dict, 'B86')
	
	plt.subplot(313)
	cc_plot(trial_dict, 'B47')
	#pylab.savefig(path+"analysis/CC.png")
	##plt.xlim([4,10])
	
	##### plotting FWHM data
	i = 2
	
	plt.figure(i)
	
	plt.xlabel('Applied Pressure (mm H2O)')
	
	
	plt.subplot(311)
	plot_and_gaussian(trial_dict, 'B121', 0, 5,True)
		
	plt.subplot(312)
	plot_and_gaussian(trial_dict, 'B86', 0, 5, True)
	
	plt.subplot(313)
	plot_and_gaussian(trial_dict, 'B47', 0, 5, True)
	
	#pylab.savefig(path+"analysis/FWHM.png")
	
	i=i+1
	for trial in trialspst05V:
		plt.figure(i)
		
		##mean
		plt.subplot(411)
		plot_and_gaussian(trial_dict, trial, 0, 1, False)
				
		plt.subplot(412)
		plot_and_gaussian(trial_dict, trial, 0, 10, False)
		
		plt.subplot(413)
		plot_and_gaussian(trial_dict, trial, 0, 25, False)
		
		plt.subplot(414)
		plot_and_gaussian(trial_dict, trial, 0, 50, False)
		
		plt.xlabel('Applied Pressure (mmH2O)')
		#pylab.savefig(path+"analysis/"+trial+"_means.png") #save fig for each analysis as .png
		
		i = i+1	
		
		plt.figure(i)
		
		###median
		plt.subplot(411)
		plot_and_gaussian(trial_dict, trial, 1, 1, False)
				
		plt.subplot(412)
		plot_and_gaussian(trial_dict, trial, 1, 10, False)
		
		plt.subplot(413)
		plot_and_gaussian(trial_dict, trial, 1, 25, False)
		
		plt.subplot(414)
		plot_and_gaussian(trial_dict, trial, 1, 50, False)
		
		plt.xlabel('Applied Pressure (mmH2O)')
		#pylab.savefig(path+"analysis/"+trial+"_medians.png") #save fig for each analysis as .png
		i = i+1	
		
		plt.figure(i)
	
		##mode
		plt.subplot(411)
		plot_and_gaussian(trial_dict, trial, 2, 1, False)
				
		plt.subplot(412)
		plot_and_gaussian(trial_dict, trial, 2, 10, False)
		
		plt.subplot(413)
		plot_and_gaussian(trial_dict, trial, 2, 25, False)
		
		plt.subplot(414)
		plot_and_gaussian(trial_dict, trial, 2, 50, False)
		
		plt.xlabel('Applied Pressure (mmH2O)')
		#pylab.savefig(path+"analysis/"+trial+"_modes.png") #save fig for each analysis as .png
		i = i+1	
	
	
	#statistics analysis on averaging methods
	for trial in trialspst05V:
		plt.figure(i)
						
		#mean
		plt.subplot(311)
		plot_and_gaussian(trial_dict, trial, 0, 5, False)
		
		#median		
		plt.subplot(312)
		plot_and_gaussian(trial_dict, trial, 1, 5, False)
		
		#mode
		plt.subplot(313)
		plot_and_gaussian(trial_dict, trial, 2, 5, False)
	
		plt.xlabel('Applied Pressure (mmH2O)')
		#pylab.savefig(path+"analysis/"+trial+"_stats.png") #save fig for each analysis as .png
		i = i+1	
	##############writing the resuts to a .png file##########
	
	pylab.show()
	
	plt.close()
	
if __name__ == '__main__': 
	main() 



