#!/usr/bin/python

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
import scipy 
import fnmatch



def create_library(parent, pattern):
	datafile_library = []
	for path, dirs, files in os.walk(os.path.abspath(parent)):
		for filename in fnmatch.filter(files, pattern):
			datafile_library.append(os.path.join(path, filename))
	print datafile_library
	return datafile_library
	
def find_nearest(array,value):
	    idx=(numpy.abs(array-value)).argmin()
	    return array[idx]


### plot and fit Scurve ##
def scurvecubicfit(app_pressure, cumulative_counts):
	cubiccoef= numpy.polyfit(app_pressure,cumulative_counts,3) # evaluate polynomial
	cubicfit= numpy.polyval(cubiccoef,app_pressure)
	inflectioncubic = (-cubiccoef[1])/(3*cubiccoef[0]) #inflection location from cubic
	return cubicfit, inflectioncubic

def scurveparabolicfit(app_pressure, cumulative_counts,x_inflectioncubic, discard):
	x_list = app_pressure.tolist()# convert array tolist
	x_discard = x_inflectioncubic+discard
	x_discardsplit = find_nearest(app_pressure,x_discard)  #find existing x value closest to 'discard'mm around cubic inflection point 
	x_linesplit = x_list.index(x_discardsplit) #search x to find index of matching value (only works if x is list -  not if x is an array)
	
	if discard < 0:
		x_frominflection =  app_pressure[(x_linesplit):]#split data and discard
		y_frominflection =  cumulative_counts[(x_linesplit):]
	elif discard > 0:
		x_frominflection =  app_pressure[:(x_linesplit)]
		y_frominflection =  cumulative_counts[:(x_linesplit)]
	else:
		x_frominflection =  app_pressure[:]
		y_frominflection =  cumulative_counts[:]
		
	pbfit= numpy.polyfit(x_frominflection,y_frominflection,2)#positive parabolic fit
	pbfitplot = numpy.polyval(pbfit,x_frominflection) #plot to checkpyplot
		
	inflectionpbfit = (-pbfit[1])/(2*pbfit[0]) #inflection location from positive parabola

	return pbfitplot, inflectionpbfit, x_frominflection, y_frominflection
	
def scurveplot(app_pressure, cumulative_counts,cubicfit_data,x_inflectioncubic,x_aboveinflection, y_aboveinflection,pstpbfit_data, x_belowinflection, y_belowinflection,ngtpbfit_data,x_inflectionpstvepbfit, x_inflectionngtvepbfit):
	plt.plot(app_pressure,cumulative_counts,linestyle='-.', color='k',)
	plt.plot(app_pressure,cubicfit_data) #display fit on plot
	
	plt.plot(x_aboveinflection, y_aboveinflection, marker='.', color='b') #plot split data as a check
	plt.plot(x_belowinflection, y_belowinflection, marker='.', color='b') #plot split data as a check
	plt.plot(x_aboveinflection, pstpbfit_data, color='r')
	plt.plot(x_belowinflection, ngtpbfit_data, color='r')
	
	plt.ylabel('Cumulative Counts')
	#plt.legend(('Raw Data','Cubic Fit','Filtered Data', 'Filtered Data', 'Parabolic Fits'), loc=2)
	plt.axvline(x=x_inflectioncubic, ymin=0.01, ymax=0.99, color='g',linestyle=':') #green vertical line at cubic inflection
	plt.axvline(x=x_inflectionpstvepbfit, ymin=0.25, ymax=0.75, color='r',linestyle=':') #red vertical line at pstve parabolic inflection
	plt.axvline(x=x_inflectionngtvepbfit, ymin=0.25, ymax=0.75, color='r',linestyle=':') #red vertical line at ngtve parabolic inflection
	plt.xlim([-37,30])
	
	roundedcubicinflection = round(x_inflectioncubic, 2)
	roundedpositiveinflection = round(x_inflectionpstvepbfit, 2)
	roundednegativeinflection = round(x_inflectionngtvepbfit,2)
	
	plotinputcubicinflection = "$P_{0,cubic}$ = "+str(roundedcubicinflection)
	plotinputpstinflection = "$P_{0,+ve}$ = "+str(roundedpositiveinflection)
	plotinputngtinflection = "$P_{0,-ve}$ = "+str(roundednegativeinflection)
	
	plt.text(0.75*min(app_pressure), 0.2*max(cumulative_counts), plotinputpstinflection )
	plt.text(0.75*min(app_pressure), 0.1*max(cumulative_counts), plotinputngtinflection )
	plt.text(0.75*min(app_pressure), 0.3*max(cumulative_counts), plotinputcubicinflection )
	
### plot and fit FWHM ##
def splitting(data,number):
	#creating averaging variables
	average_number = number
	data_list = data.tolist()#convert arrays to indexed lists
	grouped_data_slices = zip(*(iter(data_list),) * average_number)#split app_pressure_list of tuples of length 'average_number'
	av_data_list_mean = [] #define emptylists
	av_data_list_median = []
	av_data_list_mode = []
	#calculates the average from each group of tuples and places it in the list 'av_app_pressure'
	for i in grouped_data_slices:
		values = i
		av_data_list_mean.append(numpy.average(values))
		av_data_list_median.append(numpy.median(values))
		av_data_list_mode.append(scipy.stats.mode(values)[0][0])
	return av_data_list_mean, av_data_list_median, av_data_list_mode

def gaussianfit(dependant_v,independant_v):
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
	                                     
	corrfit = fitfunc(p1, dependant_v)# compute the best fit function from the best fit parameters
	#statfit = fitfunc(p0, dependant_v)# estimate best fit paramters from statistical analysis
	return corrfit, p1[1], p1[2]
	

def fwhmplot(pressure,FWHM, corrfit, mu, sigma, FWHM_average_number):
	##############################################FWHM plotting time! ####################3	
	#plotting original FWHM data
	plt.plot(pressure,FWHM,marker='o', color='k')
	plt.plot(pressure,corrfit, color='r')
	
	plt.ylabel('Average FWHM Duration (mS)')
	#plt.legend(('Averaged Data, n = '+str(average_number),'Fitted Gaussian'), loc=2)
	plt.xlim([-37,30])
	
	roundedmu = round(mu, 2)
	roundedsigma = round(sigma, 2)
	#roundedFWHM = round(p_equ,2)
	plotinputmu = "$\mu$ (n = "+str(FWHM_average_number)+") = "+str(roundedmu)
	plotinputsigma = "$\sigma$ (n = "+str(FWHM_average_number)+") = %s" % (roundedsigma)  
	#plotinputFWHMmax  = "P @ $FWHM_{max}$ = "+str(roundedFWHM)
	
	plt.text(0.75*min(pressure), 0.6*max(FWHM), plotinputmu )
	plt.text(0.75*min(pressure), 0.5*max(FWHM), plotinputsigma )
	#plt.text(0.7*min(pressure), 0.7*max(FWHM), plotinputFWHMmax )
	plt.axvline(x=roundedmu, ymin=0.01, ymax=0.99, color='b',linestyle=':') #blue vertical line at FWHM mean from fit
	
### plot and fit current trace ##
def baselineplot(app_pressure, trace,av_norm_pressure, av_blockade_magnitude,average_blockade_magnitude, average_trace):
	
	ax1 = plt.subplot(313)
	ax1.plot(av_norm_pressure, av_blockade_magnitude, 'r.')
	ax1.set_xlabel('Average Applied Pressure (mmH20)')
	ax1.set_ylabel('Blockade Height')
	ax1.set_ylim(0,0.4)
	
	geometry_indicator = (average_blockade_magnitude/average_trace)
	geometry_text = '${\Delta I}/{I}$ = '+'%.3e' % (geometry_indicator)
	
	plt.text(-17, 0.35, geometry_text )
	
	ax2 = twinx()
	ax2.plot(app_pressure, trace, 'k-')
	ax2.set_ylabel('Baseline current (nA)')
		
	if min(trace) < 0:
		plt.ylim([1.01*min(trace), 0.94*min(trace)]) ##use this for negative voltage plots
	else:
		plt.ylim([0.94*max(trace),1.01*max(trace)]) ##use this for positive voltage plots
	plt.xlim([-37,30])
	return geometry_indicator
	
def error_analysis(mean_roundedmu,median_roundedmu,mode_roundedmu,x_inflectioncubic, x_inflectionpstvepbfit, x_inflectionngtvepbfit):
	##########################error analysis (or estimate of variance between different fitting methods##############
	#set of all balance points
	balance_points = [mean_roundedmu,median_roundedmu,mode_roundedmu, x_inflectioncubic, x_inflectionpstvepbfit, x_inflectionngtvepbfit]
	balance_point = numpy.average(balance_points)
	ranges = (max(balance_points)-min(balance_points))/2
	if ranges > 4:
		print ranges
		print "you've got a disagreement here boss!"
	else:
		print ranges
		print "agreement between different methods"
	return balance_point, ranges
	
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
	
def statisticsplot(mean_FWHM, median_FWHM, mode_FWHM,av_app_pressure):
	plt.subplot(311)
	plt.plot(av_app_pressure,mean_FWHM,marker='o', color='k')
	mean_fit = gaussianfit(av_app_pressure,mean_FWHM)[0]
	mean_roundedmu = round(gaussianfit(av_app_pressure,mean_FWHM)[1], 2)
	mean_plotinputmu = "$\mu_{mean}$ = "+str(mean_roundedmu)
	plt.text(0.75*min(av_app_pressure), 0.6*max(median_FWHM), mean_plotinputmu )
	plt.axvline(x=mean_roundedmu, ymin=0.01, ymax=0.99, color='b',linestyle=':') #blue vertical line at FWHM max from fit
	plt.plot(av_app_pressure,mean_fit,marker='-', color='r')
	plt.ylabel('Mean FWHM (mS)')
	
	plt.subplot(312)
	plt.plot(av_app_pressure,median_FWHM,marker='o', color='k')
	median_fit = gaussianfit(av_app_pressure,median_FWHM)[0]
	median_roundedmu = round(gaussianfit(av_app_pressure,median_FWHM)[1], 2)
	median_plotinputmu = "$\mu_{median}$ = "+str(median_roundedmu)
	plt.text(0.75*min(av_app_pressure), 0.6*max(median_FWHM), median_plotinputmu )
	plt.axvline(x=median_roundedmu, ymin=0.01, ymax=0.99, color='b',linestyle=':') #blue vertical line at FWHM max from fit
	plt.plot(av_app_pressure,median_fit,marker='-', color='r')
	plt.ylabel('Median FWHM (mS)')
	
	plt.subplot(313)
	plt.plot(av_app_pressure,mode_FWHM,marker='o', color='k')
	mode_fit = gaussianfit(av_app_pressure,mode_FWHM)[0]
	mode_roundedmu = round(gaussianfit(av_app_pressure,mode_FWHM)[1], 2)
	mode_plotinputmu = "$\mu_{mode}$ = "+str(mode_roundedmu)
	plt.text(0.75*min(av_app_pressure), 0.6*max(mode_FWHM), mode_plotinputmu )
	plt.axvline(x=mode_roundedmu, ymin=0.01, ymax=0.99, color='b',linestyle=':') #blue vertical line at FWHM max from fit
	plt.plot(av_app_pressure,mode_fit,marker='-', color='r')
	plt.ylabel('Mode FWHM (mS)')
	plt.xlabel('Applied Pressure (mmH2O)')
	return mean_roundedmu, median_roundedmu, mode_roundedmu
	
	

def main(): 
	parent = "/home/jeldridge/current_data/syringe_pump_returns/"
	pattern = "*.csv"
	datafile_library = create_library(parent, pattern)
	
	#creating averaging variables
	FWHM_average_number = 5
	blockade_average_number = 50
	
	output=open(parent+"plots/"+"complete_results_summary.csv","a")  #save results to output file
	output.write('%s , %s , %s , %s , %s , %s , %s , %s , %s, %s , %s, %s , %s ,%s ,%s\n'%('Trial Name','Inflection point from cubic fit=','Inflection point from positive parabola fit=','Inflection point from negative parabola fit=', "FWHM fit Mean","FWHM fit Std Dev","FWHM fit Median","FWHM fit Mode", "Averaged over","Av Balance Point","Range/Error in points", "pH", "particle charge","av blockade height","deltaI over I")) # print the header
	
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
		av_norm_pressure = splitting(app_pressure,blockade_average_number)[0]
		av_app_pressure = splitting(app_pressure,FWHM_average_number)[0]
		av_blockade_magnitude = splitting(blockade_magnitude,blockade_average_number)[0]
		mean_FWHM = splitting(FWHM,FWHM_average_number)[0]
		median_FWHM = splitting(FWHM,FWHM_average_number)[1]
		mode_FWHM = splitting(FWHM,FWHM_average_number)[2]
		av_normalised_blockades = splitting(normalised_blockades,blockade_average_number)[0]
		
		plt.title(trial)
		plt.figure(1)
		
		plt.subplot(311)
		
		scurvecubic_result = scurvecubicfit(app_pressure, cumulative_counts)
		cubicfit_data = scurvecubicfit(app_pressure, cumulative_counts)[0]
		x_inflectioncubic = scurvecubicfit(app_pressure, cumulative_counts)[1]
		
		pstpbfit_data = scurveparabolicfit(app_pressure, cumulative_counts,x_inflectioncubic,5)[0]
		x_inflectionpstvepbfit = scurveparabolicfit(app_pressure, cumulative_counts,x_inflectioncubic,5)[1]
		x_aboveinflection = scurveparabolicfit(app_pressure, cumulative_counts,x_inflectioncubic,5)[2]
		y_aboveinflection = scurveparabolicfit(app_pressure, cumulative_counts,x_inflectioncubic,5)[3]
		
		ngtpbfit_data = scurveparabolicfit(app_pressure, cumulative_counts,x_inflectioncubic,-5)[0]
		x_inflectionngtvepbfit = scurveparabolicfit(app_pressure, cumulative_counts,x_inflectioncubic,-5)[1]
		x_belowinflection = scurveparabolicfit(app_pressure, cumulative_counts,x_inflectioncubic,-5)[2]
		y_belowinflection = scurveparabolicfit(app_pressure, cumulative_counts,x_inflectioncubic,-5)[3]
				
		scurveplot(app_pressure, cumulative_counts,cubicfit_data,x_inflectioncubic,x_aboveinflection, y_aboveinflection,pstpbfit_data, x_belowinflection, y_belowinflection,ngtpbfit_data, x_inflectionpstvepbfit, x_inflectionngtvepbfit)
		
		plt.subplot(312)
		gaussian_fit_result = gaussianfit(av_app_pressure,mean_FWHM)
		fwhmplot(av_app_pressure,mean_FWHM,gaussian_fit_result[0],gaussian_fit_result[1], gaussian_fit_result[2], FWHM_average_number)
		
		plt.subplot(313)
		geometry_indicator = baselineplot(app_pressure, trace,av_norm_pressure, av_blockade_magnitude,average_blockade_magnitude, average_trace)
		
		#plt.show()
		pylab.savefig(trial+".png") #save fig for each analysis as .png
		plt.close()
		
		plt.figure(2)
		mean_roundedmu = statisticsplot(mean_FWHM, median_FWHM, mode_FWHM,av_app_pressure)[0]
		median_roundedmu = statisticsplot(mean_FWHM, median_FWHM, mode_FWHM,av_app_pressure)[1]
		mode_roundedmu = statisticsplot(mean_FWHM, median_FWHM, mode_FWHM,av_app_pressure)[2]
		#plt.show()
		pylab.savefig(trial+"_stats.png") 
		plt.close()
				
		balance_point = error_analysis(mean_roundedmu,median_roundedmu,mode_roundedmu,x_inflectioncubic, x_inflectionpstvepbfit, x_inflectionngtvepbfit)[0]
		ranges = error_analysis(mean_roundedmu,median_roundedmu,mode_roundedmu,x_inflectioncubic, x_inflectionpstvepbfit, x_inflectionngtvepbfit)[1]
		
		trial_name = correct_name(trial)[0]
		particle_charge = correct_name(trial)[1]
		pH = correct_name(trial)[2]
		
		output.write('%s , %f , %f , %f , %f , %f , %f , %f ,%f , %f , %f , %s , %s, %s, %s\n'%(trial_name,x_inflectioncubic,x_inflectionpstvepbfit,x_inflectionngtvepbfit,  mean_roundedmu, gaussian_fit_result[2], median_roundedmu, mode_roundedmu, FWHM_average_number,balance_point, ranges, pH, particle_charge,average_blockade_magnitude, geometry_indicator)) # print the data
			
	output.close	

if __name__ == '__main__': 
	main() 

