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

def extract_variables(path, file_name):
	file = os.path.join(path, file_name)
	#print file_name		
	csvdata = numpy.loadtxt(file, delimiter=',', usecols = (3,5,8,14,17)) #loading data #3=FWHM, 5=FWQM,8 = baseline,14 = cumulative counts, 15 = mass outflow,16 = height change,17= app pressure
	app_pressure = csvdata[:,4]
	cumulative_counts = csvdata[:,3]
	FWHM = csvdata[:,0]
	FWQM = csvdata[:,1]
	normalised = cumulative_counts/(max(cumulative_counts))
	return app_pressure, cumulative_counts, FWHM, FWQM, normalised

def cubic_fit_data(app_pressure, normalised):
	cubiccoef= numpy.polyfit(app_pressure,normalised,3) # evaluate polynomial
	cubicfit= numpy.polyval(cubiccoef,app_pressure)
	x_inflectioncubic = (-cubiccoef[1])/(3*cubiccoef[0]) #inflection location from cubic
	return cubicfit, x_inflectioncubic

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
	
def lognormal_fit(data):
	shape, loc, scale = stats.lognorm.fit(data)
	distfl = stats.lognorm(shape, loc, scale)
	mu_l = np.log(scale)
	sigma_l = shape
	
	#check for validity of lognormal fitting
	#transformed_data = np.exp(data)
	#p_value = stats.normaltest(transformed_data)
	#print 'normaltest teststat = %6.3f pvalue = %6.4f' % p_value
	#'''
	#"rejects the null hypothesis" when the p-value is less than the significance level,
	 #which is often 0.05 (5%)or 0.01(1%). When the null hypothesis is rejected, the result is said to be statistically significant.
	#'''
	##evaluating normalness of transformed data using 5% cutoff
	#if p_value < 0.05:
		#print "Non-lognormally distributed" #"very small chance that distribution is normal (rejects the null hypothesis)"
		#normalcy = True
	#else:
		#print "Lognormally distributed" #"cannot reject the null hypothesis (i.e. that data is NOT normally distributed)"
		#normalcy = False
	return mu_l,sigma_l, distfl, #normalcy, p_value,
	
def splitting_stats(data,number):
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
	return av_data_list_mean, av_data_list_median, av_data_list_mode, 
	
def splitting_lognormal(data,number):
	#creating averaging variables
	average_number = number
	data_list = data.tolist()#convert arrays to indexed lists
	grouped_data_slices = zip(*(iter(data_list),) * average_number)#split app_pressure_list of tuples of length 'average_number'
	av_data_list_lognormal = []
	#calculates the average from each group of tuples and places it in the list 'av_app_pressure'
	for i in grouped_data_slices:
		values = i
		#print min(values)
		#av_data_list_lognormal.append(lognormal_fit(values))
	return av_data_list_lognormal


def stats_plot(dict_of_trials,trial_name, stats_method):#stats_method 0 = mean, 1 = median, 2=mode, 3=lognormal
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
	
	pstfit = gaussianfit(dict_of_trials[trial_name+'av_app_pressure'][0],dict_of_trials[trial_name+'av_FWHM'][stats_method])
	ngtfit = gaussianfit(dict_of_trials[trial_name+'returnav_app_pressure'][0],dict_of_trials[trial_name+'returnav_FWHM'][stats_method])
	
	plt.plot(dict_of_trials[trial_name+'av_app_pressure'][0],dict_of_trials[trial_name+'av_FWHM'][stats_method],marker='o', color='b')
	plt.plot(dict_of_trials[trial_name+'av_app_pressure'][0],pstfit[0],marker='-', color='b')
	
	plt.plot(dict_of_trials[trial_name+'returnav_app_pressure'][0],dict_of_trials[trial_name+'returnav_FWHM'][stats_method],marker='o', color='r')
	plt.plot(dict_of_trials[trial_name+'returnav_app_pressure'][0],ngtfit[0],marker='-', color='r')
	
	pstroundedmu = round(pstfit[1], 2)
	pstplotinputmu = "$\mu_{"+str(stats_name)+"}$ = "+str(pstroundedmu)
	
	ngtroundedmu = round(ngtfit[1], 2)
	ngtplotinputmu = "$\mu_{"+str(stats_name)+"}$ = "+str(ngtroundedmu)
	
	plt.text(0.75*min(dict_of_trials[trial_name+'av_app_pressure'][0]), 0.6*max(dict_of_trials[trial_name+'av_FWHM'][0]), pstplotinputmu )
	plt.axvline(x=pstroundedmu, ymin=0.01, ymax=0.99, color='b',linestyle=':') #blue vertical line at FWHM max from fit of 0.5V
	
	plt.text(0.75*max(dict_of_trials[trial_name+'returnav_app_pressure'][0]), 0.6*max(dict_of_trials[trial_name+'returnav_FWHM'][0]), ngtplotinputmu )
	plt.axvline(x=ngtroundedmu, ymin=0.01, ymax=0.99, color='r',linestyle=':') #blue vertical line at FWHM max from fit of 0.5V return
	
	plt.title(trial_name)
	plt.ylabel(stats_name+' FWHM (mS)')
	
	return pstroundedmu, ngtroundedmu
	
def main(): 
	trials = {
	'B47':"B47_pH8_0.5V_19.4.2012_techtest_trial2.csv",
	'B47return' : "B47_pH8_-0.5V_19.4.2012_techtest_trial2return.csv",
	'B86':"B86_pH8_0.5V_19.4.2012_techtest.csv",
	'B86return' : "B86_pH8_-0.5V_19.4.2012_techtest_return.csv",
	'B121':"B121_pH8_0.5V_19.4.2012_techtest.csv",
	'B121return': "B121_pH8_-0.5V_19.4.2012_techtestreturn.csv"
	}
	
	trialspst05V = {
	'B47':"B47_pH8_0.5V_19.4.2012_techtest_trial2.csv",
	'B86':"B86_pH8_0.5V_19.4.2012_techtest.csv",
	'B121':"B121_pH8_0.5V_19.4.2012_techtest.csv",
	}
	
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
	
	#generate data averages
	averaging_number = 5
	trial_averages = {}
	for trial in trials:
		#print trial
		trial_averages[trial+'av_app_pressure'] = splitting_stats(trial_dict[trial+'_app_pressure'],averaging_number )
		trial_averages[trial+'av_cumulative_counts'] = splitting_stats(trial_dict[trial+'_cumulative_counts'],averaging_number )
		trial_averages[trial+'av_FWHM']= splitting_stats(trial_dict[trial+'_FWHM'],averaging_number )
		#trial_averages[trial+'av_FWHM_log']= splitting_lognormal(trial_dict[trial+'_FWHM'],averaging_number )
		trial_averages[trial+'av_FWQM']= splitting_stats(trial_dict[trial+'_FWQM'],averaging_number )
		#trial_averages[trial+'av_FWQM_log']= splitting_lognormal(trial_dict[trial+'_FWQM'],averaging_number )
		trial_averages[trial+'av_normalised']= splitting_stats(trial_dict[trial+'_normalised'],averaging_number )
		
	#generate all fits
	trial_fits = {}
	for trial in trials:
		print trial
		trial_fits[trial+'_fit'] = cubic_fit_data(trial_dict[trial+'_app_pressure'], trial_dict[trial+'_normalised'])
		#Raw FWHM plot and fit
		#trial_fits[trial+'_gaussian_fit'] = gaussianfit(trial_dict[trial+'_app_pressure'], trial_dict[trial+'_FWHM'])
		#averaged plot and fit
		trial_fits[trial+'av_gaussian_fit'] = gaussianfit(trial_averages[trial+'av_app_pressure'][0], trial_averages[trial+'av_FWHM'][0])
					
	#plotting cumulative counts
	plt.figure(1)
	
	plt.title("Normalised counts, at pH 8")
	plt.ylabel('Number of Events')
	plt.xlabel('Applied Pressure (mm H2O)')
	
	#trying to annotate plots with inflections etc
	#roundedcubicinflection = round(x_inflectioncubic, 2)
	#plotinputcubicinflection = "$P_{0,cubic}$ = "+str(roundedcubicinflection)
	#plt.text(0.75*min(app_pressure), 0.2*max(cumulative_counts), plotinputpstinflection )
		
	plt.subplot(311)
	a1= plt.plot(trial_dict['B121_app_pressure'],trial_dict['B121_normalised'], 'bo')
	a2= plt.plot(trial_dict['B121return_app_pressure'],trial_dict['B121return_normalised'], 'ro')
	plt.axvline(x=trial_fits['B121_fit'][1], ymin=0.6, ymax=0.99, color='b',linestyle=':')
	plt.axvline(x=trial_fits['B121return_fit'][1], ymin=0.4, ymax=0.8, color='r',linestyle=':')
	plt.ylim([0,1.0])
	plt.legend([a1,a2],["B121 0.5V","B121 -0.5V"])
	
	plt.subplot(312)
	b1= plt.plot(trial_dict['B86_app_pressure'],trial_dict['B86_normalised'], 'bo')
	b2= plt.plot(trial_dict['B86return_app_pressure'],trial_dict['B86return_normalised'], 'ro')
	plt.axvline(x=trial_fits['B86_fit'][1], ymin=0.6, ymax=0.99, color='b',linestyle=':')
	plt.axvline(x=trial_fits['B86return_fit'][1], ymin=0.4, ymax=0.8, color='r',linestyle=':')
	plt.ylim([0,1.0])
	plt.legend([b1,b2],["B86 0.5V","B86 -0.5V"])
	
	plt.subplot(313)
	c1= plt.plot(trial_dict['B47_app_pressure'],trial_dict['B47_normalised'], 'bo')
	c2= plt.plot(trial_dict['B47return_app_pressure'],trial_dict['B47return_normalised'], 'ro')
	plt.axvline(x=trial_fits['B47_fit'][1], ymin=0.6, ymax=0.99, color='b',linestyle=':')
	plt.axvline(x=trial_fits['B47return_fit'][1], ymin=0.4, ymax=0.8, color='r',linestyle=':')
	plt.ylim([0,1.0])
	
	plt.legend([c1,c2],["B47 0.5V", "B47 -0.5V"])
	pylab.savefig(path+"analysis/CC.png")
	##plt.xlim([4,10])
	
	##### plotting FWHM data
	
	plt.figure(2)
	
	plt.title("Mean FWHM, at pH 8")
	plt.ylabel('FWHM Duration (mS)')
	plt.xlabel('Applied Pressure (mm H2O)')
	
	plt.subplot(311)
	d1= plt.plot(trial_averages['B121av_app_pressure'][0],trial_averages['B121av_FWHM'][0], 'bo')
	d2= plt.plot(trial_averages['B121returnav_app_pressure'][0],trial_averages['B121returnav_FWHM'][0], 'ro')
	d3= plt.plot(trial_averages['B121av_app_pressure'][0],trial_fits['B121av_gaussian_fit'][0], 'b-')
	d4= plt.plot(trial_averages['B121returnav_app_pressure'][0],trial_fits['B121returnav_gaussian_fit'][0], 'r-')
	#plt.ylim([0,2.0])
	plt.legend([d1,d2],["B121 0.5V","B121 -0.5V"])
	
	plt.subplot(312)
	e1= plt.plot(trial_averages['B86av_app_pressure'][0],trial_averages['B86av_FWHM'][0], 'bo')
	e2= plt.plot(trial_averages['B86returnav_app_pressure'][0],trial_averages['B86returnav_FWHM'][0], 'ro')
	e3= plt.plot(trial_averages['B86av_app_pressure'][0],trial_fits['B86av_gaussian_fit'][0], 'b-')
	e4= plt.plot(trial_averages['B86returnav_app_pressure'][0],trial_fits['B86returnav_gaussian_fit'][0], 'r-')
	#plt.ylim([0,2.0])
	plt.legend([e1,e2],["B86 0.5V","B86 -0.5V"])
	
	plt.subplot(313)
	f1= plt.plot(trial_averages['B47av_app_pressure'][0],trial_averages['B47av_FWHM'][0], 'bo')
	f2= plt.plot(trial_averages['B47returnav_app_pressure'][0],trial_averages['B47returnav_FWHM'][0], 'ro')
	f3= plt.plot(trial_averages['B47av_app_pressure'][0],trial_fits['B47av_gaussian_fit'][0], 'b-')
	f4= plt.plot(trial_averages['B47returnav_app_pressure'][0],trial_fits['B47returnav_gaussian_fit'][0], 'r-')
	#plt.ylim([0,2.0])
	plt.legend([f1,f2],["B47 0.5V", "B47 -0.5V"])
	pylab.savefig(path+"analysis/FWHM.png")
	#statistics analysis on averaging methods
	
	i= 3
	
	for trial in trialspst05V:
		plt.figure(i)
		
		plt.xlabel('Applied Pressure (mm H2O)')
		##mean
		plt.subplot(311)
		mean_plot = stats_plot(trial_averages,trial, 0)
				
		plt.subplot(312)
		median_plot = stats_plot(trial_averages,trial, 1)
		
		plt.subplot(313)
		mode_plot = stats_plot(trial_averages,trial, 2)
	
		plt.xlabel('Applied Pressure (mmH2O)')
		pylab.savefig(path+"analysis/"+trial+"_stats.png") #save fig for each analysis as .png
		i = i+1	
	#############writing the resuts to a .png file##########
	
	#pylab.show()
	
	plt.close()

	
if __name__ == '__main__': 
	main() 

