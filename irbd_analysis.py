#!/usr/bin/python

import re, os
import numpy as np
import pylab
from pylab import *
import matplotlib.pyplot as plt
from scipy import optimize
import scipy 
import csv
import sys
import math
import itertools
import time

time_i = time()
path = '/home/jeldridge/irbd_test.dat'
trace_averaging_window_mS = 200.0
datapoints_in_trace_av_window = (trace_averaging_window_mS/1000)*20000
print datapoints_in_trace_av_window
threshold_percentage = 0.05 #%

"""
Loading binary
"""
full_data_file = open(path, "r")

loadarray = np.empty([1,2])

#define data in array with time = array[i,0] and trace = array[i,1] 
for line in full_data_file:
	linesplt = line.rstrip()
	cols = linesplt.split()
	loadarray = np.append(loadarray, [[float(cols[0]), float(cols[1])]], axis=0)

#deletes first 'empty' row
loadarray = np.delete(loadarray, 0, 0)

"""
Describing discard region
"""
print "No. of discarded datapints at ends"
d = int(0.5*datapoints_in_trace_av_window)
print d
print "Total no. of datapoints"
t = int(len(loadarray[:,0]))
print t

#discards end data points from loadarray
ies = np.arange(d, (t-d))
print ies

"""
Creates an array with the characteristics: time, trace, moving average, difference 
"""
#creates an array with the dimensions required for polulation
array = np.empty([1,8])
print array.shape

#time = array[i,0]
#trace = array[i,1]
#moving average = array[i,2]
#RMS = array[i,3]

for i in ies: #calculated moving average of points on either side of i 
	average = sum(loadarray[i-d:i+d,1])/datapoints_in_trace_av_window
	RMS_noise = abs(math.sqrt(sum(loadarray[i-d:i+d,1]**2)/datapoints_in_trace_av_window) - average)
	array = np.append(array, [[loadarray[i,0], loadarray[i,1], average,RMS_noise, 0, 0, 0, 0]], axis=0)

#discards 'empty' first datapoint
array = np.delete(array, 0, 0)

#append difference between background and trace = array[i,4]
#append trace data for threshold breaching points = array[i,5]
for line in array:
	line[4] = (abs(line[1]-line[2]))
	if 	line[4] >= (threshold_percentage/100)*(line[2]):
		line[5] = line[1]
	else:
		line[5] = 0


"""
Extracting Event Data
"""
def find_first(list, number, threshold):
	for i in list:
		if abs(i-number) <= threshold:
			return list.index(i)

#scans through array[i,5] and creates array[i,6] which contains more info on each event
j = 0

###
#This step could be improved by changing the 'return' requirement that stipulates the end of events to an average returning to the baseline rather than individual events
###

for line in array:
	if j+1 == len(array[:,5]):
		break
	elif abs(line[5] - array[(j-1),5]) >= 100.0: #if at the beginning of an event
		print 'extending event on line'+str(j)
		index = find_first(list(reversed(array[0:j,4])),0.0, 0.001) #checks for closest proceeding point within 0.001 of the baseline
		print index
		for i in range(index):
			array[(j-i),6] = array[(j-i),1] #copies datapoint to array[i,6]
	elif abs(line[5] - array[(j+1),5]) >= 100.0: #if at the end of an event
		print 'extending event on line'+str(j)
		index = find_first(list(array[j:,4]),0.0, 0.001)  #checks for closest following point within 0.001 of the baseline
		print index
		for i in range(index):
			array[(j+i),6] = array[(j+i),1] #copies datapoint to array[i,6]
	j=j+1

for line in array:
	if not line[6] == 0:
		line[7] = line[0] 

def isplit(iterable,splitters):
    return [list(g) for k,g in itertools.groupby(iterable,lambda x:x in splitters) if not k]

def extract_events(list):
	#separates between zero entires and individual events
	events = isplit(list,(0,))
	return events

event_points = extract_events(array[:,6])
time_points = extract_events(array[:,7])
print "Detected "+str(len(event_points))+" events"

#examine individual events
#plt.figure(1)
#plt.plot(time_points[1], event_points[1], 'ko-')
#plt.show()

"""
Calculating event characteristics
"""
Output_array = np.empty([1,8])
#Maximum
#FWHM
#Start time
#Baseline duration
#Asymmetry
#Translocation direction
#FWQM(?)


"""
Writing results to CSV
"""

##output.csv
#ofile = open(path+'.csv', 'wb')
#writer = csv.writer(ofile)
#writer.writerow( ('Time (s)', 'Current (nA)', 'Averaged Baseline Current(nA)' ) )

#ofile.close()	

"""
Plotting
"""
plt.figure(2)

#plots trace 
plt.plot(loadarray[:,0], loadarray[:,1], 'k-')
#plots threshold points
plt.plot(array[:,0], array[:,5], 'ro')

#plots preceeding threshold points
plt.plot(array[:,0], array[:,6], 'bo')

#plots threshold cut-off region
plt.plot(array[:,0], (1.0+threshold_percentage/100)*array[:,2], 'g:')
plt.plot(array[:,0], (1.0-threshold_percentage/100)*array[:,2], 'g:')

#plots baseline current
plt.plot(array[:,0],array[:,2], 'r:')

plt.axhline(y=1.001*array[0,1], xmin=0.1*array[-1,0], xmax=(0.1*array[-1,0]+trace_averaging_window_mS/1000), color='g',linestyle='-')

plt.ylim(130, 135)
plt.ylabel('Current (nA)')
plt.xlabel('Time (s)')
plt.show()

time_f = time()

print 'Calculation time = '+str(time_f-time_i)+' s'



'''
additional crap
'''
#b = 0
#for line in full_data_file:
	#list_of_traces = []
	#j = datapoints_in_trace_av_window
	#a = b
	#i = 0
	#for line in full_data_file:
		#if i < b:
			#i = i+1
		#elif a < (j+b):
			#linesplt = line.rstrip()
			#cols = linesplt.split()
			#list_of_traces.append(float(cols[1]))
			#a = a+1
			#i = i+1
		#else:
			#break
		
	#linesplt = line.rstrip()
	#cols = linesplt.split()
	#cols.append(average(list_of_traces))
	#b=b+1
	#writer.writerow(cols)

