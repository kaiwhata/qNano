
#Author: Elf Eldridge
#email: kaiwhata@gmail.com
#Date: 16/12/2013
#License: None
#version: 0.0.1 (alpha)

#Notes - this is custom software for use in analysing odd event shapes in .irbd files. 
#It will likely fail if asked to analyse typical event traces

from Tkinter import *
import tkFileDialog
import tkMessageBox
import numpy

#setting window parameters
window = Tk() 
window.title('Irbd analysis tool')
window.minsize(width=1000, height=500)

#defining load options
def doOpen():
	file = tkFileDialog.askopenfile(mode='r')
	#fileContents = doRead#this will need to be modified
	
	#read the filename and location
	fileTitle = file.name	
	
	#CHECK FILENAME to see if has been given an .irbd
	if not fileTitle[-5:] == '.irbd':
		tkMessageBox.showwarning('Warning','The file you have selected aint an .irbd ')	
	else:
		inpt_fln.set(fileTitle)
		inpt_dat.set(fileTitle[:-5]+".dat")

def doSample():
	#setting preview window parameter
	window_preview = int(preview_time.get()) #in no of data pointa at 50kHz
	import subprocess
	
	sample=[]
	def run_command(command):
   		p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True) #stderr=subprocess.STDOUT
   		return iter(p.stdout.readline, b'')

	command = "od -A d -t f8 -w32 -N"+str(window_preview*32)+" "+str(inpt_fln.get()) #32 bytes per output line

	#load first Xs of data
	print "Running the subprocess %s" % command

	for line in run_command(command):
	    line = line.rstrip()
	    cols = line.split()
	    if len(cols) ==5: #ckecks line is complte prior to appending it
		sample.append([float(k) for k in cols[1:3]])
	    else:
		print line
   	
	sample_arr = numpy.array(sample)
	sample_plot(sample_arr)


def doConvert():
	import os
	from time import time
	time_i=time()
	statinfo = os.stat(inpt_fln.get()) 
	print "filesize is (in MB):" 
	sz = (float(statinfo.st_size)/1032.0**2)
	print sz
	alert = "this process is likely to take about %s seconds" % round(0.25*sz, 1)
	tkMessageBox.showwarning('Warning',alert)
	os.system("od -A d -t f8 -w32 "+str(inpt_fln.get())+" > "+str(inpt_fln.get())[:-5]+".dat")
	time_f=time()
	tkMessageBox.showinfo('Info','System call complete. The .dat file should now be ready for acessing')	
	print "System call complete. The .dat file should now be ready for acessing"
	print "time taken = %s seconds" % round((time_f-time_i), 2)

	statinfo_dat = os.stat(str(inpt_fln.get())[:-5]+".dat") 
	print "output filesize is (in MB):" 
	sz_dat = (float(statinfo_dat.st_size)/1032.0**2)
	print sz_dat
	if sz_dat >= 1500.0: #checks if input file size is larger than 1G 
		alert = "The produced file is more than 1G in size (%s MB). This will likely cause memory issues when you run the analysis. I suggest using the split command to split files first." % round(sz, 1)
		tkMessageBox.showwarning('Warning',alert)

	inpt_dat.set(str(inpt_fln.get())[:-5]+".dat")



def read_dat():
	import numpy

	try:	
		#name = str(inpt_fln.get())[:-5]+".dat"
		name = inpt_dat.get()
		ifile = open(name, 'rb')
		print "reading file to memory"
		data = []
		i=0
		for row in ifile:
		    row = row.rstrip()
		    row=row.split()
		    i=i+1
		    if len(row) == 5:
			data.append([float(row[1]),float(row[2])]) 
		    else:
			print row

		ifile.close
		print "input file closed"
		arr = numpy.array(data)
		print arr.shape
		return arr
	except:
		Message = "I cant find the .dat file. Please make sure you have Coverted from the binary one (with the convert button) and that the .dat file is still in the folder this script is running from. It will have the same name as your input file"
		tkMessageBox.showwarning('Warning',Message)	
def doRead():
	
	import os
	from time import time
	time_i=time()
	statinfo = os.stat(str(inpt_fln.get())[:-5]+".dat") 
	print "filesize is (in MB):" 
	sz = (float(statinfo.st_size)/1032.0**2)
	print sz
	alert = "This is a highly memory intensive process. This process is likely to take about %s seconds. Check the Command Line to see if it has stalled or not (it shouldnt)" % round(0.37*sz, 1)
	tkMessageBox.showwarning('Warning',alert)
	

	#setting preview window parameter
	window_preview = int(preview_time.get()) #in no of data pointa at 50kHz

	time_i = time()

	#progress = ttk.Progressbar(window,length=7)
	#progress.pack(side=RIGHT)

	print "input file closed"
	arr = read_dat()
	#update progressbar	
	#progess.step(1)
	#progress.update()	
	
	print "now im calculating a moving average based on the window value %s" % av_window.get()
	window_array = mv_av(arr, float(av_window.get()))
	print window_array.shape
	arr = []

	#update progressbar	
	#progess.step(1)
	#progress.update()	
	
	print "now i'll generate an example output for the first %s seconds of the file" % (window_preview/50000.0)
	#plot(window_array)

	print "now I'm beginning to check for points above and below the threshold values"
	event_array = event_detection(window_array)
	
	#update progressbar	
	#progess.step(1)
	#progress.update()	
	
	print "now im cross-referencing that against your estimate of event duration"
	final_event_array = event_check(event_array, window_array, float(event_window.get()))

	#update progressbar	
	#progess.step(1)
	#progress.update()	
	
	event_plot(window_array, final_event_array)
	print "Now im writing the detected files out to the .csv file you specified as %s" % filename.get()

	#update progressbar	
	#progess.step(1)
	#progress.update()	
	
	import csv
	ofile = open(str(filename.get()), 'w')
	writer = csv.writer(ofile, delimiter=',')
	b=0
	writer.writerow(["Mid Time","Start Time","Max","Max_t", "Max_diff", "Min","Min_t", "Min_diff", "Av baseline", "max_min_time difference", "FWHM_i_start","FWHM_i_end", "FWHM_i", "FWHM_f_start","FWHM_f_end", "FWHM_f" ])
#mid_time, start_time, Max, Max_t, Max_diff, Min, Min_t, Min_diff, av_baseline, max_min_duration, FWHM_i_start,FWHM_i_end, FWHM_i, FWHM_f_start,FWHM_f_end, FWHM_f 
	for line in final_event_array:
		writer.writerow(line)
		b=b+1	
	ofile.close()	
	time_f = time()
	#update progressbar	
	#progess.step(1)
	#progress.update()	
	
	#progress.destroy()

	completion_string = "Processing Complete! "+str(b)+" detected events have been written to "+filename.get()
	tkMessageBox.showwarning('Warning', completion_string)	
	print "time taken"	
	print (time_f-time_i)

def doSplit():
	import os, math
	file_to_split = str(inpt_fln.get())[:-5]+".dat" 
	statinfo = os.stat(file_to_split) 
	print "filesize is (in MB):" 
	sz = (float(statinfo.st_size)/1032.0**2)
	target_size = float(split_size.get())
	print sz
	print target_size
	alert = "The input file is "+str(round(sz, 1))+" MBs. Im going to split it into "+str(math.ceil(sz/target_size))+" files"
	tkMessageBox.showwarning('Warning',alert)

	#500MB Should be (at about 9830 lines per MB) about 4915500 lines
	print math.ceil(sz/target_size)
	print range(int(math.ceil(sz/target_size)))
	section_ends = [j*(9830*int(target_size)) for j in range(int(math.ceil(sz/target_size)))]
	print section_ends
	import numpy, csv
	batch_filename = ''
	q=0	
	for z in section_ends:
		try:				
			ifile = open(file_to_split, 'rb')
			ofile = open(file_to_split[:-4]+"_"+str(q)+".dat", 'w')
			writer = csv.writer(ofile, delimiter=' ')

			i=0
			for row in ifile:
				if i >=z and i < (z+(9830*int(target_size))):  
					row = row.rstrip()
		    			row=row.split()				    
					writer.writerow(row) 
				i+=1
			ifile.close
			ofile.close
			batch_filename+file_to_split[:-4]+"_"+str(q)+".dat,"
		except:
			Message = "I cant find the .dat file. Please make sure you have Coverted from the binary one (with the convert button) and that the .dat file is still in the folder this script is running from. It will have the same name as your input file"
			tkMessageBox.showwarning('Warning',Message)
		q+=1	
	
	inpt_dat.set(batch_filename)	
	

#there are now the analysis methods
def runningMeanFast(x, N):
	import numpy as np
	weightings = np.repeat(1.0, N) / N
	truncated = np.convolve(x, weightings)[N-1:-(N-1)]
	#creating false ends
	return np.convolve(x, weightings)

def mv_av(input_array, averaging_window):
	import numpy as np
	b = runningMeanFast(input_array[:,1], averaging_window)
	data_2 = []
	for i in range(len(input_array[:,0])):
		row = input_array[i,:].tolist()
		row.append(b[i]) #with moving average with difference from av
        	row.append(b[i]-input_array[i,1])
		data_2.append(row)
	out_array = np.array(data_2)
	return out_array

def event_detection(window_array):
	import numpy
	avg_window = float(av_window.get())
	threshold_value = float(threshold.get())

	event_times = []
	#we set up the system to go into an array with [:,0] = dips and [:,1] = peaks
	j=0
	for row in window_array:
		if not j<avg_window or j>(len(window_array[:,0])-avg_window):
            #ignores the ending events from averaging window
            #we want to detect something that dips below the average >0.1% and then immediately goes bove the average by >0.1% within the period of 0.001 seconds
			if row[3]/row[2] >= threshold_value/100: #looks for dips
				event_times.append([row[0],row[3]])
			if row[3]/row[2] <= -threshold_value/100: #looks for peaks
				event_times.append([row[0],row[3]])
        	j=j+1
   
	print len(event_times)
	event_arr = numpy.array(event_times)
	print event_arr.shape
	return event_arr

def find_nearest(array,value):
	import numpy
	if not len(array)==0:
    		idx = (numpy.abs(array-value)).argmin()
    		return idx
	else:
		return 0

def event_check(event_array, window_array, event_window):
	import os
	if event_csvs.get() == True:
		#creates directory for output files
		os.system("mkdir "+str(inpt_fln.get())[:-5])

	import numpy    
	output_list = []
	for u in range(len(event_array)-1):
		if (event_array[(u+1),0]-event_array[u,0]) >= (float(event_window)/1000.0):
	            output_list.append(event_array[u,:].tolist())
    
    	print 'refined event number'        
    	print len(output_list)
    	#here we search the window_array time column for matching time values
    	#then we take a sectoin of the data the event width on either side of it
    	#and look for the maxima and minima in these regions
    	output_form = []
	o=0
    	for event in output_list:
        	itemindex=numpy.where(window_array[:,0]==event[0])
        	indint = int(itemindex[0])
        	#convert event window into number of samples @ 50 kHz
        	samples = (event_window/1000.0)*(50000.0)   
      		
        	section = window_array[(indint-samples):(indint+samples),1]
		section_times = window_array[(indint-samples):(indint+samples),0]
	
        	Max = max(section)
        	Min = min(section)

		Max_t = window_array[(int(numpy.where(section==Max)[0])+(indint-samples)),0]
		Min_t = window_array[(int(numpy.where(section==Min)[0])+(indint-samples)),0]
        
        	av_baseline = numpy.average(window_array[(indint-samples):(indint+samples),2])
        
        	start_time = window_array[(indint-samples),0]
		mid_time = window_array[indint,0]

		max_min_duration =  abs(Max_t-Min_t)
		
		diff_from_baseline = numpy.subtract(window_array[(indint-samples):(indint+samples),2],window_array[(indint-samples):(indint+samples),1])
		Max_diff = max(diff_from_baseline)
		Min_diff = min(diff_from_baseline)
			
		if attempt_FWHM.get() == True:

			MaxIndx = int(numpy.where(diff_from_baseline==Max_diff)[0])
			MinIndx = int(numpy.where(diff_from_baseline==Min_diff)[0])

			#search split arrays for nearest values
			half_max_ind_1 = find_nearest(diff_from_baseline[:MaxIndx], Max_diff/2)
			half_max_ind_2 = find_nearest(diff_from_baseline[MaxIndx:], Max_diff/2)+MaxIndx
		
			half_min_ind_1 = find_nearest(diff_from_baseline[:MinIndx], Min_diff/2)
			half_min_ind_2 = find_nearest(diff_from_baseline[MinIndx:], Min_diff/2)+MinIndx
		
			#
			FWHM_i_start = section_times[half_max_ind_1]
			FWHM_i_end = section_times[half_max_ind_2]
			FWHM_f_start = section_times[half_min_ind_1]
			FWHM_f_end = section_times[half_min_ind_2]
			#
			FWHM_i = FWHM_i_end-FWHM_i_start
			FWHM_f = FWHM_f_end-FWHM_i_start
			
			output_form.append([mid_time, start_time, Max, Max_t, Max_diff, Min, Min_t, Min_diff, av_baseline, max_min_duration, FWHM_i_start,FWHM_i_end, FWHM_i, FWHM_f_start,FWHM_f_end, FWHM_f ])
        	else:
			output_form.append([mid_time, start_time, Max, Max_t, Max_diff, Min, Min_t, Min_diff, av_baseline, max_min_duration, 0, 0, 0, 0, 0, 0])
		o=o+1
		#checking for output to csv checkbox
		if event_csvs.get() == True:
			import csv
			#full event data
			full_event_data = window_array[(indint-samples):(indint+samples),:2]
			#write to file			
			oi_file = open(str(inpt_fln.get())[:-5]+'/ event'+str(o)+'.csv', 'w')
			writer = csv.writer(oi_file, delimiter=',')
			writer.writerow(["Time (s)","Current (nA)"])
			for line in full_event_data:
				writer.writerow(line)
			oi_file.close()			
    	final_output = numpy.array(output_form)
    	return final_output

def event_plot(window_array, event_array):
	#setting preview window parameter
	window_preview = int(preview_time.get()) #in no of data pointa at 50kHz
	from numpy import average
	import matplotlib.pyplot as plt
	from matplotlib.figure import Figure
	from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

	fig = Figure(figsize=(18,4), dpi=100)
	a = fig.add_subplot(111)
	
	a.plot(window_array[:window_preview,0], window_array[:window_preview,1], 'k-', label='Trace')
	a.plot(window_array[int(av_window.get()):window_preview,0], window_array[int(av_window.get()):window_preview,2], 'r-', label='Baseline Average')

	#adding detected events
	print numpy.average(window_array[:window_preview,1])
	a.plot(event_array[:,0], [numpy.average(window_array[:window_preview,1]) for i in range(len(event_array[:,0]))], 'go', label='Detected Events')
	
	#adding minima and maxima
	a.plot(event_array[:,0],event_array[:,2] , 'ro', label='Min/Max')
	a.plot(event_array[:,0],event_array[:,5] , 'ro' )

	#adding FWHM data
	a.plot(event_array[:,10],[(event_array[k,8]-0.5*event_array[k,4]) for k in range(len(event_array[:,8]))] , 'm.', label='FWHM')
	a.plot(event_array[:,11],[(event_array[k,8]-0.5*event_array[k,4]) for k in range(len(event_array[:,8]))] , 'm.')

	a.plot(event_array[:,13],[(event_array[k,8]-0.5*event_array[k,7]) for k in range(len(event_array[:,8]))] , 'm.')
	a.plot(event_array[:,14],[(event_array[k,8]-0.5*event_array[k,7]) for k in range(len(event_array[:,8]))] , 'm.')


	threshold_value = float(threshold.get())
	#setting threshold values
	a.plot(window_array[int(av_window.get()):window_preview,0],[k*(1-(threshold_value/100)) for k in window_array[int(av_window.get()):window_preview,2]],'b-',label='Threshold')
	a.plot(window_array[int(av_window.get()):window_preview,0],[k*(1+(threshold_value/100)) for k in window_array[int(av_window.get()):window_preview,2]],'b-')

	a.set_ylim(min(event_array[:,5])*(1-(threshold_value/100)), max(event_array[:,2])*(1+(threshold_value/100)))
	a.set_xlim(0.0, max(window_array[:window_preview,0]))
	
	a.set_xlabel('Time (s)')
	a.set_ylabel('current (nA)')
	a.legend(loc='best')
	
	dataPlot = FigureCanvasTkAgg(fig, master=f5i)
	dataPlot.draw()	
	dataPlot.show()
	dataPlot.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)


def plot(window_array):
	#setting preview window parameter
	window_preview = int(preview_time.get()) #in no of data pointa at 50kHz

	import matplotlib.pyplot as plt
	from matplotlib.figure import Figure
	from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

	fig = Figure(figsize=(18,4), dpi=100)
	a = fig.add_subplot(111)
	
	a.plot(window_array[:window_preview,0], window_array[:window_preview,1], 'k-', label='Trace')
	a.plot(window_array[:window_preview,0], window_array[:window_preview,2], 'r-', label='Baseline Average')
	a.set_ylim(min(window_array[100:window_preview,1]), max(window_array[100:window_preview,1]))
	
	a.set_xlabel('Time (s)')
	a.set_ylabel('current (nA)')
	a.legend(loc='best')
	


	dataPlot = FigureCanvasTkAgg(fig, master=f5i)
	dataPlot.draw()	
	dataPlot.show()
	dataPlot.get_tk_widget().pack(side=TOP)
	

def sample_plot(sample_array):

	import matplotlib.pyplot as plt
	from matplotlib.figure import Figure
	from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

	fig_s = Figure(figsize=(18,4), dpi=100)
	a = fig_s.add_subplot(111)
	
	a.plot(sample_array[:,0], sample_array[:,1], 'k-', label='Trace')
	a.set_ylim(min(sample_array[:,1]), max(sample_array[:,1]))
	
	a.set_xlabel('Time (s)')
	a.set_ylabel('Current (nA)')
	a.legend(loc='best')
	
	dataPlot = FigureCanvasTkAgg(fig_s, master=f5i)
	dataPlot.draw()	
	dataPlot.show()
	dataPlot.get_tk_widget().pack(side=TOP)



#erase button

def doErase():
	global f5i
	#print dataPlot
	print f5i
	c=0
	for item in f5i.pack_slaves():
		if c==0:
			item.destroy()
			c=c+1
#setting up file menu

menubar = Menu(window)
filemenu = Menu(menubar, tearoff=0)
filemenu.add_command(label='Open', command=doOpen, accelerator="Ctrl+O")
menubar.add_cascade(label='File', menu=filemenu)
window.config(menu=menubar)

#frame 0
f0 = Frame(window)
inpt_fln_L = Label(f0, text='Input (.irbd) File location')
inpt_fln = StringVar()
inpt = Entry(f0,textvariable=inpt_fln)
inpt_fln_L.pack(side=LEFT)
inpt.pack(side=LEFT)

samp = Button(f0, text='Sample File', command=doSample) 
samp.pack(side=RIGHT,padx=10)

convert = Button(f0, text='Convert', command=doConvert) 
convert.pack(side=RIGHT,padx=10)

f0.pack()

#frame 0a
f0a = Frame(window)
preview_time_L = Label(f0a, text='Preview time window')
preview_time = StringVar()
preview_time.set('25000')
prv = Entry(f0a,textvariable=preview_time)
preview_time_L.pack(side=LEFT)
prv.pack(side=LEFT)
f0a.pack()

#frame 0b
f0b = Frame(window)
split_L = Label(f0b, text='Split file size (in MB)')
split_size = StringVar()
split_size.set('500')
split_bs = Entry(f0b,textvariable=split_size)
split_L.pack(side=LEFT)
split_bs.pack(side=LEFT)

split = Button(f0b, text='Split', command=doSplit) 
split.pack(side=RIGHT,padx=10)


f0b.pack()

#frame 1a
f1a = Frame(window)
av_L = Label(f1a, text='Baseline Averaging Window (ms)')
av_window = StringVar()
av_window.set('1000')
av_bs = Entry(f1a,textvariable=av_window)
av_L.pack(side=LEFT)
av_bs.pack(side=LEFT)
f1a.pack()

#frame 1b
f1b = Frame(window)
inpt_dat_L = Label(f1b, text='.dat File location')
inpt_dat = StringVar()
dat = Entry(f1b,textvariable=inpt_dat)
inpt_dat_L.pack(side=LEFT)
dat.pack(side=LEFT)
f1b.pack()


#frame 1
f1 = Frame(window)
win_L = Label(f1, text='Estimate event duration (ms)')
event_window = StringVar()
event_window.set('10')
e = Entry(f1,textvariable=event_window)
win_L.pack(side=LEFT)
e.pack(side=LEFT)

attempt_FWHM = IntVar()
FWHM_d = Checkbutton(f1, text="Attempt FWHM analysis?", variable=attempt_FWHM)
FWHM_d.pack(side=RIGHT)

f1.pack()

#frame 2
f2 = Frame(window)
thr_L = Label(f2, text='Threshold Value (% of baseline)')
threshold = StringVar()
threshold.set('0.04')
thr = Entry(f2,textvariable=threshold)
thr_L.pack(side=LEFT)
thr.pack(side=LEFT)

discard = IntVar()
cb_d = Checkbutton(f2, text="Discard multi-peak events?", variable=discard)
cb_d.pack(side=RIGHT)

f2.pack()

#frame 3
f3 = Frame(window)
fln_L = Label(f3, text='Output filename')
filename = StringVar()
filename.set('Output.csv')
fln = Entry(f3,textvariable=filename)
fln_L.pack(side=LEFT)
fln.pack(side=LEFT)

event_csvs = IntVar()
cb = Checkbutton(f3, text="Output events to csvs?", variable=event_csvs)
cb.pack(side=RIGHT)

f3.pack()

#frame 4
f4 = Frame(window)

go = Button(f4, text='Analyse', command=doRead) 
go.pack(side=LEFT,padx=10)

rm = Button(f4, text='Erase Plot', command=doErase) 
rm.pack(side=LEFT,padx=10)

f4.pack()

#frame 5

f5 = Frame(window)
f5i = Frame(f5)
f5i.pack()
f5.pack()

# packing default plot into f5i
#import matplotlib.pyplot as plt
#from matplotlib.figure import Figure
#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

#fig = Figure(figsize=(5,4), dpi=100)

#dataPlot = FigureCanvasTkAgg(fig, master=f5i)
#dataPlot.get_tk_widget().pack(side=TOP)
#dataPlot.show()

window.mainloop() 
