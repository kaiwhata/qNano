#!/usr/bin/python

##truncated cone particle resistance calcuation methods
import math
from scipy.integrate import odeint
from scipy.integrate import quad
from scipy import real
import numpy as np
import matplotlib.pyplot as plt
from sympy import limit
from pylab import twinx
import sys

"""
Pore parameters
                a0
           |----------|           _      _
-----------\          /-----------|      |        
            \    a   /            | p    |     
             \|----|/             _      |      
             /      \                    |             
            /        \                   | d       
           /          \                  |       
          /            \                 |     
         /              \                |   
        /                \               |       
-------/                  \-------       _
       |------------------|
                 b
"""

def find_nearest(array,value):
	    idx=(np.abs(array-value)).argmin()
	    return idx


"""
####DEFINE Parameters
"""
#SET PARAMETERS FOR SPECIFIC RUN 
V0 = 0.5 #V 
kB = 1.38e-23 #Joules per Kelvin
Temp = 293 #Kelvin
permettivity = 7.1e-10 #Coulomb^2 N^-1 m^-1
viscosity = 0.001002 #eta #Pa s; Reference 408; my p. 325, 337; 15C = 1140 mPa s; 25C = 894 mPa s
conductivity = 0.75 #rho #CHECK THIS FOR 0.1M KCL with 0.01M Tris/HEPBS

P0 = 50 #pa coressponded to 4.7mm pressure head of water. 9.81 Pa per 1 mm water; 4.7 mm is the likely 'zero' pressure head

#-----
#scaling for electric field 'cone' outside pore ... equals 1.25 to get the
#'0.8a' end effects factor
ends=1.25

'''Effective charge per surface area--------POORLY KNOWN PARAMETER-------
# Ref 780 For 100 nm carboxylated latex nanospheres: "effective charges
# measured by static light scattering and from electrophoretic light 
# scattering do not exceed a value of 500 e for carboxylate latex"
# ~reasonably consistent with Stober et al. 844 suggest -2.6e-3 for 1.1 micron radius spheres in NaCl'''
EffCh=-1.54e-3; #Coulombs per metres sq; negative means negative surface charge

'''Pore wall surface potential--------POORLY KNOWN PARAMETER-------
 Not used in this model at present: assumption that electroosmosis is not
 significant to transport or current'''
psi0=0.075; #Volt

#DEFINE PORE GEOMETRY

a=200.0e-9  #pore constriction radius (m) 
b=5.0e-6	#large pore opening radius (m)
d=160.0e-6  #pore length (m)
p =5.0e-6 #depth of constriction beneath a0 (m)
a0 = a+p*(b-a)/(d-p) #radius of small opening (m)
print 'small pore opening = '+str(a0/1e-9)+'nm'

	##Check baseline current calculations
	
I0 = 130e-9 #input from experiment 	
print 'I0 = '+str(I0)+'nA'
Icheck=(math.pi*V0/conductivity)*(a0*a*b/(p*b+(0.8*a*(a0+b))+a0*(d-p)))
print 'I calculated = '+str(Icheck)+'nA'

R0 = V0/I0
print 'R0 = '+str(R0)+' Ohms'

##DEFINE PARTICLE
#Sphere
r=100e-9 #sphere radius (nm) 
DiffCoeff=kB*Temp/(6*math.pi*viscosity*r)
SurfArea=4*math.pi*(r**2)

##idiot check for paticle size compared to pore
if r >= a:
	print "Particle is larger than pore opening........Jackass"
	sys.exit()
else:
	print 'Paticle fits'
	
#DEFINE PARTICLE LOCATION IN PORE
'''The position of the particle, off centre, which would be a function of
height %--------TO DEVELOP/STUDY------------------
Present assumptions are such that this affects pressure-driven transport but not
resistance.'''
rCent = 0


##DEFINE INITIAL CONDITIONS

#limit for integral (far away from pore), should be infinity "inf"
lim=np.inf

#Note convergence errors if edge of particle is coincident with a boundary
init=1e-6 #initial distance of the pore above upper surface (m)
tlim= 0.01 # length of time for simulation
steplength=2e-6 #time resolution (s)

tspan = np.linspace(0.0, tlim, (tlim/steplength)) #This determines timespan of generated data

zinit = d+init #Initial position of the particle


"""
######DEFINE calculation methods
"""

#ELECTRIC FIELD WITHIN PORE
def E(z): #Electric field in the pore
	return I0*conductivity/(math.pi*(rPore(z)**2))
	

##### PORE geometry
def rPore(z):
	if z>d: #above pore
		r = a0+16*(z-d)/(27*math.pi)
	elif z>=(d-p):#within top cone	
		r = a + ((b-a)/(d-p))*(z-(d-p))
	elif z>=0: #within bottom cone
		r = b-((b-a)/(d-p))*z
	elif z<0: #below pore
		r = b-16*(z)/(27*math.pi)
	return r

## idiot check that pore radius calculations are correct
print 'Pore geometry checks d and 0 respecively'
if not abs(rPore(0.0)-b) <= 1e-9:
	print "Pore geometry is calculating incorrectly at large opening"
	print 'rPore(0) ='+str(rPore(0.0))
	print 'b ='+str(b)
	sys.exit()
else:
	print 'Pore geometry correct at base'

if not abs(rPore((d-p))-a) <= 1e-9:
	print "Pore geometry is calculating incorrectly at constriction"
	print 'rPore(d) ='+str(rPore(d-p))
	print 'a ='+str(a)
	sys.exit()
else:
	print 'Pore geometry correct at constriction'

####Flow profile
def pressure_profile(z):  #integral for pressure with general pore profile
    return 8.*viscosity/(math.pi*rPore(z)**4)

#calculate flow profile U
def U(): #combining end effects with integral of pressures within pore
	'''Assumes: 
	(i)pressure drop outside ends equivalent to Sampson result ref 916
	Would be better to use the semi-infinite solution:
	Wang-Yi and Skalak, Applied Mathematics and Mechanics 6 #1 1985 p. 9
	(ii)Pouiseille flow (parabolic profile) internal to pore; ~OK because
	 pore radius changes slowly with pore length'''
	return P0/((9*viscosity/(2*(a0**3)))+(9*viscosity/(2*b**3))+quad(pressure_profile, 0., d)[0])

print 'Pressure_profile integration from 0 to d'
print quad(pressure_profile, 0., d)[0]
print 'U integration and evaluation test'
print U()

#####DEFINE specific flow regimes

def pressure_flow_only(z, t): #Four terms: diff, eph, eos, pdf
	
	termdiff= 0 #Ensure sign is correct
	termeos= 0
	termeph= 0
	
	if z > d: #particle above 'a' end of pore
		termpdf= -2*U()*((rPore(z)**2-rCent**2)/(math.pi*rPore(z)**4))
	elif z < 0: #particle below pore
		termpdf = -2*U()*((rPore(z)**2-rCent**2)/(math.pi*rPore(z)**4)) 
	elif z <= d: #particle within pore
		termpdf = -2*U()*((rPore(z)**2-rCent**2)/(math.pi*rPore(z)**4))	
	else:
		print "Where the fuck is the particle then??"
	
	dzdt = termdiff + termeph + termeos + termpdf;
	return dzdt

def eo_flow_only(z, t): #Four terms: diff, eph, eos, pdf
	#Function describing electroosmotic flow
	A=1
	
	termdiff= 0 #Ensure sign is correct
	termeph= 0
	termpdf= 0
	
	if z > d: #particle above 'a' end of pore
		termeos = 0
	elif z < 0: #particle below pore
		termeos = 0
	elif z <= d: #particle within pore
		termeos= -permettivity*psi0*A*E(z)/(4*math.pi*viscosity) #Sign correct?	
	else:
		print "Where the fuck is the particle then??"
	
	dzdt = termdiff + termeph + termeos + termpdf;
	return dzdt

def eph_flow_only(z, t): #Four terms: diff, eph, eos, pdf
	#Function describing electroosmotic flow
	A=1
	
	termdiff= 0 #Ensure sign is correct
	termeos= 0
	termeph= 0
	termpdf= 0

	if z > d: #particle above 'a' end of pore
		termeph= (SurfArea*EffCh*DiffCoeff/(kB*Temp))*(I0*conductivity/(math.pi*((a0+ends*(z-d))**2)))#Sign correct
	elif z < 0: #particle below pore
		termeph= (SurfArea*EffCh*DiffCoeff/(kB*Temp))*(I0*conductivity/(math.pi*((b-ends*z)**2)))#Sign correct
	elif z <= d: #particle within pore
		termeph= SurfArea*EffCh*DiffCoeff*E(z)/(kB*Temp)#Sign correct				
	else:
		print "Where the fuck is the particle then??"
	
	dzdt = termdiff + termeph + termeos + termpdf;
	return dzdt

def ek_flow_only(z, t): #Four terms: diff, eph, eos, pdf
	#Original name from m script: dzdt = zfnt(t,z)
	
	#Function describing electroosmotic flow
	A=1
				
	termdiff= 0 
	termeos=0
	termeph= 0
	termpdf= 0
	
	if z > d: #particle above 'a' end of pore
		termeph= (SurfArea*EffCh*DiffCoeff/(kB*Temp))*(I0*conductivity/(math.pi*((a0+ends*(z-d))**2)))#Sign correct
		termeos=0
	elif z < 0: #particle below pore
		termeph= (SurfArea*EffCh*DiffCoeff/(kB*Temp))*(I0*conductivity/(math.pi*((b-ends*z)**2)))#Sign correct
		termeos=0
	elif z <= d: #particle within pore
		termeph= SurfArea*EffCh*DiffCoeff*E(z)/(kB*Temp)#Sign correct
		termeos= -permettivity*psi0*A*E(z)/(4*math.pi*viscosity) #Sign correct?				
	else:
		print "Where the fuck is the particle then??"
	
	dzdt = termdiff + termeph + termeos + termpdf;
	return dzdt

#####DEFINE Summative regime
def dzdt(z,t): #Four terms: diff, eph, eos, pdf
	#Original name from m script: dzdt = zfnt(t,z)
	
	#Function describing electroosmotic flow
	A=1
				
	termdiff= 0 #Ensure sign is correct
	termeos= 0
	termeph= 0
	
	##particle moving from a to b side assumed.....but not important(??)
	if z > d: #particle above 'a' end of pore
		#print z
		#print 'particle above pore'
		termpdf= -2*U()*((rPore(z)**2-rCent**2)/(math.pi*rPore(z)**4))
		termeph= (SurfArea*EffCh*DiffCoeff/(kB*Temp))*(I0*conductivity/(math.pi*((a0+ends*(z-d))**2)))#Sign correct
	elif z < 0: #particle below pore
		#print z
		#print 'particle below pore'
		termpdf = -2*U()*((rPore(z)**2-rCent**2)/(math.pi*rPore(z)**4)) 
		termeph= (SurfArea*EffCh*DiffCoeff/(kB*Temp))*(I0*conductivity/(math.pi*((b-ends*z)**2)))#Sign correct
	elif z <= d: #particle within pore
		#print z
		#print 'particle within pore'
		termpdf = -2*U()*((rPore(z)**2-rCent**2)/(math.pi*rPore(z)**4))# P0 is overpressure from top to bottom page 782
		termeph= SurfArea*EffCh*DiffCoeff*E(z)/(kB*Temp)#Sign correct
		termeos= -permettivity*psi0*A*E(z)/(4*math.pi*viscosity) #Sign correct?				
	else:
		print "Where the fuck is the particle then??"
	
	dzdt = termdiff + termeph + termeos + termpdf;
	return dzdt

###DEFINE SOLUTIONS DEPENDING ON REGIME

def pore_with_particle(z): #numerical solution within the pore, particle present #f1
	return 1./(rPore(z)**2-r**2)

def pore_no_particle(z): #numerical solution within the pore, no particle #f2
	return 1./(rPore(z)**2)

def above_pore_with_particle(z): #numerical solution above the pore, particle present #f3a
	return 1./((a0+ends*(z-d))**2-r**2)

def below_pore_with_particle(z):  #numerical solution below the pore, particle present #f3b
	return 1./((b-ends*z)**2-r**2)

def above_pore_no_particle(z): #analytic solution above the pore, no particle #f4a
	return -1./(ends*(a0+ends*(z-d)))

def below_pore_no_particle(z): #analytic solution below the pore, no particle #f4b
	return 1./(ends*(b-ends*z))

#### define Runge-Kutta 4th order numerical integration method

def rk4( f, x0, t ):
    """Fourth-order Runge-Kutta method to solve x' = f(x,t) with x(t[0]) = x0.

    USAGE:
        x = rk4(f, x0, t)

    INPUT:
        f     - function of x and t equal to dx/dt.  x may be multivalued,
                in which case it should a list or a NumPy array.  In this
                case f must return a NumPy array with the same dimension
                as x.
        x0    - the initial condition(s).  Specifies the value of x when
                t = t[0].  Can be either a scalar or a list or NumPy array
                if a system of equations is being solved.
        t     - list or NumPy array of t values to compute solution at.
                t[0] is the the initial condition point, and the difference
                h=t[i+1]-t[i] determines the step size h.

    OUTPUT:
        x     - NumPy array containing solution values corresponding to each
                entry in t array.  If a system is being solved, x will be
                an array of arrays.
    """

    n = len( t )
    x = np.array( [ x0 ] * n )
    for i in xrange( n - 1 ):
        h = t[i+1] - t[i]
        k1 = h * f( x[i], t[i] )
        k2 = h * f( x[i] + 0.5 * k1, t[i] + 0.5 * h )
        k3 = h * f( x[i] + 0.5 * k2, t[i] + 0.5 * h )
        k4 = h * f( x[i] + k3, t[i+1] )
        x[i+1] = x[i] + ( k1 + 2.0 * ( k2 + k3 ) + k4 ) / 6.0

    return x


"""
Testing that variables act as they should in the right directions etc
"""
#FLOW TESTS 
expt_z = np.linspace(-10e-6,d+10e-6,1000) #dummy data
expt_z_microns = expt_z/1e-6
#total flow
dzdt_points = []
for i in expt_z:
	dzdt_points.append(dzdt(i,0))
#radius
radius_points = []
for i in expt_z:
	radius_points.append(rPore(i)/1e-6)
#Pressure Flow Only
pressure_points = []
for i in expt_z:
	pressure_points.append(pressure_flow_only(i,0))
#Eo Flow Only
Eo_points = []
for i in expt_z:
	Eo_points.append(eo_flow_only(i,0))
#Eph Flow Only
Eph_points = []
for i in expt_z:
	Eph_points.append(eph_flow_only(i,0))
#Ek Flow Only
Ek_points = []
for i in expt_z:
	Ek_points.append(ek_flow_only(i,0))


"""
Calculations
"""


#python ode solver format = odeint(function, i.c., points)
#solution = odeint(dzdt, zinit, tspan)
solution = rk4(dzdt, zinit, tspan) #original = dzdt
print 'Length of solution = '+str(len(solution))
print 'Final solution = '+str(solution[-1])
print 'Final time = '+str(tspan[-1])


#The resistance partitioned depending on sphere location
#q = quad(fun,a,b) tries to approximate the integral of function fun from a to b 
#to within an error of 1e-6 using recursive adaptive Simpson quadrature. 

#data output with header
output_format = np.array([[zinit,0.0, Icheck, 0.0]])

#DEFINING PARTICLE REGIMES

for i in solution:
	h = i # h is the current instatiation of z
	#quad_options = 'epsabs = 1e-7, epsrel = 1e-4 '
	#'InitialStep',steplength,'MaxStep',steplength
	
	if (h-r) >= d: #particle above pore
		Within=(conductivity/math.pi)*quad(pore_no_particle,0,d,epsabs = 1e-9)[0] #numerical within pore no particle
		Below_particle=(conductivity/math.pi)*(above_pore_no_particle(h-r)-above_pore_no_particle(d)) #analytic from top of membrane to base of particle
		Particle=(conductivity/math.pi)*quad(above_pore_with_particle,h-r,h+r,epsabs = 1e-9)[0] #numerical above pore particle present
		Above_particle=(conductivity/math.pi)*(above_pore_no_particle(lim)-above_pore_no_particle(h+r)) #analytic from top of particle to infinity 
		Below = (conductivity/math.pi)*(below_pore_no_particle(0)-below_pore_no_particle(-lim)) #analytic below pore no particle
		R = (Above_particle+Below_particle+Particle+Below+Within) 
		#print 'particle above pore'
				
	elif (h-r) >= (d-2*r): #transition between above to within pore
		Within_to_particle = quad(pore_no_particle,0,h-r,epsabs = 1e-9)[0] #numerical from base of membrane to bottom of particle
		Particle_to_top = quad(pore_with_particle,h-r,d,epsabs = 1e-9)[0] #numerical from base of particle to top of membrane
		Above_particle = quad(above_pore_with_particle,d,h+r,epsabs = 1e-9)[0] #numerical from top of membrane to top of particle
		Above = (above_pore_no_particle(lim)-above_pore_no_particle(h+r)) #analytic from top of particle to infinity
		Below = (below_pore_no_particle(0)-below_pore_no_particle(-lim)) #analytic from infinity to base of membrane without particle
   		R = (conductivity/math.pi)*(real(Within_to_particle+Particle_to_top+Above_particle)+Above+Below)  
		#print 'transition between above to within'

	elif (h-r) >= 0: #particle within pore
		Above_particle=(conductivity/math.pi)*quad(pore_no_particle,0,h-r,epsabs = 1e-9)[0] #numerical within pore no particle
		Particle=(conductivity/math.pi)*quad(pore_with_particle,h-r,h+r,epsabs = 1e-9)[0] #numerical within pore particle present (over area containing particle)
		Below_particle=(conductivity/math.pi)*quad(pore_no_particle,h+r,d,epsabs = 1e-9)[0] #numerical within pore no particle
		Above_and_below=(conductivity/math.pi)*((above_pore_no_particle(lim)-above_pore_no_particle(d))+(below_pore_no_particle(0)-below_pore_no_particle(-lim))) #analytic above and below no particle
		R = (Above_particle+Particle+Below_particle+Above_and_below)
		#print 'particle within pore'
	
	elif (h-r) > (0-2*r): #transition between below to within pore
		Within_to_particle = quad(pore_no_particle,0.,h+r)[0] #numerical from base to upper edge of particle
		Above_particle = quad(pore_no_particle,h+r,d)[0] #numerical from top of particle to top of pore
		Particle_to_bottom = quad(below_pore_with_particle,h-r,0)[0] #numerical from membrane base to lower edge of particle
		Above = above_pore_no_particle(lim)-above_pore_no_particle(d) #analytic above pore without particle
		Below = below_pore_no_particle(h-r)-below_pore_no_particle(-lim) #  analytic from base of particle to infinity
		R = (conductivity/math.pi)*(real(Within_to_particle+Above_particle+Particle_to_bottom)+(Above)+(Below))
		#print 'transition between below to within'
		
	elif (h+r) < 0: #particle below pore
		Below = (below_pore_no_particle(h-r)-below_pore_no_particle(-lim)) #analytic below particle to ininity
		Particle = quad(below_pore_with_particle,h-r,h+r)[0] #numerical below pore over diameter of particle
		Below_to_particle = (below_pore_no_particle(0)-below_pore_no_particle(h+r)) #analytic from membrane base to top of particle
		Within = quad(pore_no_particle,0,d)[0] #numerical solution within pore no particle
		Above = (above_pore_no_particle(lim)-above_pore_no_particle(d)) #analytic above pore no particle
		R = (conductivity/math.pi)*real(Below+Particle+Below_to_particle+Within+Above)
		#print 'particle below pore'
	else:
		print "Shit i've lost the particle!"
		sys.exit()
	I = V0/R
	dI = I-Icheck
	output_format = np.append(output_format,[[h,R,I, dI]],axis=0)
	
	#discards first row
output_format = np.delete(output_format, 0, 0)

print 'Solution Matrix'
print output_format.shape
R_list = output_format[:,1]
z_list = output_format[:,0]
I_list = output_format[:,2]
Delta_I = output_format[:,3]

"""
Plotting
"""
j = 1 #j= variable to keep track of figure numbers

def flow_overlaid(j):
	#Pore geometry and total flow velocity check
	plt.figure(j)
	ax1 = plt.subplot(111)
	one = ax1.plot(expt_z_microns,dzdt_points)
	two = ax1.plot(expt_z_microns,pressure_points)
	three = ax1.plot(expt_z_microns,Eo_points)
	four = ax1.plot(expt_z_microns,Eph_points)
	five = ax1.plot(expt_z_microns,Ek_points)
	ax1.legend([one,two,three,four,five],['Total','P @ '+str(P0),'Eo','Eph','Ek @ '+str(V0)], loc=6)
	ax1.set_xlabel('z (micron)')
	ax1.set_ylabel('Flow Velocity, dz/dt (m/s)')
	
	ax2 = twinx()
	ax2.plot(expt_z_microns,radius_points, 'k-')
	ax2.axvline(x=0, ymin=0.01, ymax=0.99, color='b',linestyle=':')
	ax2.axvline(x=(d/1e-6), ymin=0.01, ymax=0.99, color='b',linestyle=':')
	ax2.set_ylabel('Radius (micron)')
	
	plt.title('Total Flow Velocity with z and pore geometry')
		 
flow_overlaid(j)
j=j+1

def flow_plot(ax1,expt_points,data_points, regime):
	plt.figure(j)
	ax1.plot(expt_points,data_points,'k-')
	#ax1.set_xlabel('z (microm)')
	ax1.set_ylabel(regime)
	ax2 = twinx()
	ax2.plot(expt_points,radius_points, 'r-')
	ax2.axvline(x=0, ymin=0.01, ymax=0.99, color='b',linestyle=':')
	ax2.axvline(x=(d/1e-6), ymin=0.01, ymax=0.99, color='b',linestyle=':')
	plt.show

def constituents(j):	
	#Constituent flow characteristics plot overlaid on radius
	plt.figure(j)
	plt.title('Flow Velocities with z and pore geometry')
	plt.ylabel('Flow Velocity')
	plt.xlabel('z (um)')
	#pressure
	ax3 = plt.subplot(411)
	flow_plot(ax3,expt_z_microns,pressure_points,'P')
	#Eo
	ax4 = plt.subplot(412)
	flow_plot(ax4,expt_z_microns,Eo_points, 'Eo')
	#Eph
	ax5 = plt.subplot(413)
	flow_plot(ax5,expt_z_microns,Eph_points, 'Eph')
	#Ek
	ax6 = plt.subplot(414)
	flow_plot(ax6,expt_z_microns,Ek_points, 'Ek')
	plt.show()

def z_solution(j):
	plt.figure(j)
	
	plt.subplot(211)
	plot_tspan = tspan/1e-6
	plot_solution = solution/1e-6
	 
	plt.plot(plot_tspan,plot_solution, 'g-')  
	plt.axhline(y=d/1e-6, xmin=0.01, xmax=0.99, color='r',linestyle=':')
	               

	plt.ylabel('Position(z) (um)'); 
	plt.title('Z and Trace with Time, dz/dt numerical solution')
	######zplot finished
	plt.subplot(212)
	
	plot_I_list = Delta_I/1e-9
	plot_t_list = tspan/1e-6
	
	plt.plot(plot_t_list, plot_I_list, 'k-')
	plt.ylabel('Current (nA)')
	plt.xlabel('Time (ms)'); 	
	#Find maximum Delta I and FWHM of what was just calcuated
	max_Delta_I = min(Delta_I)
	Delta_I_list = Delta_I.tolist()
	index_max = Delta_I_list.index(max_Delta_I)
	tAtMax=tspan[index_max]
	
	plt.axvline(x=tAtMax/1e-6, ymin=0.01, ymax=0.99, color='b',linestyle=':')  
	
	##split list at max
	below_max = output_format[:index_max,3]
	above_max = output_format[index_max:,3]
	
	FWHM1 = tspan[find_nearest(below_max,max_Delta_I/2)]
	FWHM2 = tspan[len(below_max)+find_nearest(above_max,max_Delta_I/2)]
	
	plt.axvline(x=FWHM1/1e-6, ymin=0.01, ymax=0.99, color='r',linestyle=':')
	plt.axvline(x=FWHM2/1e-6, ymin=0.01, ymax=0.99, color='r',linestyle=':')    
	
	#calculate event asymmetry - as defined by Willmott and Parry 2011
	F = ((tAtMax-FWHM1)/(FWHM2-FWHM1))
	
	plottext = "Asymmetry = "+str(F)
	plt.text(0.5*max(tspan/1e-6),0.6*min(plot_I_list), plottext )
	
	print 'FWHM = '+str(FWHM2-FWHM1)+'s'
	print 'Asymmetry "F" = '+str(F)
	return 
	
z_solution(j)
j=j+1

def resistance_plot(j):
	plt.figure(j)
	plot_R_list = R_list/1e-6
	plot_z_list = z_list/1e-6
	plt.plot(plot_z_list, plot_R_list)
	plt.axvline(x=(d/1e-6), ymin=0.01, ymax=0.99, color='b',linestyle=':')
	plt.xlabel('Z (um)')
	plt.ylabel('Resistance (Ohms)')


resistance_plot(j)
j=j+1


plt.show()

#ADDITIONAL NOTES FOR MODEL
'''
%End effects at the small end are crucial.
%If geometry is to change, then rPore, Ez and flow terms should be checked.

% The present general approach is as follows:
% The ELECTRIC FIELD assumes uniform resistivity.
% - Resistance within the pore is determined by the geometry of the pore.
% - Resistance beyond each end is well-known: proportional to 0.8 times
% the opening radius (e.g. Jeans paper; deBlois and Bean (1970); Hall 913; Stober 844)
% - the electric field beyond the membrane is modelled as a cone with pitch so that the resistance is correct 
% - this method is imperfect but more practical than others; E(z) does not depend on off-axis geometry
% - The blockade size is determined purely geometrically by integrating the electric field flux through the unblocked area along z.
% The PRESSURE DRIVEN FLOW beyond the membrane also assumes an artificial cone.
% - The extrenal pressure benchmark is due to Sampson 916; this might be improved using the semi-infinite cylinder solution: Wang-Yi and Skalak, Applied Mathematics and Mechanics 6 #1 1985 p. 9
% - The artificial cones at the cylinder ends have a pitch such that the external pressure is applied between the membrane surface and (infinity).
% - Both internal and external to the pore, planar Pouseuille flow is assumed, parabolic across the pore / cone. This assumption is OK as long as radii are slowly varying with z.
% - The pressure gradient is unique at any value of z, and allows calculation of vz(r) at that point.
'''



