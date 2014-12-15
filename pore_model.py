def find_nearest(array,value):
	    import numpy as np
	    idx=(np.abs(array-value)).argmin()
	    return idx

def pore_calculate(a,b,d,rxy,rz,P0,V0, psi0, theta_d, zeta_part):
	
	import math
	from scipy.integrate import odeint
	from scipy.integrate import quad
	from scipy import real
	import numpy as np
	#import matplotlib.pyplot as plt
	#from sympy import limit
	from pylab import twinx
	import sys
	
	###OBLATE/PROLATE Particle Characteristics
		
	def particle_eccentricity(rz,rxy): #SOURCE:wikipedia
		if rz<rxy:#spheroid
			print 'Oblate Spheroid'
			e = np.sqrt(1-((rz**2)/(rxy**2))) 
		elif rxy<rz:#prolate
			print 'Prolate Spheroid'
			e = np.sqrt(1-((rxy**2)/(rz**2))) 
		else:
			print 'Sphere'
			e = 0 #sphere
		return e
		
	def particle_surface_area(rz,rxy):  #SOURCE:wikipedia
		E = particle_eccentricity(rz,rxy)
		if rz<rxy: #oblate
			#print 'Oblate Spheroid'
			SurfArea = 2*math.pi*rxy**2*(1+((1-E**2)/E)*np.arctanh(E))
		elif rxy<rz:#prolate
			#print 'Prolate Spheroid'
			SurfArea = 2*math.pi*rxy**2*(1+(rz/(rxy*E))*np.arcsin(E)) 
		else:#sphere
			 #print 'Sphere'
			 SurfArea= 4*math.pi*(rz**2) 
		return SurfArea

	def particle_diffusion(rz, rxy):     #SOURCE:Happel and Brenner
		kB = 1.38e-23 #Joules per Kelvin
		Temp = 293 #Kelvin
		viscosity = 0.001002 #eta #Pa s; Reference 408; my p. 325, 337; 15C = 1140 mPa s; 25C = 894 mPa s

		#Elliptical particles: diffusion coefficient and surface area
		#For DiffCoeff see Happel and Brenner 5.11.18 and 5.11.20
		#For SurfArea see wikipedia
		#Each of these has an 'oblate' vs 'prolate' case 
		v = rz/rxy #ratio
		print v
		if (rxy > rz): #Oblate
			#print 'Oblate Spheroid'
			DiffCoeff = kB*Temp/(6*math.pi*viscosity*np.real((8*rxy/3)*(1/(-2*(v)/((v)**2-1)+((2*((v)**2)-1)/(((v)**2-1+0j)**1.5))*np.log(((v)+((v)**2-1+0j)**0.5)/((v)-((v)**2-1+0j)**0.5)))))) 
		elif(rxy < rz): #Prolate
			#print 'Prolate Spheroid'
			DiffCoeff = kB*Temp/(6*math.pi*viscosity*np.real((8*rxy/3)*(1/(2*(v)/(1-(v)**2)+(2*(1-2*((v)**2))/((1-(v)**2+0j)**1.5))*np.arctan(((1-(v)**2+0j)**0.5)/(v)))))) 
		else: #Sphere
			#print 'Sphere'
			DiffCoeff=kB*Temp/(6*math.pi*viscosity*rz)	
		#Happel and Brenner table 5-11.1 p.223
		#For oblate ellipsoids with length/diameter 0.1 to 0.5, the reduction in drag is 71.9% to 87.6% when turned on its side.
		#The change in drag due to the walls could increase by up to a half in the orifice on-axis. Worse off-axis.
		#See table 7-6.1 p. 341; table 7-3.1 p. 309; Eq 7-2.15 p. 291
		return DiffCoeff

	def particle_volume(rz, rxy):     #SOURCE:http://en.wikipedia.org/wiki/Spheroid
		Vol = (4*math.pi*rxy**2*rz)/3 
		print 'Particle Volume is %s m^3' % Vol
		return Vol
    
	def Eff_ch(zeta_particle, permettivity, rz, rxy):
        #commented lines are Elf's incorrect Eff_ch calculation
		#DiffCoeff=kB*Temp/(6*math.pi*viscosity*radius)
		#eff_ch = (4*permettivity*math.pi*zeta_particle*1e-3*radius)/(DiffCoeff)
        ###############
		#def particle volume - use volume to calculate effective radius
		eff_radius = math.pow((3*particle_volume(rz, rxy))/(4*math.pi), 1.0/3)      
		eff_ch = 3*permettivity*zeta_particle/(2*eff_radius**2) #Geoff's edits
		print 'Eff Charge is %s C/m^2' % eff_ch
		return eff_ch	
		
	##define surface function of spheroid
	def rSpheroid(z,h,theta,rz,rxy): #equivalent radius of an angled ellipsoid 
		Cz = z-h
		Cx = 2*Cz*np.cos(theta)*np.sin(theta)*(rxy**2-rz**2)/(2*rxy**2*np.sin(theta)**2+2*rz**2*np.cos(theta)**2)
		EA=rxy**2*rz**2/rxy**2
		EB=(rxy**2.*np.sin(theta)**2)+(rz**2.*np.cos(theta)**2)
		EC1=(rxy**2.*rz**2)
		EC2=+((rxy**2-rz**2)*(2.*Cx*Cz*np.sin(theta)*np.cos(theta)))
		EC3=-1*Cz**2.*(rz**2.*np.sin(theta)**2+rxy**2.*np.cos(theta)**2)
		EC4=-1*(rxy**2.*np.sin(theta)**2+rz**2.*np.cos(theta)**2)*Cx**2
		EC=EC1+EC2+EC3+EC4
		return np.sqrt(EC/np.sqrt(EA*EB))
	
	###define vertical projection	
	def rProjetion(theta,rz,rxy): #calculates vertical length of particle at angle theta
		#import numpy as np
		r_p = rxy*np.sin(theta)**2+rz*np.cos(theta)**2
		return r_p	
		
	########
	def E(z): #ELECTRIC FIELD 
		if z > d:#Above pore
			e = I0*conductivity/(math.pi*((a+ends*(z-d))**2))
		elif z < 0:#Below pore
			e = I0*conductivity/(math.pi*((b-ends*z)**2))
		else:#Electric field in the pore
			e = I0*conductivity/(math.pi*(rPore(z)**2))
		return e
	
	##### PORE geometry
	def rPore(z):
		if z>d: #above pore
			r = a+16*(z-d)/(9*math.pi)
		elif z>=0: #within pore
			r = b-((b-a)/d)*z
		elif z<0: #below pore
			r = b-16*(z)/(9*math.pi)
		return r
	
	####Flow profile
	def pressure_profile(z):  #integral for pressure with general pore profile
	    return 8.*viscosity/(math.pi*rPore(z)**4)
	
	#calculate flow profile U
	def U(): #combining end effects with integral of pressures within pore
        	#print "Integral P profile = %s"%(quad(pressure_profile, 0., d)[0])
        	q = P0/((3*viscosity/(2*(a**3)))+(3*viscosity/(2*b**3))+quad(pressure_profile, 0., d)[0])
        	#print "Complete pressure effect = %s"%(q)
		return q
    
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
			termeph= (SurfArea*EffCh*DiffCoeff/(kB*Temp))*E(z)#Sign correct
		elif z < 0: #particle below pore
			#print z
			#print 'particle below pore'
			termpdf = -2*U()*((rPore(z)**2-rCent**2)/(math.pi*rPore(z)**4)) 
			termeph= (SurfArea*EffCh*DiffCoeff/(kB*Temp))*E(z)#Sign correct
		elif z <= d: #particle within pore
			#print z
			#print 'particle within pore'
			termpdf = -2*U()*((rPore(z)**2-rCent**2)/(math.pi*rPore(z)**4))# P0 is overpressure from top to bottom page 782
			termeph= SurfArea*EffCh*DiffCoeff*E(z)/(kB*Temp)#Sign correct
			termeos= -permettivity*psi0*A*E(z)/(4*math.pi*viscosity) #Sign correct?				
		else:
			print "Where the fuck is the particle then??"
		
		dzdt = termdiff + termeph + termeos + termpdf
		return dzdt
	
	###DEFINE SOLUTIONS DEPENDING ON REGIME

	def pore_with_particle(z): #numerical solution within the pore, particle present #f1
		return 1./(rPore(z)**2-rSpheroid(z,h,theta,rz,rxy)**2)

	def pore_no_particle(z): #numerical solution within the pore, no particle #f2
		return 1./(rPore(z)**2)

	def above_pore_with_particle(z): #numerical solution above the pore, particle present #f3a
		return 1./((a+ends*(z-d))**2-rSpheroid(z,h,theta,rz,rxy)**2)

	def below_pore_with_particle(z):  #numerical solution below the pore, particle present #f3b
		return 1./((b-ends*z)**2-rSpheroid(z,h,theta,rz,rxy)**2)

	def above_pore_no_particle(z): #analytic solution above the pore, no particle #f4a
		return -1./(ends*(a+ends*(z-d)))

	def below_pore_no_particle(z): #analytic solution below the pore, no particle #f4b
		return 1./(ends*(b-ends*z))
	
	def rk4( f, x0, t ):
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
	    
	def calculate_FWHM(tspan, Delta_I):
		#Find maximum Delta I and FWHM of what was just calcuated
		max_Delta_I = min(Delta_I)
		print max_Delta_I
		Delta_I_list = Delta_I.tolist()
		index_max = Delta_I_list.index(max_Delta_I)
		tAtMax=tspan[index_max]
		
		##split list at max
		below_max = Delta_I[:index_max]
		above_max = Delta_I[index_max:]
		
		FWHM1 = tspan[find_nearest(below_max,max_Delta_I/2)]
		print FWHM1
		FWHM2 = tspan[len(below_max)+find_nearest(above_max,max_Delta_I/2)]
		print FWHM2
				      
		#calculate event asymmetry - as defined by Willmott and Parry 2011
		F = ((tAtMax-FWHM1)/(FWHM2-FWHM1))
		return (FWHM2-FWHM1), F
    
    	def magn(Delta_I):
       		max_Delta_I = max(Delta_I)
        	min_Delta_I = min(Delta_I)
        	if abs(max_Delta_I)>=abs(min_Delta_I):
            		magn = abs(max_Delta_I)
        	else:
            		magn = abs(min_Delta_I)
        	return magn
    
	#SET PARAMETERS FOR SPECIFIC RUN 
	kB = 1.38e-23 #Joules per Kelvin
	Temp = 293 #Kelvin
	permettivity = 7.1e-10 #Coulomb^2 N^-1 m^-1
	viscosity = 0.001002 #eta #Pa s; Reference 408; my p. 325, 337; 15C = 1140 mPa s; 25C = 894 mPa s
	conductivity = 0.75 #rho #CHECK THIS FOR 0.1M KCL with 0.01M Tris/HEPBS
	ends=1.25
	#EffCh=-1.54e-3; #Coulombs per metres sq; negative means negative surface charge
	zeta_particle = zeta_part/1000 #mV
	
	
		##Check baseline current calculations
	R0 = (conductivity*(d+0.8*(a+b)))/(math.pi*a*b)#baseline pore resistance without particle or end effects	
	
	I0=V0/R0 #130e-9 
	
	Icheck=(math.pi*V0/conductivity)*(a*b/(0.8*(a+b)+d)) 
	
	#calculating general properties of oblate/spheroid particle
	SurfArea = particle_surface_area(rz,rxy)
	DiffCoeff = particle_diffusion(rz, rxy)
	EffCh=Eff_ch(zeta_particle, permettivity, rz, rxy)
	
	###for now leave theta=0
	#Angle of ellipsoid to z-axis
	theta = theta_d 
	#thetadeg=theta*180/math.pi; #In degrees
	
	##idiot check for paticle size compared to pore
	if max([rz,rxy]) >= a:
		print "Particle is larger than pore opening........Jackass"
		sys.exit()
	else:
		print 'Particle fits'
	
	#checking pressure calculations are correct
	print "Integral P profile = %s"%(quad(pressure_profile, 0., d)[0])
    	print "Complete pressure effect = %s"%(P0/((3*viscosity/(2*(a**3)))+(3*viscosity/(2*b**3))+quad(pressure_profile, 0., d)[0]))
    
	#DEFINE PARTICLE LOCATION IN PORE
	rCent = 0
	
	#calculate J at d
	J_at_d = dzdt(d,0)
	
	##DEFINE INITIAL CONDITIONS
	lim=np.inf
	init=10e-6 #initial distance of the pore above upper surface (m)
	tlim= 0.01 # length of time for simulation
	steplength=2e-6 #time resolution (s) - matches output of qnano
	tspan = np.linspace(0.0, tlim, (tlim/steplength)) #This determines timespan of generated data
	zinit = d+init #Initial position of the particle
	
	solution = rk4(dzdt, zinit, tspan) #original = dzdt
	
	#data output with header
	output_format = np.array([[zinit,0.0, Icheck, 0.0]])
	
	#DEFINING PARTICLE REGIMES
	for i in solution:
		h = i # h is the current instatiation of z
		#quad_options = 'epsabs = 1e-7, epsrel = 1e-4 '
		#'InitialStep',steplength,'MaxStep',steplength
		
		if (h-rProjetion(theta_d,rz,rxy)) >= d: #particle above pore
			Within=(conductivity/math.pi)*quad(pore_no_particle,0,d,epsabs = 1e-9)[0] #numerical within pore no particle
			Below_particle=(conductivity/math.pi)*(above_pore_no_particle(h-rProjetion(theta_d,rz,rxy))-above_pore_no_particle(d)) #analytic from top of membrane to base of particle
			Particle=(conductivity/math.pi)*quad(above_pore_with_particle,h-rProjetion(theta_d,rz,rxy),h+rProjetion(theta_d,rz,rxy),epsabs = 1e-9)[0] #numerical above pore particle present
			Above_particle=(conductivity/math.pi)*(above_pore_no_particle(lim)-above_pore_no_particle(h+rProjetion(theta_d,rz,rxy))) #analytic from top of particle to infinity 
			Below = (conductivity/math.pi)*(below_pore_no_particle(0)-below_pore_no_particle(-lim)) #analytic below pore no particle
			R = (Above_particle+Below_particle+Particle+Below+Within) 
			#print 'particle above pore'
					
		elif (h-rProjetion(theta_d,rz,rxy)) >= (d-2*rProjetion(theta_d,rz,rxy)): #transition between above to within
			Within_to_particle = quad(pore_no_particle,0,h-rProjetion(theta_d,rz,rxy),epsabs = 1e-9)[0] #numerical from base of membrane to bottom of particle
			Particle_to_top = quad(pore_with_particle,h-rProjetion(theta_d,rz,rxy),d,epsabs = 1e-9)[0] #numerical from base of particle to top of membrane
			Above_particle = quad(above_pore_with_particle,d,h+rProjetion(theta_d,rz,rxy),epsabs = 1e-9)[0] #numerical from top of membrane to top of particle
			Above = (above_pore_no_particle(lim)-above_pore_no_particle(h+rProjetion(theta_d,rz,rxy))) #analytic from top of particle to infinity
			Below = (below_pore_no_particle(0)-below_pore_no_particle(-lim)) #analytic from infinity to base of membrane without particle
	   		R = (conductivity/math.pi)*(real(Within_to_particle+Particle_to_top+Above_particle)+Above+Below)  
			#print 'transition between above to within'
			
		elif (h-rProjetion(theta_d,rz,rxy)) >= 0: #particle within pore
			Above_particle=(conductivity/math.pi)*quad(pore_no_particle,0,h-rProjetion(theta_d,rz,rxy),epsabs = 1e-9)[0] #numerical within pore no particle
			Particle=(conductivity/math.pi)*quad(pore_with_particle,h-rProjetion(theta_d,rz,rxy),h+rProjetion(theta_d,rz,rxy),epsabs = 1e-9)[0] #numerical within pore particle present (over area containing particle)
			Below_particle=(conductivity/math.pi)*quad(pore_no_particle,h+rProjetion(theta_d,rz,rxy),d,epsabs = 1e-9)[0] #numerical within pore no particle
			Above_and_below=(conductivity/math.pi)*((above_pore_no_particle(lim)-above_pore_no_particle(d))+(below_pore_no_particle(0)-below_pore_no_particle(-lim))) #analytic above and below no particle
			R = (Above_particle+Particle+Below_particle+Above_and_below)
			#print 'particle within pore'
			
		elif (h-rProjetion(theta_d,rz,rxy)) > (0-2*rProjetion(theta_d,rz,rxy)): #transition between below to within
			Within_to_particle = quad(pore_no_particle,0.,h+rProjetion(theta_d,rz,rxy))[0] #numerical from base to upper edge of particle
			Above_particle = quad(pore_no_particle,h+rProjetion(theta_d,rz,rxy),d)[0] #numerical from top of particle to top of pore
			Particle_to_bottom = quad(below_pore_with_particle,h-rProjetion(theta_d,rz,rxy),0)[0] #numerical from membrane base to lower edge of particle
			Above = above_pore_no_particle(lim)-above_pore_no_particle(d) #analytic above pore without particle
			Below = below_pore_no_particle(h-rz)-below_pore_no_particle(-lim) #  analytic from base of particle to infinity
			R = (conductivity/math.pi)*(real(Within_to_particle+Above_particle+Particle_to_bottom)+(Above)+(Below))
			#print 'transition between below to within'
			
		elif (h+rProjetion(theta_d,rz,rxy)) < 0: #particle below pore
			Below = (below_pore_no_particle(h-rProjetion(theta_d,rz,rxy))-below_pore_no_particle(-lim)) #analytic below particle to ininity
			Particle = quad(below_pore_with_particle,h-rProjetion(theta_d,rz,rxy),h+rProjetion(theta_d,rz,rxy))[0] #numerical below pore over diameter of particle
			Below_to_particle = (below_pore_no_particle(0)-below_pore_no_particle(h+rProjetion(theta_d,rz,rxy))) #analytic from membrane base to top of particle
			Within = quad(pore_no_particle,0,d)[0] #numerical solution within pore no particle
			Above = (above_pore_no_particle(lim)-above_pore_no_particle(d)) #analytic above pore no particle
			R = (conductivity/math.pi)*real(Below+Particle+Below_to_particle+Within+Above)
			#print 'particle below pore'
		else:
			print "Shit i've lost the particle!"
			sys.exit()
		I = V0/R
		dI = I-I0
		output_format = np.append(output_format,[[h,R,I, dI]],axis=0)
	
	#discards first row
	output_format = np.delete(output_format, 0, 0)
	
	R_list = output_format[:,1]
	z_list = output_format[:,0] #particle position with time
	I_list = output_format[:,2]
	Delta_I = output_format[:,3]
	
	Duration = calculate_FWHM(tspan, Delta_I)
	Magnitude = magn(Delta_I)
	
	return tspan, Delta_I, Icheck, Duration, J_at_d, Magnitude, z_list

#So for a non-rotated oblate spheroid
#the command would be
#pore_calculate(a,b,d,rxy,rz,P0,V0, psi0, theta_d, zeta_part)
#pore_calculate(6.0e-6,40e-6,200e-6,3.35e-6, 1.04e-6, 15*9.8,0.15, 0.013, 0.0, -0.02)
