### Calculate Optical Path Delays, Projected Baselines, and Noise Models for all pairs of telescopes
### Made by Sahar Nikkhah
### Modified by Mike Lisa and Mackenzie Scott
### Modified by Josie Rose to add uv

### josie adapted 3 nov 2023 to run with input parameters from analysis macro

### Import Packages
from astropy.coordinates import EarthLocation
from astropy.time import Time, TimeDelta
from astropy.coordinates import AltAz
from astropy.coordinates import SkyCoord
from astropy.coordinates import ITRS
import numpy as np
#import matplotlib.pyplot as plt
import math
from datetime import datetime
import subprocess
import sys # allows to read in args when executing script

# read in parameters from command line - for run length, source, and local start time

longestRun = sys.argv[1]
source = sys.argv[2]
locStartDate = sys.argv[3]
locStartTime = sys.argv[4] # these are output as a single parameter from the root macro, but must be read in as 2 because they have a space

print("inputs", longestRun, source, locStartDate, locStartTime)

#Configure Source and Time  
TargetName = source[:3] + ' ' + source[3:]
start_observing_time = Time(locStartDate + ' ' + locStartTime)
start_observing_time += TimeDelta(25200,format='sec') # Convert to UTC time 
dt = TimeDelta(float(longestRun), format='sec')
end_observing_time = start_observing_time + dt 
Target = SkyCoord.from_name(TargetName)

print() 
print(TargetName, "target coords", Target) 
print("start time", start_observing_time, "end time", end_observing_time, "dt", dt)

### Parameters for noise analysis 
NoiseSourceLocation = [43.5, 55.5, 8.5]  # [0,0,0] would be the "middle" of the VERITAS site.
WavelengthOfNoise   = 1.75765            # in meters
CameraRadius        = 14                 # the distance of the camera from the pivot point of the telescope, in meters

### Location of VERITAS, telescopes, and Length of cables
Veritas      = EarthLocation(lat='31.675', lon='-110.952', height='1270')  #location of veritas
# telLocs      = np.array([[135.48, -8.61, 12.23], [44.1, -47.7, 4.4], [29.4, 60.1, 9.8], [-35.9, 11.3, 7.0]]) # old McGill/ASIIP coordinates
# telLocs      = np.array([[135.49, -8.72, 7.23], [45.145, -49.115, -0.94], [28.86, 60.73, 4.51], [-36.112, 11.696, 1.63]]) # these are the 2011 database coordinates
telLocs        = np.array([[135.48, -8.61, 12.23], [44.836, -49.601, 5.102], [29.335, 60.022, 10.636], [-35.885, 11.742, 6.417]]) # measured by Dave 2023, to cm precision
CableDelays  = [676.8e-9, 585.0e-9, 955.0e-9, 1063.7e-9] #Cable lengths for all 4 telescopes 
bucketLength = 4e-9

# There are six PAIRS of telescopes, and I index them as follows:
# index 0=(1,2); 1=(1,3); 2=(1,4); 3=(2,3); 4=(2,4); 5=(3,4);

# --> Important: the first telescope of a pair is considered to be the origin, so these are the angles relative to that
Dvec	= np.empty([6,3])  # The Distance between the telescopes as a 3D vector 
Dmag	= np.empty([6])    # The Magnitude of the Distance between the telescope pairs. 
Theta	= np.empty([6])    # Horizonial angle realtive to north
Beta	= np.empty([6])    # Angle realtive to the horizon 
Tpair   = np.empty([6,2])  # identifies which is the first and the second telescope


##Determine the distance and angle between all 6 pairs of telscopes 
iwhich = 0
for t1 in range(0,4):
	for t2 in range(t1+1,4):
		Tpair[iwhich,0] = t1
		Tpair[iwhich,1] = t2
		for icomp in range(0,3):
			Dvec[iwhich,icomp] = telLocs[t2,icomp] - telLocs[t1,icomp]
		Dmag[iwhich]	= math.sqrt(math.pow(Dvec[iwhich,0],2) + math.pow(Dvec[iwhich,1],2) + math.pow(Dvec[iwhich,2],2))
		Theta[iwhich]	= math.atan2(Dvec[iwhich,0],Dvec[iwhich,1])	  # note the ordering of indices is intentional.  Theta (and Azimuth) are measured relative to North
		Beta[iwhich]	= math.asin(Dvec[iwhich,2]/Dmag[iwhich])
		iwhich += 1
	
print("Dvec: " , Dvec) 

#Delays and baselines between the pairs
times = np.empty([0]) # time of day (in seconds)
delay	= [np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0])] # time delays between the six pairs of telescopes. (in seconds)
baseline = [np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0])] # next is projected baseline VECTOR (in meters)
baselinemag = [np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0])] #Magnitude of the baselines for each pair
altitude = np.empty([0]) # altitude angle in a numpy array
azimuth = np.empty([0]) # azimuth
local_ha = np.empty([0]) # Nolan's hour angle from ITRS coords - these are LOCAL ra/dec!!
local_dec = np.empty([0]) # Nolan's declination from ITRS

DistFromNoiseSource = np.empty([2])  # just used temporarily in loop
cameraPosition = np.empty([3])       # position of the CAMERA, which is at the end of a lever arm from the (fixed) telescope base position.  used temporarily in loop
DistDiff = [np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0])]
RelNoisePhase = [np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0]), np.empty([0])]


#---------------------------------------------------
# --------------------------------------------------
#           Iterate over the whole run 
#---------------------------------------------------
#---------------------------------------------------

observing_time = start_observing_time
while observing_time < end_observing_time:
    #Configure altitude and azimuth of star at given time 
	AltAndAz = AltAz(location=Veritas, obstime=observing_time) 
	ITRScoords = Target.transform_to(ITRS(obstime=observing_time))
	#print("ITRS coords: ", ITRScoords)
	alt=Target.transform_to(AltAndAz).alt.rad
	az=Target.transform_to(AltAndAz).az.rad
	altitude = np.append(altitude, alt)
	azimuth = np.append(azimuth, az)
	hrang = Veritas.lon.rad - ITRScoords.spherical.lon.rad
	if abs(hrang) > 2*math.pi:
		if hrang > 2*math.pi:
			hrang = hrang - (2*math.pi)
		else: 
			hrang = hrang + (2*math.pi)
	#print("ha:", hrang)
	decl = ITRScoords.spherical.lat.rad
        
	#print("dec: ", decl)
	local_ha = np.append(local_ha, hrang)
	local_dec = np.append(local_dec, decl)
	
	# add in an if statement to subtract 360 (abs val) if angle is larger than 360 for ha
	
	# no, this is COMPLETELY WRONG!!  see page 186 in logbook.  S = np.array([-math.cos(alt)*math.sin(az),-math.cos(alt)*math.sin(az),-math.sin(az)])	# light Poynting vector	
	S = np.array([-math.cos(alt)*math.sin(az),-math.cos(alt)*math.cos(az),-math.sin(alt)])	# light Poynting vector	
	for iwhich in range(0,6):
		delta		        = -Dmag[iwhich]*(math.cos(alt)*math.cos(Beta[iwhich])*math.cos(Theta[iwhich]-az) + math.sin(alt)*math.sin(Beta[iwhich]))
		delay[iwhich]	    = np.append(delay[iwhich],delta/3e8)
		baseline[iwhich]    = np.append(baseline[iwhich], [Dvec[iwhich,0]-delta*S[0], Dvec[iwhich,1]-delta*S[1], Dvec[iwhich,2]-delta*S[2]])
		baselinemag[iwhich] = np.append(baselinemag[iwhich],math.sqrt(math.pow(Dvec[iwhich,0]-delta*S[0],2)+math.pow(Dvec[iwhich,1]-delta*S[1],2)+math.pow(Dvec[iwhich,2]-delta*S[2],2)))
		# now some noise phase calculation...
		for TelInPair in range(0,2):
			cameraPosition[0] = telLocs[int(Tpair[iwhich,TelInPair]),0] + CameraRadius*math.cos(alt)*math.sin(az)
			cameraPosition[1] = telLocs[int(Tpair[iwhich,TelInPair]),1] + CameraRadius*math.cos(alt)*math.cos(az)
			cameraPosition[2] = telLocs[int(Tpair[iwhich,TelInPair]),2] + CameraRadius*math.sin(alt)
			SquaredDistFromSource = 0
			for icomp in range (0,3):
				SquaredDistFromSource += math.pow(cameraPosition[icomp]-NoiseSourceLocation[icomp],2)
			DistFromNoiseSource[TelInPair] = math.sqrt(SquaredDistFromSource)
		DistDiff[iwhich] = np.append(DistDiff[iwhich],DistFromNoiseSource[0]-DistFromNoiseSource[1])  	#DistDiff[iwhich] = np.append(DistDiff[iwhich],math.fmod(DistFromNoiseSource[0]-DistFromNoiseSource[1],WavelengthOfNoise))
		rnp = 2.0*math.pi*math.fmod(DistFromNoiseSource[0]-DistFromNoiseSource[1],WavelengthOfNoise)/WavelengthOfNoise
		rnp = rnp*(-1.0)       ###### THIS IS A TOTALLY UNJUSTIFIED SIGN CHANGE!  MACKENZIE AND I ARE WORKING ON THIS!!!
		if rnp < 0:
			rnp += 2.0*math.pi
		RelNoisePhase[iwhich] = np.append(RelNoisePhase[iwhich],rnp)   
	sectime = observing_time
	sectime.format = 'cxcsec'
	sectime.out_subfmt = 'decimal'
	times = np.append(times,sectime.value)
	observing_time += TimeDelta(10.0,format='sec')


### The pairs and time will be used in plot titles 
pairLabel = ['T1T2','T1T3','T1T4','T2T3','T2T4','T3T4']
ttt = Time(start_observing_time, format='iso', out_subfmt='date_hms')
timeString = ttt.value

# --------------------------------------------------
#             Baseline and OPD plots 
#---------------------------------------------------

#fig = plt.figure(figsize=(14, 8))
#fig.suptitle(TargetName+' '+timeString, fontsize=20)

  #Save optical path delays to a file 
for iwhich in range(0,6):
    outfile = open('delays/'+pairLabel[iwhich]+'.delay','w')
    for ipt in range(0,times.size):
        outfile.write('%d 	 %f\n' % (times[ipt], (1.0e9*delay[iwhich][ipt])))
    outfile.close()
	
	#Save projected baselines to a file 
    outfile = open('delays/'+pairLabel[iwhich]+'.baseline','w')
    for ipt in range(0,times.size):
        outfile.write('%d   %f\n' % (times[ipt],baselinemag[iwhich][ipt]))
    outfile.close()
    
    #Save projected baselines to a file 
    inc = 0
    outfile = open('delays/'+pairLabel[iwhich]+'.baselinecoord','w')
    for ipt in range(0,times.size):
        outfile.write('%d   %f    %f   %f\n' % (times[ipt], baseline[iwhich][inc], baseline[iwhich][inc+1], baseline[iwhich][inc+2]))
        inc += 3
    outfile.close()
	
	#Calculate Average Offset for correlations
    average = int(((CableDelays[int(Tpair[iwhich,1])] - CableDelays[int(Tpair[iwhich,0])]) - 0.5*(delay[iwhich][0]+delay[iwhich][times.size-1])) / bucketLength) #### HOW SHOULD WE KEEP DOING THIS IN THE FUTURE!!!!!!!!!
    outfile = open('delays/'+pairLabel[iwhich]+'aveNumberBucketsOffset','w')
    outfile.write('%d\n' % average)
    outfile.close()
    
outfile = open('delays/all.altaz','w')
outfile.write(str(times.size))
outfile.write("\n")
for ipt in range(0,times.size):
    outfile.write('%d   %f    %f   %f   %f\n' % (times[ipt], altitude[ipt], azimuth[ipt], local_ha[ipt], local_dec[ipt]))
outfile.close()

print("python all done!")
