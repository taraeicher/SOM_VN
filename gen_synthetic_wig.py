"""
Generate wig the size of a single chromasome with the following characteristics:
1. A resolution of desired size (taken as input)
2. Regions with sharp peaks (simulated by Gaussian) with some randomness
3. Regions with broad peaks (simulated by Gaussian) with some randomness
4. Bimodal Gaussian peaks with some randomness
5. Skewed peaks
6. Large regions with multiple peaks close together
7. Highly erratic large regions
8. A background distribution near zero for 90% of the data
"""

from scipy import stats
from scipy import signal
import sys
import numpy as np
import math

#Parallel processing
from joblib import Parallel, delayed
import multiprocessing

def main():
	region_length = 640
	resolution_str = sys.argv[1] #resolution of file
	file_name = sys.argv[2] #name of output file
	snr_str = sys.argv[3] #signal-to-noise ratio
	signal_p_str = sys.argv[4] #percentage of regions containing only noise
	mean_noise_str = sys.argv[5] #mean noise level
	resolution = int(resolution_str)
	pos_count = (int)(62036316 / 24)
	snr = int(snr_str)
	signal_p = int(signal_p_str)
	mean_noise = int(mean_noise_str)
	names = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
	
	#Divide the signal into chromosomes for faster processing.
	Parallel(n_jobs=24)(delayed(generate_chromosome)(pos_count, region_length, names[i], resolution, snr, mean_noise, signal_p, file_name) for i in range(0, 24))

	
#Fill in the file data for a given chromosome.
#pos_count: number of positions to fill in
#chrom_name: name of chromosome to print to each line
#length: length of each region.
#resolution: resolution of the WIG (used in generating positions)
#snr: signal-to-noise ratio
#mean_noise: mean level of the noise region.
#signal_p: percentage of regions containing only noise
#file_name: base file name with path
def generate_chromosome(pos_count, length, chrom_name, resolution, snr, mean_noise, signal_p, file_name):
	#Open the file. Set the window count and figure out how much noise to add at the end.
	file = open(file_name + ".chr" + chrom_name + ".wig", 'w')
	window_count = int(pos_count / length)
	remaining_noise = pos_count % length;
	chrom_pos = 1

	#Add the first two lines.
	file.write("track type=wiggle_0 name=\"for chromosome " + chrom_name + "\" description=\"for chromosome " + chrom_name + "\"\n")
	file.write("variableStep chrom=" + chrom_name + " span=50\n")

	#Add data at the largest window resolution.
	for i in range(window_count):
		#Sample from a Gaussian. If the sample is within a z score of
		#1.645 (alpha = 0.05 for a two-tailed test), then impute one
		#of the specialized regions. Else, generate low signal data.
		decider = np.random.normal(0, 1, 1)
		sig_score = stats.norm.ppf(float(signal_p) / 100)
		if decider >= sig_score:
			#Randomly sample from a uniform distrbution and decide which
			#region to insert.
			region_decider = np.random.uniform(0, 1, 1)
			if region_decider < float(1) / 6:
				#Generate a sharp peak.
				chrom_pos = make_sharp_peak_region(length, chrom_pos, resolution, file, snr, mean_noise)
				
			elif region_decider < float(2) / 6:
				#Generate a broad peak.
				chrom_pos = make_broad_peak_region(length, chrom_pos, resolution, file, snr, mean_noise)
				
			elif region_decider < float(3) / 6:
				#Generate a left-skewed peak.
				chrom_pos = make_left_skewed_peak_region(length, chrom_pos, resolution, file, snr, mean_noise)
				
			elif region_decider < float(4) / 6:
				#Generate a right-skewed peak.
				chrom_pos = make_right_skewed_peak_region(length, chrom_pos, resolution, file, snr, mean_noise)
			
			elif region_decider < float(5) / 6:
				#Generate a bimodal peak.
				chrom_pos = make_bimodal_peak_region(length, chrom_pos, resolution, file, snr, mean_noise)
				
			elif region_decider <= float(6) / 6:
				#Generate a group of peaks.
				chrom_pos = make_grouped_region(length, chrom_pos, resolution, file, snr, mean_noise)
				
		else:
			chrom_pos = make_low_signal_region(length, chrom_pos, resolution, file, mean_noise)
	#Add noise at end to fill in the last part of the file.
	chrom_pos = make_low_signal_region(remaining_noise, chrom_pos, resolution, file, mean_noise)
		
#Create a region with a single sharp peak at a random location.
def make_sharp_peak_region(dim, chrom_pos, resolution, file, snr, mean_noise):

	#Make the region 1/32 of the region.
	width = int(dim / 32)
	sig = 2 * snr * mean_noise
	scl = mean_noise / 4
	height = abs(np.random.normal(loc = sig, scale = scl, size = 1))

	#Choose a random position to start the region.
	#Surround it by a low-signal area.
	position = int(np.random.uniform(low = 0, high = dim - width, size = 1))
	if(position > 0):
		chrom_pos = make_low_signal_region(position, chrom_pos, resolution, file, mean_noise)

	#Create the window and print it out.
	window = signal.gaussian(width, std = width / 6)
	for i in window:
		file.write(str(chrom_pos) + "\t" + str(float(i * height)) + "\n")
		chrom_pos += resolution

	#Insert low-signal area after window.
	if(dim - (position + width) > 0):
		chrom_pos = make_low_signal_region(dim - (position + width), chrom_pos, resolution, file, mean_noise)

	#Return chromosome position.
	return chrom_pos
		
#Create a region with a single broad peak at a random location.
def make_broad_peak_region(dim, chrom_pos, resolution, file, snr, mean_noise):

	#Make the region 1/4 of the region.
	width = int(dim / 4)
	sig = snr * mean_noise
	scl = mean_noise / 4
	height = abs(np.random.normal(loc = sig, scale = scl, size = 1))
	
	#Choose a random position to start the region.
	#Surround it by a low-signal area.
	position = int(np.random.uniform(low = 0, high = dim - width, size = 1))
	if(position > 0):
		chrom_pos = make_low_signal_region(position, chrom_pos, resolution, file, mean_noise)

	#Create the window and print it out.
	window = signal.gaussian(width, std = width / 6)
	for i in window:
		file.write(str(chrom_pos) + "\t" + str(float(i * height)) + "\n")
		chrom_pos += resolution

	#Insert low-signal area after window.
	if(dim - (position + width) > 0):
		chrom_pos = make_low_signal_region(dim - (position + width), chrom_pos, resolution, file, mean_noise)

	#Return chromosome position.
	return chrom_pos
		
#Create a region with a single skewed peak (poisson distribution)
def make_right_skewed_peak_region(dim, chrom_pos, resolution, file, snr, mean_noise):

	#Make the region 1/4 of the region.
	width = int(dim / 8)
	sig = 2 * snr * mean_noise
	scl = mean_noise / 4
	height = abs(np.random.normal(loc = sig, scale = scl, size = 1))
	
	#Choose a random position to start the region.
	#Surround it by a low-signal area.
	position = int(np.random.uniform(low = 0, high = dim - width, size = 1))
	if(position > 0):
		chrom_pos = make_low_signal_region(position, chrom_pos, resolution, file, mean_noise)

	#Create the window and print it out.
	for i in range(width):
		val = stats.gamma.pdf(i / 8, 2)
		file.write(str(chrom_pos) + "\t" + str(float(val * height)) + "\n")
		chrom_pos += resolution

	#Insert low-signal area after window.
	if(dim - (position + width) > 0):
		chrom_pos = make_low_signal_region(dim - (position + width), chrom_pos, resolution, file, mean_noise)

	#Return chromosome position.
	return chrom_pos
		
#Create a region with a single skewed peak (reverse poisson dist)
def make_left_skewed_peak_region(dim, chrom_pos, resolution, file, snr, mean_noise):

	#Make the region 1/4 of the region.
	width = int(dim / 8)
	sig = 2 * snr * mean_noise
	scl = mean_noise / 4
	height = abs(np.random.normal(loc = sig, scale = scl, size = 1))
	
	#Choose a random position to start the region.
	#Surround it by a low-signal area.
	position = int(np.random.uniform(low = 0, high = dim - width, size = 1))
	if(position > 0):
		chrom_pos = make_low_signal_region(position, chrom_pos, resolution, file, mean_noise)

		#Create the window and print it out.
	for i in range(width):
		val = stats.gamma.pdf((width - (i + 1)) / 8, 2)
		file.write(str(chrom_pos) + "\t" + str(float(val * height)) + "\n")
		chrom_pos += resolution

	#Insert low-signal area after window.
	if(dim - (position + width) > 0):
		chrom_pos = make_low_signal_region(dim - (position + width), chrom_pos, resolution, file, mean_noise)

	#Return chromosome position.
	return chrom_pos
		
#Create a region with a bimodal peak at a random location.
#Let 1/3 of the region be overlapped.
def make_bimodal_peak_region(dim, chrom_pos, resolution, file, snr, mean_noise):

	#Make the window 1/8 of the region.
	width_single = int(dim / 8)
	width_total = int((5 * width_single) / 3)
	width_actual = 0 #Needed to keep track of integer rounding
	sig = 1.5 * snr * mean_noise
	scl = mean_noise / 4
	height1 = abs(np.random.normal(loc = sig, scale = scl, size = 1))
	height2 = abs(np.random.normal(loc = sig, scale = scl, size = 1))
	
	#Choose a random position to start the region.
	#Surround it by a low-signal area.
	position = int(np.random.uniform(low = 0, high = dim - width_total, size = 1))
	if(position > 0):
		chrom_pos = make_low_signal_region(position, chrom_pos, resolution, file, mean_noise)

	#Create the window and print it out.
	#Print the non-intersecting part of the first gaussian
	window1 = signal.gaussian(width_single, std = width_single / 3)
	window2 = signal.gaussian(width_single, std = width_single / 3)
	for i in range(int((width_single * 2) / 3)):
		file.write(str(chrom_pos) + "\t" + str(float(window1[i] * height1)) + "\n")
		chrom_pos += resolution
	width_actual += int((width_single * 2) / 3)

	#Print the intersecting part of the peak.
	for i in range(int((width_single * 2) / 3) + 1, width_single):
		sig = max(window1[i] * height1, window2[i - int((width_single * 2) / 3) + 1] * height2)
		file.write(str(chrom_pos) + "\t" + str(float(sig)) + "\n")
		chrom_pos += resolution
	width_actual += width_single - int(((width_single * 2) / 3) + 1)
	
	#Print the non-intersecting part of the second gaussian
	for i in range(int(width_single / 3), int(width_single)):
		file.write(str(chrom_pos) + "\t" + str(float(window2[i] * height2)) + "\n")
		chrom_pos += resolution
	width_actual += width_single - int(width_single / 3)
	if(position + width_actual + dim - (position + width_actual) != 640):
		print(str(position) + " " + str(width_actual) + " " + str(dim - (position + width_actual)))
	#Insert low-signal area after window.
	if(dim - (position + width_actual) > 0):
		chrom_pos = make_low_signal_region(dim - (position + width_actual), chrom_pos, resolution, file, mean_noise)

	#Return chromosome position.
	return chrom_pos
		
#Create a region with a bimodal peak at a random location.
#Let 1/3 of the region be overlapped.
def make_grouped_region(dim, chrom_pos, resolution, file, snr, mean_noise):

	#Make the window 1/8 of the region.
	width = int(dim / 8)
	spacing1 = int(dim / 8)
	spacing2 = int(dim / 2)
	width_total = (4 * width) / 3
	sig1 = 1.5 * snr * mean_noise
	sig2 = snr * mean_noise
	sig3 = 2 * snr * mean_noise
	scl = mean_noise / 4
	height1 = abs(np.random.normal(loc = sig1, scale = scl, size = 1))
	height2 = abs(np.random.normal(loc = sig2, scale = scl, size = 1))
	height3 = abs(np.random.normal(loc = sig3, scale = scl, size = 1))
	
	#Start the region with the first peak.
	window1 = signal.gaussian(width, std = width / 6)
	for i in range(width):
		file.write(str(chrom_pos) + "\t" + str(float(window1[i] * height1)) + "\n")
		chrom_pos += resolution

	#Make the region between the first and second peaks.
	chrom_pos = make_low_signal_region(spacing1, chrom_pos, resolution, file, mean_noise)

	#Create the second peak.
	window2 = signal.gaussian(width, std = width / 6)
	for i in range(width):
		file.write(str(chrom_pos) + "\t" + str(float(window2[i] * height2)) + "\n")
		chrom_pos += resolution

	#Make the region between the first and second peaks.
	chrom_pos = make_low_signal_region(spacing2, chrom_pos, resolution, file, mean_noise)

	#Create the third peak.
	window3 = signal.gaussian(width, std = width / 6)
	for i in range(width):
		file.write(str(chrom_pos) + "\t" + str(float(window3[i] * height3)) + "\n")
		chrom_pos += resolution

	#Return chromosome position.
	return chrom_pos
		
#Create a region that is mostly low-signal and random.
#Choose points from a normal distribution.
def make_low_signal_region(dim, chrom_pos, resolution, file, mean_noise):

	scl = float(mean_noise) / 4
	for i in range(dim):
		sig = abs(np.random.normal(loc = mean_noise, scale = scl, size = 1))
		file.write(str(chrom_pos) + "\t" + str(float(sig)) + "\n")
		chrom_pos += resolution
		
	#Return chromosome position.
	return chrom_pos
		
if __name__ == "__main__":
	main()