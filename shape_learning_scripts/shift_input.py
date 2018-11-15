import numpy as np
import sys
import math
from joblib import Parallel, delayed
import multiprocessing
import sys
import os
sys.path.append(os.path.abspath("../common_scripts"))
import wig_and_signal_utils as wsu

"""
Input regions will be 1 1/2 times the expected size. Shift them to find the best
representation of the region. There will be a margin of 1/4 the region size on
either side. Keep the shifted version and print it out.

"Best Representation" means that it contains a large portion of its RPKM signal
near the center and that its edges contain little signal. We seek to maximize
a weighted sum of signals, where the weights drop off as the distance from the
center increases and become slightly negative at the edges.
"""
def main():

    file_path = sys.argv[1]
    output_path = sys.argv[2]
    bin_size = int(sys.argv[3])
    region_size = int(sys.argv[4])
    wig_file = open(sys.argv[5], 'r')
    percentile = 0.995
    fine = bool(sys.argv[6])
    minimum_intensity = float(sys.argv[7])
    
    threshold = wsu.get_intensity_percentile(percentile, wig_file, minimum_intensity, fine)
    print(str(threshold))
    
    shiftRegions(file_path, output_path, bin_size, region_size, threshold, minimum_intensity)
    
    #Print message to user.
    print("Shifting complete for all windows.")

#Find the best representation of each region by maximizing
#the sum of RPKM signals multiplied by the weight vector.
def shiftRegions(file_path, output_path, bin_size, region_size, threshold, minimum_intensity):

    #Print message to user.
    print("Shifting regions")

    #Open the file containing the regions and the output file.
    #Read the first line.
    input = open(file_path, 'r')
    output = open(output_path, 'w')
    next_line = input.readline()
    
    #Initialize the weight vector to use in finding the best representation.
    dim = int(region_size / bin_size)
    weightVector = np.zeros(dim)
    factor = 3 * dim / 8
    for i in range(0, dim):
        distance = min((abs(i - dim / 2)), abs(i - dim / 2 + 1))
        weightVector[i] = (factor - distance) / factor
    ctr = 0
    #For each region, print out its best representation.
    while next_line:
        #Obtain the labels and the signal intensities.
        split_line = next_line.split(",")
        labels = split_line[0:3]
        intensityStr = split_line[3:len(split_line)]
        intensities = [(float(i) - minimum_intensity) for i in intensityStr]
        old_labels = labels[:]
        
        #Find the best representation and print it out with the labels.
        finalIntensities, d = findBestRep(intensities, weightVector)
        finalIntensities = np.asarray(finalIntensities) + minimum_intensity
        
        #Get number of crossings across threshold.
        crossings = wsu.find_crossing_count(finalIntensities, threshold, minimum_intensity)
        
        #Update the labels.
        labels[1] = str(int(labels[1]) + d * bin_size)
        labels[2] = str(int(labels[2]) + d * bin_size)
        
        #Print shifted region.
        for label in labels:
            output.write(label)
            output.write(",")
        output.write(str(crossings))
        output.write(",")
        for val in range(0, dim):
            output.write(str(finalIntensities[val]))
            if val < dim - 1:
                output.write(",")
        output.write("\n")
        
        #Get the next line.
        next_line = input.readline()
        ctr += 1
    
#Find the best representation of each region by maximizing
#the sum of RPKM signals multiplied by the weight vector.
def findBestRep(sigs, w):

    #Store optimal delay and maximum so far.
    bestD = 0;
    maxSum = 0;
    
    #For each delay, compute its weighted sum and update max.
    for d in range(0, len(sigs) - len(w) + 1):
        sigVec = sigs[d:d + len(w)]
        sum = np.sum(sigVec * w)
        
        #If the weighted sum is the greatest so far, update.
        if sum > maxSum:
            bestD = d
            maxSum = sum

    #Return optimal vector.
    return sigs[bestD: bestD + len(w)], bestD

if __name__ == "__main__":
    main()