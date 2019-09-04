from collections import namedtuple
import numpy as np
import wig_and_signal_utils as wsu
import math

"""
This class holds a region taken directly from the WIG file
with no alterations.
"""
class Region:
    # Status constants
    NONE = 0
    REGION_COMPLETED_WITH_ZEROS = 1
    REGION_COMPLETED_WITH_SIGNAL = 2
    REGION_INCOMPLETE = 3
    REGION_SIGNAL_OVERFLOW = 4
    
    def __init__(self, start, end, chrom, size):
        self.start = start
        self.end = end
        self.signals = self.make_default_signals(int(size))
        self.chromosome = chrom
        self.region_size = int(size)
        self.signal_pos = 0

    """
    Initialize signal list to contain all negatives,
    which will never be found in a real data set.
    """
    def make_default_signals(self, region_size):
        signals = np.ones(region_size) * (-1.0)
        return signals
        
    """
    Add the signal that was just read at the given position.
    If there was a skip but the signal just read was within
    the region boundaries, impute with the last signal.
    If there was a skip outside of the range of the signal,
    fill in zeros for the rest of the region.
    """
    def add_signal(self, raw_position, raw_position_end, signal, resolution):
        
        add_status = Region.NONE
        # If we have passed the end, fill in the last signal that we saw.
        if raw_position >= self.end:
            for i in range(self.signal_pos, self.region_size):
                self.signals[i] = 0.0
            add_status = Region.REGION_COMPLETED_WITH_ZEROS
                
        # Otherwise, fill in the signal.
        # If the end of the bin extends past the end of the region,
        # return an appropriate status.
        elif raw_position_end > self.end:
            position = int(math.floor((raw_position_end - self.start) / resolution))
            position_start = int(math.floor((raw_position - self.start) / resolution))
            while self.signal_pos < position_start:
                self.signals[self.signal_pos] = 0.0
                self.signal_pos = self.signal_pos + 1 
            while self.signal_pos < self.region_size:
                self.signals[self.signal_pos] = signal
                self.signal_pos = self.signal_pos + 1 
            add_status = Region.REGION_SIGNAL_OVERFLOW
        else:
            # Compute the position of the signal within the list.
            # Fill it in as needed.
            position = int(math.floor((raw_position_end - self.start) / resolution))
            position_start = int(math.floor((raw_position - self.start) / resolution))
            while self.signal_pos < position_start:
                self.signals[self.signal_pos] = 0.0
                self.signal_pos = self.signal_pos + 1 
            while self.signal_pos < position:
                self.signals[self.signal_pos] = signal
                self.signal_pos = self.signal_pos + 1
            
            if self.signal_pos == self.region_size:
                add_status = Region.REGION_COMPLETED_WITH_SIGNAL
            else:
                add_status = Region.REGION_INCOMPLETE
            
        return add_status
   
"""
This class holds a region that has been altered
according to optimization criteria.
"""   
class Shifted_Region:
 
    def __init__(self, region, region_size, factor, threshold):
        self.start = region.start
        self.end = region.end
        self.region_size = int(region_size)
        [self.signals, self.crossings] = self.shift_signal(region.signals, self.region_size, factor, threshold)
        self.chromosome = region.chromosome
    
    """
    Find the best representation of each region by maximizing
    the sum of RPKM signals multiplied by the weight vector.
    """
    def shift_signal(self, signals, dim, factor, threshold):
       
        # Initialize the weight vector to use in finding the best representation.
        weightVector = np.zeros(dim)
        for i in range(0, dim):
            distance = min((abs(i - dim / 2)), abs(i - dim / 2 + 1))
            factor_scaled_distance = factor * dim
            weightVector[i] = (factor_scaled_distance - distance) / factor_scaled_distance
        
        # Find the best representation.
        final_region = self.find_best_representation(signals, weightVector)
        final_region = np.asarray(final_region)

        # Get number of crossings across threshold.
        crossings = wsu.find_crossing_count(final_region, threshold)
        return [final_region, crossings]

    """
    Find the best representation of a region by maximizing
    the sum of RPKM signals multiplied by the weight vector.
    """
    def find_best_representation(self, sigs, w):

        #Store optimal delay and maximum so far.
        best_delay = 0;
        max_sum = 0;
        
        #For each delay, compute its weighted sum and update max.
        for d in range(0, len(sigs) - len(w) + 1):
            signal_subvector = sigs[d:d + len(w)]
            sum = np.sum(signal_subvector * w)
            
            #If the weighted sum is the greatest so far, update.
            if sum > max_sum:
                best_delay = d
                max_sum = sum

        #Return optimal vector.
        return sigs[best_delay: best_delay + len(w)]
        
"""
This class holds a learned shape.
"""   
class Shape:
    def __init__(self, name, mapped_count, signals):
        self.signals = signals
        self.name = name
        self.mapped_count = mapped_count
        
"""
This class holds a learned shape with its percentage distribution.
"""   
class Shape_Association:
    def __init__(shape, promoter, enhancer, repressor, weak):
        self.shape = shape
        self.promoter_percentage = promoter
        self.enhancer_percentage = enhancer
        self.repressor_percentage = repressor
        self.weak_percentage = weak
        
"""
This class holds a matched region.
"""   
class Matched_Region:
    def __init__(self, region, shape, crosscorr):
        self.start = region.start
        self.end = region.end
        self.signals = region.signals
        self.crossings = region.crossings
        self.chromosome = region.chromosome
        self.shape = shape
        self.crosscorr = crosscorr
        