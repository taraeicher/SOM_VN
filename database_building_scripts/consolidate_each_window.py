from joblib import Parallel, delayed
import multiprocessing
import sys
import glob, os
"""
Take a BED file with annotations and scores and consolidate the regions so that none overlap.
We want to maximize the total sum of region scores.
"""
def main():
    
    #Input, output, and window count.
    input_base = sys.argv[1]
    output_base = sys.argv[2]
    
    #List of lists for storing best set information. Each list corresponds to one window size.
    #List contains a vector denoting where to switch one annotation trajectory to another.
    do_switch = []
    
    #For each window size, find the optimal set of annotations.
    build_best_set(do_switch, input_base)
    write_best_set(do_switch, input_base, output_base)

"""
Given the scores in the file, find the optimal set of regions so that none overlap.
Here, optimal is defined as the set of regions with the highest sum of scores.
Because the overlap is half the region size, we have two "trajectories": the set of
odd scores and the set of even scores. The indices in switch_vec tell us when to
switch trajectories.

This is a dynamic programming solution that runs in linear time.
"""
def build_best_set(switch_vec, input_name):
    #Open the input file.
    input = open(input_name, "r")
    
    #For each region in the file, store the optimal set of regions until that region.
    #Save trajectory changes.
    sum_until = []
    [region_string, score] = get_next_region(input)
    i = 0
    while region_string:

        #Check whether this score + previous score in the trajectory is the best.
        #If not, switch. Update the highest sum in the vector.

        #Case 1: It's past at least the third score
        if len(sum_until) >= 2 :
            #Case 1a: The current trajectory is best.
            if score + sum_until[i - 2] >= sum_until[i - 1]:
                score += sum_until[i - 2]
                sum_until.append(score)
            #Case 1b: We should switch trajectories.
            else:
                sum_until.append(sum_until[i - 1])
                switch_vec.append(i)
        #Case 2: It's the second score.
        elif len(sum_until) == 1:
            #Case 2a: The current trajectory is best.
            if score > sum_until[0]:
                sum_until.append(score)
            #Case 2b: We should switch trajectories.
            else:
                sum_until.append(sum_until[0])
                switch_vec.append(1)
        #Case 3: It's the first score.
        else:
            sum_until.append(score)
            
        i += 1
        [region_string, score] = get_next_region(input)
        
    #Close the file.
    input.close()

"""
Return the next region in the file. This includes both the string and the score.
"""
def get_next_region(file):
    
    #Get the next line.
    next_line = file.readline()

    #Get the score.
    split_line = next_line.split()
    score = -1
    if len(split_line) == 6:
        score = float(split_line[4])

    #Return data.
    return [next_line, score]

"""
Given the list of switch locations, write the optimal set of regions.
""" 
def write_best_set(switch_vec, input_name, output_name):

    #Open the input file.
    input = open(input_name, "r")
    
    #Open the output file.
    output = open(output_name, "w")
    
    #Open the cluster files.
    #cluster_file = open(cluster_input_name, "r")
    #cluster_out = open(cluster_output_name, "w")
    
    #Find the switching position to determine which trajectory to start on.
    #Read lines until that point, then skip that line, then repeat.
    i = 0
    pos = -1
    [region_string, score] = get_next_region(input)

    for j in range(0, len(switch_vec)):
    
        #Print lines until reaching pos.
        old_pos = pos
        pos = switch_vec[j]

        while i < pos:
            
            #If pos is odd, write only even regions. If pos is even, write only odd regions.
            if (pos - i) % 2 == 1 and not i == old_pos:
                output.write(region_string)

            #Get next line.
            [region_string, score] = get_next_region(input)
            i += 1

    #Move to the end using the trajectory after the last switch.
    last_switch = -1
    if len(switch_vec) > 0:
        last_switch = switch_vec[len(switch_vec) - 1]

    while last_switch >= 0 and region_string:
        #Write in the same trajectory as the last switch.
        if i - last_switch % 2 == 0 and not i == last_switch:
            output.write(region_string)
        
        #Get next line.
        [region_string, score] = get_next_region(input)
        i +=1
    
    #Close files.
    input.close()
    #cluster_file.close()
    #cluster_out.close()
    output.close()
        
if __name__ == "__main__":
    main()