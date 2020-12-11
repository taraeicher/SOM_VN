import sys
import numpy as np
import pandas as pd

"""
Find the range of for each type of shape.
"""
def main():

    #Read in the shape file.
    shapes = pd.read_csv(sys.argv[1], sep = "\t", header = None)
    enhancer_range = get_maximum_range("Enhancer", shapes)
    if enhancer_range[0] <= enhancer_range[1]:
        print("Enhancer Min: " + str(enhancer_range[0]))
        print("Enhancer Max: " + str(enhancer_range[1]))
    promoter_range = get_maximum_range("Promoter", shapes)
    if promoter_range[0] <= promoter_range[1]:
        print("Promoter Min: " + str(promoter_range[0]))
        print("Promoter Max: " + str(promoter_range[1]))
    weak_range = get_maximum_range("Weak", shapes)
    if weak_range[0] <= weak_range[1]:
        print("Weak Min: " + str(weak_range[0]))
        print("Weak Max: " + str(weak_range[1]))
    
"""
Return the range of maxima for an RE association.
"""
def get_maximum_range(re, shapes):

    #Exclude all shapes that do not contain the RE.
    shapes_re = shapes.loc[shapes[1] == re]
    
    # Find the range.
    minimum_max = 1000000
    maximum_max = 0
    for shape_idx in range(shapes_re.shape[0]):
        shape = shapes_re.iloc[shape_idx]
        current_max = get_maximum(shape)
        if current_max > maximum_max:
            maximum_max = current_max
        if current_max < minimum_max:
            minimum_max = current_max
            
    # Return the range.
    return([minimum_max, maximum_max])

        
"""
Compute the maximum for a shape.
"""
def get_maximum(shape):
    
    # Convert signals into numpy array.
    signals = shape[2]
    signal_array = np.array(signals.split(","), dtype = "float")
    
    # Return the maximum.
    maximum = np.max(signal_array)
    return(maximum)
    
    
if __name__ == "__main__":
    main()