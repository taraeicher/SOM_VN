"""
This is a tensorflow implementation of a self-organizing map with variable neighborhood size.
The neighborhood size is based on the number of crossings of a discrete signal above and
below an input threshold. This code is based on https://codesachin.wordpress.com/2015/11/28/self-organizing-maps-with-googles-tensorflow/. It
requires the following parameters:
1. A pickled set of regions, where regions are defined in region_defs.py.
2. An output directory where the shapes learned by the SOM should be stored.
3. A file containing the intensity threshold for this chromosome.
4. The size of the input region.
5. The size of the bins in the WIG file.
"""

#Import all needed packages.
import tensorflow as tf
import numpy as np
from tqdm import tqdm
import math
import os
import sys
sys.path.append(os.path.abspath("../common_scripts"))
import region_defs
import wig_and_signal_utils as wsu
import pickle as pkl
import random

"""
Main function accomplishes the following tasks:
1. Train an SOM for each window size and print out the
evaluation metrics and node centroids.
2. Cluster the SOM nodes for each window size and print
out the cluster evaluation metrics and cluster centroids.
3. Overlay the clusters with relevant annotations.
4. Output BED files containing the annotations for each region
in the original training data.
"""
def main():

    #Input, output, and regions.
    regions = pkl.load(open(sys.argv[1], 'rb'))
    output = sys.argv[2]
    region_size = int(sys.argv[3])
    alpha = float(sys.argv[4])
    sigma = float(sys.argv[5])
    grid_size = float(sys.argv[6])
    iterations = float(sys.argv[7])
    bin_size = float(sys.argv[8])
    region_bins = region_size / bin_size
    
    #Create list of shapes from the grid.
    som_shapes = []

    #Create and train the SOMs and close the input files.
    train_som(regions, som_shapes, output, region_size, alpha, sigma, grid_size, iterations)
    
    #Print message to user.
    print("Grid for " + sys.argv[1] + " is complete.")
    
"""
Method to create train SOM model for a given window size and print stats.
"""
def train_som(regions, som_shapes, out, region_size, alpha, sigma, grid_size, iterations):
    
    #If there is at least one region, learn SOM.
    if len(regions) > 0:
        #Create new SOM
        som = SOM(regions, region_size, alpha, sigma, grid_size, iterations)
        
        #Train new SOM
        som.train(som.batch_size, regions, region_size)
        
        #Save centroids
        som_shapes.append(som.get_shapes())
        pkl.dump(som_shapes, open(out, "wb"))
    
""" 
2-D Self-Organizing Map with Gaussian Neighbourhood function
and linearly decreasing learning rate.
"""
class SOM(object):

    #To check if the SOM has been trained
    trained = False
    
    """
    Initializes all necessary components of the TensorFlow
    Graph.

    m X n are the dimensions of the SOM. 'n_iterations' should
    should be an integer denoting the number of iterations undergone
    while training.
    'dim' is the dimensionality of the training inputs.
    'alpha' is a number denoting the initial time(iteration no)-based
    learning rate. Default value is 0.3
    'sigma' is the the initial neighbourhood value, denoting
    the radius of influence of the BMU while training. By default, its
    taken to be half of max(m, n).
    """
    def __init__(self, regions, dim, alpha, sigma, grid_size, iterations):
        
        #Set parameters.
        batch_size = int(math.floor(len(regions) / 10))
        self.batch_size = batch_size
        m = int(math.sqrt(grid_size))
        n = int(math.sqrt(grid_size))
        self.m = m
        self.n = n
        self.iterations = iterations
        self.weightages = []
        
        ##INITIALIZE GRAPH
        self.graph = tf.Graph()
        
        ##POPULATE GRAPH WITH NECESSARY COMPONENTS
        with self.graph.as_default():
        
            ##VARIABLES AND CONSTANT OPS FOR DATA STORAGE
            #Randomly initialized weightage vectors for all neurons,
            #stored together as a matrix Variable of size [m*n, dim]
            self.weightage_vects = tf.Variable(tf.random_uniform(shape = [m*n, dim], minval = 1, maxval = 100, dtype = tf.float32)) 

            #Matrix of size [m*n, 2] for SOM grid locations
            #of neurons
            self.location_vects = tf.constant(np.array(list(self.neuron_locations(m, n))))
            
            ##PLACEHOLDERS FOR TRAINING INPUTS
            self.batch_input = tf.placeholder(tf.float32, shape = [batch_size, dim])
            self.label = tf.placeholder(tf.string, shape = [batch_size, 3])
            self.iteration_input = tf.placeholder("float")

            """
            Get the cross-correlation distance between each region-weight pair.
            """
            def get_vectorized_crosscorr(weights_raw, inputs_raw):
            
                # Clip all values to have a very small (but non-zero) value.
                very_small = 0.001
                signal_axis = 1
                weights = np.clip(weights_raw, a_min = very_small, a_max = None)
                inputs = np.clip(inputs_raw, a_min = very_small, a_max = None)
                sum_weights = np.sum(weights, axis = signal_axis)
                sum_inputs = np.sum(inputs, axis = signal_axis)
                length = weights.shape[signal_axis]
                
                # Compute normalized input-to-weights cross-correlation value.
                R_weights_inputs = np.sum(weights * inputs, axis = signal_axis)
                avg_weight_by_input = (sum_weights * sum_inputs) / length
                R_weights_inputs_norm = R_weights_inputs - avg_weight_by_input
                
                # Compute normalized input-to-input autocorrelation value.
                R_inputs_inputs = np.sum(inputs * inputs, axis = signal_axis)
                avg_input_by_input = (sum_inputs * sum_inputs) / length
                R_inputs_inputs_norm = np.clip(R_inputs_inputs - avg_input_by_input, a_min = 0, a_max = None)
                
                # Compute normalized weight-to-weight autocorrelation value.
                R_weights_weights = np.sum(weights * weights, axis = signal_axis)
                avg_weight_by_weight = (sum_weights * sum_weights) / length
                R_weights_weights_norm = np.clip(R_weights_weights - avg_weight_by_weight, a_min = 0, a_max = None)
                
                # Compute final metric
                crosscorr = R_weights_inputs_norm / np.maximum(R_inputs_inputs_norm, R_weights_weights_norm)
                return crosscorr
            
            #To compute the Best Matching Unit given a vector
            #Basically calculates the Euclidean distance between every
            #neuron's weightage vector and the input, and returns the
            #index of the neuron which gives the least value
            self.weightage_vects_stack = tf.stack([self.weightage_vects for i in range(batch_size)], axis = 2)
            self.batch_input_stack = tf.stack([tf.transpose(self.batch_input) for i in range(m*n)])
            self.crosscorrs = tf.py_func(get_vectorized_crosscorr, [self.weightage_vects_stack, self.batch_input_stack], tf.float32)
            self.bmu_indices = tf.argmax(self.crosscorrs)
            
            #This will extract the location of the BMU based on the BMU's index)
            self.bmu_locations = tf.gather(self.location_vects, self.bmu_indices)
            self.bmu_weights = tf.gather(self.weightage_vects, self.bmu_indices)
            
            #To compute the alpha and sigma values based on iteration number
            learning_rate_op = tf.subtract(1.0, tf.divide(self.iteration_input, self.iterations))
            self.alpha_op = tf.multiply(alpha, learning_rate_op)
            sigma_op = tf.multiply(sigma, learning_rate_op)

            #Construct the op that will generate a vector with learning
            #rates for all neurons, based on iteration number and location wrt BMU.
            location_vect_stack = tf.stack([self.location_vects for i in range(batch_size)], axis = 1)
            bmu_location_stack = tf.stack([self.bmu_locations for i in range(m*n)])
            bmu_distance_squares = tf.reduce_sum(tf.pow(tf.subtract(location_vect_stack, bmu_location_stack), 2), 2)
            self.neighborhood_funcs = tf.exp(tf.negative(tf.divide(tf.cast(bmu_distance_squares, "float32"), tf.pow(sigma_op, 2))))
            self.learning_rate_ops = tf.multiply(self.alpha_op, self.neighborhood_funcs)

            #Finally, the op that will use learning_rate_op to update
            #the weightage vectors of all neurons based on a particular input
            self.neighborhood_multipliers = tf.stack([tf.tile(tf.slice(self.neighborhood_funcs, np.array([i, 0]), np.array([1, batch_size])), [dim, 1]) for i in range(m*n)])
            learning_rate_multipliers = tf.stack([tf.tile(tf.slice(self.learning_rate_ops, np.array([i, 0]), np.array([1, batch_size])), [dim, 1]) for i in range(m*n)])
            self.weightage_delta_numerator = tf.reduce_sum(tf.multiply(learning_rate_multipliers, tf.stack([tf.transpose(self.batch_input) for i in range(m*n)])), 2)
            self.denominator_mins = tf.convert_to_tensor(0.00001 * np.ones((m*n, dim)), dtype = "float32")
            self.weightage_delta_denominator = tf.maximum(tf.reduce_sum(self.neighborhood_multipliers, 2), self.denominator_mins)
            self.weightage_deltas = tf.divide(self.weightage_delta_numerator, self.weightage_delta_denominator)
            self.weightage_vects_alpha = tf.multiply(self.weightage_vects, tf.subtract(1.0, self.alpha_op))
            self.shrink_weightage = tf.assign(self.weightage_vects, self.weightage_vects_alpha, use_locking = True)
            self.add_delta = tf.assign_add(self.weightage_vects, self.weightage_deltas, use_locking = True)
            
            ##INITIALIZE SESSION
            self.sess = tf.Session()
            
            ##INITIALIZE VARIABLES
    
            init_op = tf.global_variables_initializer()
            self.sess.run(init_op)

    """
    Yields one by one the 2-D locations of the individual neurons
    in the SOM.
    """
    def neuron_locations(self, m, n):
        #Nested iterations over both dimensions
        #to generate all 2-D locations in the map
        for i in range(m):
            for j in range(n):
                yield np.array([i, j])
    
    """
    Trains the SOM.
    'input_vects' should be an iterable of 1-D NumPy arrays with
    dimensionality as provided during initialization of this SOM.
    Current weightage vectors for all neurons(initially random) are
    taken as starting conditions for training.
    """
    def train(self, batch_size, regions, dim):
        #Fill in local variable to hold placeholder values.
        self.old_weightages = []
        label_batch = []
        batch_data = []
        first_regions = []
        num_weight_crossings = []
        
        #Track count of regions mapping to each shape.
        map_count = np.zeros((self.m,self.n))
        
        #Train in mini-batches
        with self.sess:
            #Training iterations
            for epoch in range(self.iterations):
                print("Training epoch " + str(epoch))
                #Permute the data for this epoch.
                random.shuffle(regions)
                #Fill in input data and train.
                for minibatch in tqdm(range(math.floor(len(regions) / batch_size))):
                    #Get data for training.
                    [label_batch, batch_data] = self.fill_in_data(regions, minibatch, batch_size)
                    #Update location counts if we are on the last iteration.
                    if epoch == self.iterations - 1:
                        locations = self.bmu_locations.eval(feed_dict = {self.batch_input: batch_data, self.label: label_batch, self.iteration_input: float(epoch)})
                        for loc in range(batch_size):
                            map_count[locations[loc][0], locations[loc][1]] += 1
                  
                    #Train the model.
                    self.sess.run(self.shrink_weightage, feed_dict = {self.batch_input: batch_data, self.label: label_batch, self.iteration_input: float(epoch)})
                    self.sess.run(self.add_delta, feed_dict = {self.batch_input: batch_data, self.label: label_batch, self.iteration_input: float(epoch)})
                    
                    cc = self.crosscorrs.eval(feed_dict = {self.batch_input: batch_data, self.label: label_batch, self.iteration_input: float(epoch)})
                        
                    #Update weightages.
                    self.old_weightages = self.weightages
                    self.weightages = list(self.weightage_vects.value().eval()) 
                    
                # Print the change in weights.
                if epoch > 0:
                    print("Delta: " + str(np.sum(np.absolute(np.asarray(self.weightages) - np.asarray(self.old_weightages)))))
                
            #Store a centroid grid and list of counts for easy retrieval later on
            shapes = []
            self.locations = list(self.sess.run(self.location_vects))
            for i, loc in enumerate(self.locations):
                shapes.append(region_defs.Shape(str(loc[0]) + "_" + str(loc[1]), map_count[loc[0],loc[1]], self.weightages[i]))
            self.shapes = shapes
            
            self.trained = True
    """
    Fills in the label and data values needed for training.
    """
    def fill_in_data(self, regions, minibatch_idx, batch_size):
        
        #Get enough lines to fill in mini-batch for training.
        index = 0
        label_batch = []
        batch_data = []
        num_weight_crossings = []
        while index in range(batch_size):
            
            # Fill in the next region.
            region = regions[minibatch_idx * batch_size + index]
            label_batch.append([str(region.chromosome), str(region.start), str(region.end)])
            batch_data.append([region.signals])
                    
            index += 1 
            
        label_batch_arr = np.stack(label_batch)
        batch_data_arr = np.squeeze(np.stack(batch_data))
            
        return [label_batch_arr, batch_data_arr]
                
    """
    Returns a list of shapes, each one containing location in the grid, counts, and weightages.
    """
    def get_shapes(self):
        if not self.trained:
            raise ValueError("SOM not trained yet")
        return self.shapes

if __name__ == "__main__":
    main()