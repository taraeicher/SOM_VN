#Need to run the SOM, then cluster on the resulting dim-reduced space.
#From https://codesachin.wordpress.com/2015/11/28/self-organizing-maps-with-googles-tensorflow/

#Import tensorflow and numpy for SOM.
import tensorflow as tf
import numpy as np
import copy
from scipy import stats
from scipy import spatial

#Import math functions
import math

#Import glob for obtaining the files.
import glob, os

#Import sys for obtaining command line args.
import sys
import common_ops as ops

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

    #Input, output, and windows.
    input = open(sys.argv[1], 'r')
    output_dir = sys.argv[2]
    wig_file = open(sys.argv[3], 'r')
    window_size = int(sys.argv[4])
    bin_size = int(sys.argv[5])
    intensity_threshold = ops.get_intensity_percentile(0.75, wig_file)
    print("INTENSITY THRESHOLD: " + str(intensity_threshold))
    
    #Create grids for labeling SOM nodes with cluster indices.
    som_centroids = []
    som_centroid_counts = []
    cluster_names = []
    cluster_count = 0

    #Create and train the SOMs and close the input files.
    train_som(input, som_centroids, som_centroid_counts, output_dir, intensity_threshold, window_size, bin_size)

    #Close files.
    input.close()
    wig_file.close()
    
"""
Method to create train SOM model for a given window size and print stats.
"""
def train_som(file, som_centroids, som_centroid_counts, out_dir, intensity_threshold, window, bin):
    
    #If file is not empty, learn SOM.
    if os.stat(file.name).st_size != 0:
        #Create new SOM
        window_size = window / bin
        som = SOM(window_size, file, intensity_threshold)
        
        #Train new SOM
        som.train(som.batch_size, file, window_size)
        
        #Print centroids
        som_centroids.append(som.get_centroids())
        som_centroid_counts.append(som.get_centroid_counts())
        print_centroids(som_centroids, som_centroid_counts, out_dir, file, window_size)

    #Print message to user.
    print("Grid for " + file.name + " is complete.")

"""
Print SOM node centroids and other metrics.
"""
def print_centroids(som_centroids, som_centroid_counts, out_dir, in_file_name, dim):

    #Print centroids.
    file = open(out_dir + "som_centroid", "w")
    for i in range(len(som_centroids)):
        for j in range(len(som_centroids[i])):
            for k in range(len(som_centroids[i][j])):
                file.write(','.join(map(str, som_centroids[i][j][k])))
                file.write("\n")
    file.close()
    
    #Print centroid counts.
    file = open(out_dir + "som_centroid_counts", "w")
    for i in range(len(som_centroid_counts)):
        for j in range(len(som_centroid_counts[i])):
            for k in range(len(som_centroid_counts[i][j])):
                file.write(str(som_centroid_counts[i][j][k]))
                file.write("\n")
    file.close()
    
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
    def __init__(self, dim, file, intensity_threshold):
        
        #Get metadata. The fractal dimension index measures jaggedness
        #of a signal. Dist and dist_q measure percentiles of the distance
        #from one region to the next random region in the dataset.
        #Line_count is the number of regions in the dataset.
        [fdi, dist, dist_q, line_count] = self.get_file_metadata(file, dim)
        
        #Set batch size to be at most half of the data. Ideally we want a large batch.
        batch_size = min(line_count / 2, 30000)
        self.batch_size = batch_size

        #Set grid size based on variability between regions. Note: Grid should not be larger than
        #number of lines in file. This will result in not enough training data.
        grid_size = min(max(0.20 * math.pow(dist, 2), 4), batch_size)
        m = int(math.sqrt(grid_size))
        n = int(math.sqrt(grid_size))
        self.m = m
        self.n = n
        
        #Set learning rate to be small. This prevents over-influence of low-signal regions.
        alpha = 0.2
        
        #Neighborhood size should also be based on the same metric. It is a fraction of the
        #length of one dimension of the grid.
        sigma = math.sqrt(min(max(0.20 * math.pow(dist_q, 2), 4), grid_size / 2))
        
        #Set number of iterations. We need more iterations if the regions are more jagged.
        self.iterations = int(fdi)      
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
            self.input_peak_count = tf.placeholder(tf.float32, shape = [batch_size, 1])
            self.label = tf.placeholder(tf.string, shape = [batch_size, 3])
            self.iteration_input = tf.placeholder("float")
            
            """
            Estimate the number of peaks in each node in the grid.
            """
            def get_peak_count(weights):
                    
                #Obtain a crude peak count by locating peaks and troughs based on cutoffs.
                peak_count = np.apply_along_axis(ops.find_peak_count, 1, weights, intensity_threshold)
                    
                #Return the estimated peak count.
                return peak_count.astype(float)
                
            #Calculate peak counts for all weightage vects.
            peak_counts_64 = tf.py_func(get_peak_count, [self.weightage_vects], tf.float64)
            self.peak_counts = tf.cast(peak_counts_64, tf.float32)
            
            #To compute the Best Matching Unit given a vector
            #Basically calculates the Euclidean distance between every
            #neuron's weightage vector and the input, and returns the
            #index of the neuron which gives the least value
            weightage_vects_stack = tf.stack([self.weightage_vects for i in range(batch_size)], axis = 2)
            peak_count_stack = tf.stack([self.peak_counts for i in range(batch_size)], axis = 1)
            batch_input_stack = tf.stack([tf.transpose(self.batch_input) for i in range(m*n)])
            input_peak_count_stack = tf.transpose(tf.contrib.layers.flatten(tf.stack([self.input_peak_count for i in range(m*n)], axis = 1)))
            diff = tf.subtract(weightage_vects_stack, batch_input_stack)
            total_max = tf.maximum(tf.reduce_max(weightage_vects_stack, axis = 1), tf.reduce_max(batch_input_stack, axis = 1))
            peak_count_max = tf.maximum(input_peak_count_stack, peak_count_stack)
            max_scales = tf.maximum(0.01, peak_count_max)
            bmu_indices = tf.argmin(tf.multiply(max_scales, tf.reduce_sum(tf.pow(diff, 2), 1)), axis = 0)
            
            #This will extract the location of the BMU based on the BMU's index)
            self.bmu_locations = tf.gather(self.location_vects, bmu_indices)
            self.bmu_weights = tf.gather(self.weightage_vects, bmu_indices)
            self.peak_count_weights = tf.gather(self.peak_counts, bmu_indices)
            self.lambdas = tf.maximum(0.01, self.peak_count_weights)
            
            #To compute the alpha and sigma values based on iteration number
            learning_rate_op = tf.subtract(1.0, tf.divide(self.iteration_input, self.iterations))
            self.alpha_op = tf.multiply(alpha, learning_rate_op)
            sigma_op = tf.multiply(sigma, learning_rate_op)

            #Construct the op that will generate a vector with learning
            #rates for all neurons, based on iteration number and location wrt BMU.
            location_vect_stack = tf.stack([self.location_vects for i in range(batch_size)], axis = 1)
            bmu_location_stack = tf.stack([self.bmu_locations for i in range(m*n)])
            bmu_distance_squares = tf.reduce_sum(tf.pow(tf.subtract(location_vect_stack, bmu_location_stack), 2), 2)
            self.neighborhood_funcs = tf.exp(tf.negative(tf.divide(tf.cast(bmu_distance_squares, "float32"), tf.multiply(tf.pow(sigma_op, 2), self.lambdas))))
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
    def train(self, batch_size, file, dim):
        #Fill in local variable to hold placeholder values.
        label_batch = []
        batch_data = []
        first_lines = []
        peak_count = []
        iters = 0
        
        #Track count of regions mapping to each centroid.
        map_count = np.zeros((self.m,self.n))
        
        #Train in mini-batches
        with self.sess:
            #Training iterations
            while iters < self.iterations:
                try:
                    #Fill in input data and train.
                    loop_count = 0
                    while True:
                        #Get data for training.
                        self.fill_in_data(file, label_batch, batch_data, first_lines, loop_count, dim, batch_size, peak_count)
                        
                        #Update location counts if we are on the last iteration.
                        if iters == self.iterations - 1:
                            locations = self.bmu_locations.eval(feed_dict = {self.batch_input: batch_data, self.label: label_batch, self.iteration_input: float(iters), self.input_peak_count: peak_count})
                            for loc in range(batch_size):
                                map_count[locations[loc][0], locations[loc][1]] += 1

                        #Train the model.
                        #if dim == 640:
                        #   deltas = self.add_delta.eval(session = self.sess, feed_dict = {self.batch_input: batch_data, self.label: label_batch, self.iteration_input: float(iters)})
                        #   print(str(deltas))
                        self.peak = list(self.peak_counts.eval(feed_dict = {self.batch_input: batch_data, self.label: label_batch, self.iteration_input: float(iters), self.input_peak_count: peak_count}))
                        if dim == 640:
                            print(str(np.mean(np.amax(self.peak)))  + " " + str(np.mean(self.peak)))
                        self.sess.run(self.shrink_weightage, feed_dict = {self.batch_input: batch_data, self.label: label_batch, self.iteration_input: float(iters), self.input_peak_count: peak_count})
                        self.sess.run(self.add_delta, feed_dict = {self.batch_input: batch_data, self.label: label_batch, self.iteration_input: float(iters), self.input_peak_count: peak_count})
                        
                        #For every 1,000 samples, print a dot.
                        if loop_count % 1000 == 0:
                            sys.stdout.flush()
                            sys.stdout.write('.')
                        loop_count += 1
                        
                        #Re-set lists.
                        label_batch = []
                        batch_data = []
                        peak_count = []
                            
                        #Update weightages.
                        self.weightages = list(self.weightage_vects.value().eval())

                #When the end of the file has been reached, train once more with last batch if necessary.
                except EOFError:
                    pass    
        
                #Re-set lists and start at beginning of file again.
                label_batch = []
                batch_data = []
                first_lines = []
                file.seek(0)
                iters += 1
                
            #Store a centroid grid and list of counts for easy retrieval later on
            centroid_grid = [[] for i in range(self.m)]
            centroid_counts = [[] for i in range(self.m)]
            self.locations = list(self.sess.run(self.location_vects))
            for i, loc in enumerate(self.locations):
                centroid_grid[loc[0]].append(self.weightages[i])
                centroid_counts[loc[0]].append(int(map_count[loc[0],loc[1]]))
            self.centroid_grid = centroid_grid
            self.centroid_counts = centroid_counts
            
            self.trained = True
    """
    Fills in the label and data values needed for training.
    """
    def fill_in_data(self, file, labels, inputs, first_lines, loop_count, dim, batch_size, peak_count):
        #Data strings
        inputStr = []
        
        #Get enough lines to fill in mini-batch for training.
        column = 0
        while column in range(batch_size):
            next_line = file.readline()
            #If end of file is reached, throw an error.
            if not next_line:
                #Fill in rest of batch with beginning of file.
                if column > 0:
                    for k in range(column, batch_size, 1):
                        labels.append(first_lines[k][0:3])
                        peak_count.append([split_line[3]])
                        inputStr.append(first_lines[k][4:dim + 4])
                        column += 1
                else:
                    raise EOFError()
            #If not at end of file, fill in the values.
            else:
                split_line = next_line.split(",")
                labels.append(list(split_line[0:3]))
                peak_count.append([split_line[3]])
                inputStr.append(list(split_line[4:dim + 4]))
                
                #Initialize beginning of file to use when end is reached and mini-batch is not filled in completely.
                if loop_count == 0:
                    first_lines.append(split_line)
                    
            column += 1 
        #Convert input strings to floats.
        for j in range(len(inputStr)):
            inputs.append([float(i) for i in inputStr[j]])
    
    """
    Returns the Xth percentile of jaggedness for all records in the file,
    as measured by the Fractal Dimension Index (FDI). 
    """
    def get_file_metadata(self, file, dim):
        
        #Count of the number of inputs with each FDI and max.
        #Note: for window size of 10, fdi_threshold is about 6.
        #FDI scaling factor allows us to store the FDI data at
        #a finer granularity.
        max_threshold = 1000
        jag_threshold = int(math.ceil(math.sqrt(math.pow((2 * max_threshold), 2)) * dim))
        jag_counts = np.zeros(jag_threshold)
        dist_counts = np.zeros(jag_threshold)
        dist_sum = 0
        dist_count = 0
        
        file_line_count = 0
        
        #Use each entry in the file to calculate running metadata.
        next_line = file.readline()
        split_line_p = []
        val_str_p = []
        val_p = []
        prev_line = False
        
        while next_line:
            split_line = next_line.split(",")
            val_str = split_line[3:dim + 3]
            val = [float(i) for i in val_str]
            
            #Increment the count of lines in the file.
            file_line_count += 1
            
            #Calculate the jaggedness metric.
            length = 0
            prev_val = -1
            for v in val:
                if prev_val >= 0:
                    length += math.sqrt(math.pow(v - prev_val, 2))
                prev_val = v
            
            #Increment the appropriate location in the FDI array.
            jag_counts[int(math.ceil(length))] += 1
            
            #Calculate the distance from this record to the previous one.
            if prev_line:
                dist = np.linalg.norm(np.asarray(val) - np.asarray(val_p))
                dist_sum += dist
                dist_count += 1
                dist_counts[int(math.ceil(dist))] += 1
                
            #Read the next line.
            prev_line = next_line
            split_line_p = prev_line.split(",")
            val_str_p = split_line_p[3:dim + 3]
            val_p = [float(i) for i in val_str_p]
            next_line = file.readline()
            
        running_sum = 0
        i = 0
        percentile_found = False
        fdi_percentile = 0
        dist_percentile = 0
        dist_percentile_1 = 0       
        
        #Find percentile of FDIs.
        target_count = int(file_line_count * 0.95)
        while i < range(len(jag_counts) - 1) and not percentile_found:
            running_sum += jag_counts[i]
            if running_sum >= target_count:
                fdi_percentile = i
                percentile_found = True
            i += 1
        i = 0
        percentile_found = False
        running_sum = 0
            
        #Find percentile of distances.
        while i < range(len(dist_counts) - 1) and not percentile_found:
            running_sum += dist_counts[i]
            if running_sum >= target_count:
                dist_percentile = i
                percentile_found = True
            i += 1
        i = 0
        percentile_found = False
        running_sum = 0
        
        #Find percentile of distances.
        tar = int(file_line_count * 0.25)
        while i < range(len(dist_counts) - 1) and not percentile_found:
            running_sum += dist_counts[i]
            if running_sum >= tar:
                dist_percentile_1 = i
                percentile_found = True
            i += 1
        i = 0
        percentile_found = False
        running_sum = 0
            
        #Return values.
        return [float(fdi_percentile), float(dist_percentile), float(dist_percentile_1), file_line_count]
                
    """
    Returns a list of 'm' lists, with each inner list containing
    the 'n' corresponding centroid locations as 1-D NumPy arrays.
    """
    def get_centroids(self):
        if not self.trained:
            raise ValueError("SOM not trained yet")
        return self.centroid_grid
        
    """
    Returns a list of 'm' lists, with each inner list containing
    the 'n' corresponding centroid locations as 1-D NumPy arrays.
    """
    def get_centroid_counts(self):
        if not self.trained:
            raise ValueError("SOM not trained yet")
        return self.centroid_counts

if __name__ == "__main__":
    main()