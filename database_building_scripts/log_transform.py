import sys
import math

"""
Transform the WIG input into log-scaled intensities.
"""
def main():

    #Input and output.
    input = open(sys.argv[1], 'r')
    output = open(sys.argv[2], 'w')
    limit = float(sys.argv[3])
    
    #Read each line.
    next_line = input.readline()
    while next_line:
        split_line = next_line.split()
        
        #If this line has data, convert it to log-scale.
        if len(split_line) == 2:
            val = float(split_line[1])
            try:
                log_val = math.log(val)
            except:
                log_val = math.log(limit)
            
            #Output the line.
            output.write(split_line[0] + "\t" + str(log_val) + "\n")
            
        else:
            output.write(next_line)
                
        #Read the next line.
        next_line = input.readline()
        
    #Close the input and output.
    input.close()
    output.close()
            
if __name__ == "__main__":
    main()