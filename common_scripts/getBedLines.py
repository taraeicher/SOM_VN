import sys
import numpy as np

#main function
def main():

    #Get arguments.
    wig_file = open(sys.argv[1], "r")
    bed = open(sys.argv[2], "r")
    out_file = sys.argv[3]
    chrom = sys.argv[4]

    printIntervals(wig_file, bed, out_file)
    
    #Close WIG file and output file.
    wig_file.close()
    bed.close()
    
    #Notify user that process has completed.
    print("Files complete for chromosome " + str(chrom));


#Build the list of peaks from all peak callers and all chromosomes.
def printIntervals(wig_file, bed, out_file):
        
    #Variables that will be used when traversing file.
    interval_start = 0;
    interval_end = 0;
    prev_end = 0;
    out = open(out_file, "w")
    
    #We don't care about the wig header.
    junk = wig_file.readline()
    junk = wig_file.readline()
    
    #Extract the first signal.
    peak = bed.readline().split("\t")
    next_line = wig_file.readline().split("\t")
    
    #Go through until reaching the end of the file.
    #Fill in the counts for each signal level.
    bin_size = 50
    while len(next_line) > 1 and len(peak) > 1:
    
        #Get the interval and signal.
        prev_end = interval_end
        interval_start = int(peak[1])
        interval_end = int(peak[2])
        chrom_pos = int(next_line[0])

        #We only want to count regions that do not overlap.
        peak_sigs = np.zeros(int((interval_end - interval_start) / bin_size))
        if interval_start >= prev_end:

            #Search for the marker that is equivalent to the first number in the interval.
            while chrom_pos < interval_start and len(next_line) > 0:
                #Get the next line.
                next_line = wig_file.readline().split("\t")
                chrom_pos = int(next_line[0])
            
            #Fill in the signal.
            i = 0  
            init_chrom_pos = 0
            while chrom_pos < interval_end and len(next_line) > 0:
                #If i == 0, this is the first position for this peak.
                if i == 0:
                    #Extract the signal and position.
                    peak_sigs[i] = float(next_line[1])
                    init_chrom_pos = chrom_pos
                
                #Check that the position in the chromosome matches where it should be.
                elif chrom_pos == init_chrom_pos + i * bin_size:
                    #Extract the signal and position.
                    peak_sigs[i] = float(next_line[1])
               
                #Get the next line.
                next_line = wig_file.readline().split("\t")
                chrom_pos = int(next_line[0])
                i += 1
                    
            #Print out the signals.
            out.write(",".join([str(s) for s in peak_sigs]) + "\n")

        #Get next peak entry.
        peak = bed.readline().split("\t")
            
    out.close()

if __name__ == "__main__":
    main()  