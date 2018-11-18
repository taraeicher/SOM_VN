#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
#include "get_bed_lines.h"

int counter;
int gapFractionAllowed = 0;

//main function
int main(int argc, char* argv[])
{
	//Initialize random number generator.
	srand(5000);
	
	//Get arguments.
	char* wigFileName = argv[1];
	char* training = argv[2];
    char* chrom = argv[3];
	char* outFileName = argv[4];
    char* outBedName = argv[5];
	
	//Open WIG file.
	FILE* wigFile = fopen(wigFileName, "r");
    FILE* outFile = fopen(outFileName, "w");
    FILE* outBed = fopen(outBedName, "w");
    FILE* inBed = fopen(training, "r");
	
	//Output intervals to files for processing.
	PrintIntervalArgs printIntervalArgs;
	void* args3;
	printIntervalArgs.wigFile = wigFile;
	printIntervalArgs.bedInput = inBed;
	printIntervalArgs.outFile = outFile;
    printIntervalArgs.bedFile = outBed;
	printIntervalArgs.chrom = chrom;
	args3 = &printIntervalArgs;
	printIntervals(args3);
	
	//Close WIG file and output file.
	fclose(wigFile);
    fclose(outFile);
    fclose(outBed);
    fclose(inBed);
	
	//Notify user that process has completed.
	printf("Files complete for chromosome %s\n",chrom);
	return(0);
}

//Build the list of peaks from all peak callers and all chromosomes.
void* printIntervals(void* argBundle)
{	
	//Dereference arguments.
	PrintIntervalArgs* args = (PrintIntervalArgs*)argBundle;
	
	//Get the chromosome from the chromosome index.
	char* chrom = args->chrom;
	
	//Check for empty file.
	char firstChar;
	firstChar = fgetc(args->wigFile);
	ungetc(firstChar, args->wigFile);
	if(args->wigFile != NULL && firstChar != -1)
	{		
        //Get output file.
        FILE* outFile = args->outFile;
        FILE* outBed = args->bedFile;
        
        //Variables that will be used when traversing file.
        int intervalStart = 0;
        int intervalEnd = 0;
        int prevEnd = 0;
        int chromPos;
        char input[100];
        char input_bed[500];
        char* signalStr;
        float signal;
        float signalList[8000];	
        char* atEnd = chrom;
        char* atEnd_bed = chrom;
        char* sub_str1;
        char* sub_str2;
        char* annotation;

        //Navigate to the beginning of the WIG file.
        //Navigate to the end of the output file to append.
        rewind(args->wigFile);
        rewind(args->bedInput);
        fseek(outFile, 0, SEEK_END);
        
        //Check for EOF.
        char c;
        c = fgetc(args->wigFile);
        ungetc(c, args->wigFile);
        if(c == -1)
        {
            atEnd = NULL;
        }
        
        //We don't care about the wig header.
        char junk[501];
        fgets(junk, 500, args->wigFile);
        fgets(junk, 500, args->wigFile);
        
        //Extract the first signal.
        //Do not continue if the end of the file has been reached.
        atEnd_bed = fgets(input_bed, 500, args->bedInput);
        atEnd = fgets(input, 100, args->wigFile);
        int prevLine = 0;
        if(atEnd != NULL)
        {
            //Extract the first signal.
            signalStr = getTabLocation(input) + 1;
            signal = atof(signalStr);
            
            //Extract the first chromosome position.		
            (*getTabLocation(input)) = '\0';
            chromPos = atoi(input);
            
            //Get the first interval.
            sub_str1 = getTabLocation(input_bed) + 1;
            (*getTabLocation(input_bed)) = '\0';
            sub_str2 = getTabLocation(sub_str1) + 1;
            (*getTabLocation(sub_str1)) = '\0';
            annotation = getTabLocation(sub_str2) + 1;
            (*getTabLocation(sub_str2)) = '\0';
            intervalStart = atoi(sub_str1);
            intervalEnd = atoi(sub_str2);
        }
        
        //Go through until reaching the end of the file.
        //Fill in the counts for each signal level.
        int windowSize = 8000;
        int binSize = 50;
        while(atEnd_bed != NULL && atEnd != NULL)
        {
            //We only want to count regions that do not overlap and that have at least 3000 bp.
            if(intervalStart >= prevEnd)
            {  
                //Initialize signal data for region.
                int j;
                for(j = 0; j < windowSize / binSize; j++)
                {
                    signalList[j] = -1;
                }
                
                //Search for the marker that is equivalent to the first number in the interval.
                while(chromPos <= intervalStart && atEnd != NULL)
                {
                    prevLine = chromPos;
                    
                    //Get the next line.
                    //Do not continue if the end of the file has been reached.
                    atEnd = fgets(input, 100, args->wigFile);

                    if(atEnd != NULL)
                    {
                        //Extract the signal.
                        signalStr = getTabLocation(input) + 1;
                        signal = atof(signalStr);
                
                        //Extract the chromosome position.		
                        (*getTabLocation(input)) = '\0';
                        chromPos = atoi(input);
                    }
                }
                
                //If the end of the file was not reached, iteratively add every element in
                //the list until the end of the interval. Stop if the gap limit is reached or
                //if there is a gap at the beginning of the region.
                int signalPos = 0;
                while (atEnd != NULL && chromPos <= intervalEnd)
                {
                    //Fill in signal
                    signalList[signalPos] = signal;
                    prevLine = chromPos;
                    signalPos++;
                    
                    //Extract the next input and check that the end of the file
                    //has not been reached.
                    atEnd = fgets(input, 100, args->wigFile);
                    if(atEnd != NULL)
                    {
                        //Extract the signal.
                        signalStr = getTabLocation(input) + 1;
                        signal = atof(signalStr);

                        //Extract the next chromosome position.		
                        (*getTabLocation(input)) = '\0';
                        chromPos = atoi(input);
                    }
                }

                //If the end was reached, then print out signals.
                //Print out the chromosome and the interval.
                char intStart[51];
                char intEnd[51];
                char peakCount[10];
                sprintf(intStart, "%d", intervalStart);
                sprintf(intEnd, "%d", intervalEnd);

                //Print out the signals.
                for(j = 0; j < windowSize / binSize && signalList[j] != -1; j++)
                {
                    char sigVal[51];
                    sprintf(sigVal, "%f", signalList[j]);
                    fputs(sigVal, outFile);
                    if(j < windowSize / binSize - 1 && signalList[j+1] != -1)
                    {
                        fputs(",", outFile);
                    }
                }
                fputs("\n", outFile);
                
                //Print the line to the BED file.
                fputs("chr", outBed);
                fputs(chrom, outBed);
                fputs("\t", outBed);
                fputs(intStart, outBed);
                fputs("\t", outBed);
                fputs(intEnd, outBed);
                fputs("\t", outBed);
                fputs(annotation, outBed);
            }
            else{
                printf("%d %d\n", intervalStart, prevEnd);
            }
            //Increment the interval values.
            atEnd_bed = fgets(input_bed, 500, args->bedInput);
            if(atEnd_bed)
            {
                prevEnd = intervalEnd;
                sub_str1 = getTabLocation(input_bed) + 1;
                (*getTabLocation(input_bed)) = '\0';
                sub_str2 = getTabLocation(sub_str1) + 1;
                (*getTabLocation(sub_str1)) = '\0';
                annotation = getTabLocation(sub_str2) + 1;
                (*getTabLocation(sub_str2)) = '\0';
                intervalStart = atoi(sub_str1);
                intervalEnd = atoi(sub_str2);
            }
		}
    }
}

//Remove null characters that sometimes randomly end up in the file.
//Then find the tab location before the start of the signal.
char* getTabLocation(char* string)
{
	char* ptr;
	for(ptr = string; *ptr != '\n'; ptr++)
	{
		if(*ptr == '\0')
		{
			int i;
			for(i = 1; ptr + i == '\0'; i++);
			*ptr = *(ptr + i);
		}
	}
	for(ptr = string; *ptr != '\0'; ptr++)
	{
		if(*ptr == '\t')
		{
			return ptr;
		}
	}
	return NULL;
}
		