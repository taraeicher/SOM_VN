#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
#include "getFileData.h"

int counter;
int gapFractionAllowed = 0;

//main function
int main(int argc, char* argv[])
{
	//Initialize random number generator.
	srand(5000);
	
	//Get arguments.
	char* wigFileName = argv[1];
	int resolution = atoi(argv[2]);
	int cutoff = atoi(argv[3]);
	char* training = argv[4];
	char* chrom = argv[5];
	char* outFileDir = argv[6];
	
	//Open WIG file.
	FILE* wigFile = fopen(wigFileName, "r");
	
	//Initialize window sizes and output files.
	int windowSizes[] = {500, 1000, 2000, 4000, 8000, 16000, 32000};
	FILE* outputFiles[(int)(sizeof(windowSizes) / sizeof(windowSizes[0]))];
	int i;
	for(i = 0; i < (int)(sizeof(windowSizes) / sizeof(windowSizes[0])); i++)
	{
		char fileName[20];
		strcpy(fileName, outFileDir);
		strcat(fileName, "window");
		char windowSizeStr[2];
		sprintf(windowSizeStr, "%d", i);
		strcat(fileName, windowSizeStr);
		outputFiles[i] = fopen(fileName, "w");
	}
	
	//Output intervals to files for processing.
	PrintIntervalArgs printIntervalArgs;
	void* args3;
	printIntervalArgs.wigFile = wigFile;
	printIntervalArgs.resolution = resolution;
	printIntervalArgs.cutoff = cutoff;
	printIntervalArgs.training = training;
	printIntervalArgs.outFiles = &outputFiles[0];
	printIntervalArgs.windowSizeList = windowSizes;
	printIntervalArgs.numWindowSizes = sizeof(windowSizes) / sizeof(windowSizes[0]);
	printIntervalArgs.chrom = chrom;
	args3 = &printIntervalArgs;
	printIntervals(args3);
	
	//Close WIG file and output files.
	fclose(wigFile);

	for(i = 0; i < (int)(sizeof(windowSizes) / sizeof(windowSizes[0])); i++)
	{
		if(outputFiles[i] != NULL)
		{
			fclose(outputFiles[i]);
		}
	}
	
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

	//Navigate to the beginning of the file.
	rewind(args->wigFile);
	
	//Check for empty file.
	char firstChar;
	firstChar = fgetc(args->wigFile);
	ungetc(firstChar, args->wigFile);
	if(args->wigFile != NULL && firstChar != -1)
	{		
		//For each window size, traverse the entire WIG file and fill in data.
		int i;
		for(i = 0; i < args->numWindowSizes; i++)
		{
			//Get output file.
			FILE* outFile = args->outFiles[i];
			
			//Variables that will be used when traversing file.
			int windowSize = *(args->windowSizeList + i);			
			int windowGap = windowSize;
			int margin = (int)(windowSize / (4 * 50)) * 50; // 1/4 of the window size.
			int stepSize = windowSize + margin;
			int intervalStart = 0;
			int intervalEnd = 0;
			int chromPos;
			char input[100];
			char* signalStr;
			float signal;
			float signalList[windowSize];	
			char* atEnd = chrom;
			int stepsCount = 2;
			int steps;
			
			// Add margins for training data.
			if (*(args->training) == 'Y')
			{
				windowSize += 2 * margin;
			}	

			//If generating data to annotate, double the window size and step by the window size.
			else if(*(args->training) == 'N')
			{
				stepSize = windowSize;
				windowSize *= 2;
				windowGap = 0;
			}
			
			//Generate the data using a sliding window.
			for(steps = 0; steps < stepsCount; steps++)
			{
				//Navigate to the beginning of the WIG file.
				//Navigate to the end of the output file to append.
				rewind(args->wigFile);
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
					prevLine = chromPos - windowSize;
					
					//Get the first interval.
					intervalStart = chromPos;
					intervalEnd = intervalStart + windowSize;
					
					//Add the step size to the interval position if this is the second iteration.
					if(steps == 1)
					{
						intervalStart += stepSize;
						intervalEnd += stepSize;
					}
				}

				//Go through until reaching the end of the file.
				//Fill in the counts for each signal level.
				while(atEnd != NULL)
				{
					//Initialize signal data for region.
					int j;
					for(j = 0; j < windowSize / args->resolution; j++)
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
					int gapExists = 0;
					while (atEnd != NULL && chromPos <= intervalEnd && gapExists == 0)
					{
						//Fill in signal. Only fill in first slot if there is no gap.
						gapExists = 0;
						if(signalPos == 0 && (chromPos - intervalStart) / args->resolution > 1)
						{
							gapExists = 1;
						}
						else
						{
							signalList[signalPos] = signal;
							prevLine = chromPos;
							signalPos++;
						}
						
						//Extract the next input and check that the end of the file
						//has not been reached.
						atEnd = fgets(input, 100, args->wigFile);
						if(atEnd != NULL && gapExists == 0)
						{
							//Extract the signal.
							signalStr = getTabLocation(input) + 1;
							signal = atof(signalStr);

							//Extract the next chromosome position.		
							(*getTabLocation(input)) = '\0';
							chromPos = atoi(input);
						
							//If there is a gap, set gapExists to 1
							if ((chromPos - prevLine) / args->resolution > 1)
							{
								gapExists = 1;
							}
						}
					}
					
					//Check to see whether any -1 exists in the region. If so,
					//set gapExists to 1.
					int sigIndex = 0;
					for(; sigIndex < windowSize / args->resolution; sigIndex++)
					{
						if(signalList[sigIndex] < 0)
						{
							gapExists = 1;
						}
					}

					//If the end was reached, then print out non-zero signals without gaps.
					if(gapExists == 0)
					{
						//Check to see if the region contains anything above the cutoff.
						//If not, we don't need it.
						int belowCutoff = 1;
						int sig;
						for(sig = 0; sig < windowSize / args->resolution; sig++)
						{
							if(signalList[sig] > args->cutoff)
							{
								belowCutoff = 0;
							}
						}
						if(belowCutoff == 0)
						{
							//Print out the chromosome and the interval.	
							fputs(chrom, outFile);
							fputs(",", outFile);
							char intStart[51];
							char intEnd[51];
							char peakCount[10];
							sprintf(intStart, "%d", intervalStart);
							sprintf(intEnd, "%d", intervalEnd);
							fputs(intStart, outFile);
							fputs(",", outFile);
							fputs(intEnd, outFile);
							fputs(",", outFile);

							//Print out the signals.
							for(j = 0; j < windowSize / args->resolution; j++)
							{
								char sigVal[51];
								sprintf(sigVal, "%f", signalList[j]);
								fputs(sigVal, outFile);
								if(j < windowSize / args->resolution - 1)
								{
									fputs(",", outFile);
								}
							}
							fputs("\n", outFile);
						}
					}
					
					//Increment the interval values.
					intervalStart += windowSize + windowGap;
					intervalEnd += windowSize + windowGap;
				}	
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
		