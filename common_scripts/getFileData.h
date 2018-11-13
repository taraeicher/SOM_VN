#include<stdio.h>

typedef struct intervalArgs
{
	int resolution;
	float cutoff;
	char* training;
	int numWindowSizes;
	char* chrom;
	FILE** outFiles;
	FILE* wigFile;
	int* windowSizeList;
}PrintIntervalArgs;

void* printIntervals(void* args);
char* getTabLocation(char* string);