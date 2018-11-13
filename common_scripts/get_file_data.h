#include<stdio.h>

typedef struct intervalArgs
{
	int resolution;
	float cutoff;
	char* training;
	int regionSize;
    int binSize;
	char* chrom;
	FILE* outFile;
	FILE* wigFile;
}PrintIntervalArgs;

void* printIntervals(void* args);
char* getTabLocation(char* string);