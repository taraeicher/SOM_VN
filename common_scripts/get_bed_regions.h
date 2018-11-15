#include<stdio.h>

typedef struct intervalArgs
{
	int resolution;
	char* chrom;
	FILE* outFile;
    FILE* bedFile;
    FILE* bedInput;
	FILE* wigFile;
}PrintIntervalArgs;

void* printIntervals(void* args);
char* getTabLocation(char* string);