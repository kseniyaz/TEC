#ifndef _C_PARSERO_H_
#define _C_PARSERO_H_
#include <stdio.h>
#include <stdbool.h>
#include "Date.h"
#include "Vector.h"
#include <string.h>

#define GNUM 32
#define RNUM 24



typedef struct {
	double p1;
	double p2;
	double l1;
	double l2;
	double c1;
	double c2;
	bool satVisible;
} OFileParams;//структура содержимого o-файла конретоного КА

typedef struct {
	OFileParams G[GNUM + 1];
	OFileParams R[RNUM + 1];
} OFileList;//структура содержимого o-файла всех КА 

typedef struct {
	Date time;
	OFileList sat;
} OFile;

VECTOR(OFile, OFile)

void parseOFile(StrVector* oFilenameList, int fileNumber, OFileVector* o);


#endif