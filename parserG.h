#ifndef _C_PARSERG_H_
#define _C_PARSERG_H_
#include <stdio.h>
#include <stdbool.h>
#include "Date.h"
#include "Vector.h"

#define GNUM 32
#define RNUM 24

typedef struct {
	double X;
	double Y;
	double Z;
	double Vx;
	double Vy;
	double Vz;
	double jx;
	double jy;
	double jz;
	double TauN;
	double GammaN;
	double leap;
	int K;
	double TauC;
	bool satVisible;
} GFileParams;//структура содержимого g-файла конретоного КА

typedef struct {
	
	Date time;
	GFileParams R[RNUM + 1];
} GFile;//структура содержимого g-файла на конкретный момент времени

VECTOR(GFile, GFile)

void parseGFile(StrVector* gFilenameList, int fileNumber, GFileVector* g);

#endif