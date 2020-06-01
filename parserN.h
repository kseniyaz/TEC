#ifndef _C_PARSERN_H_
#define _C_PARSERN_H_

#include <stdio.h>
#include <stdbool.h>
#include "Date.h"
#include "Vector.h"

#define GNUM 32
#define RNUM 24


typedef struct {
	double af0, af1, af2;
	double Crs;
	double delta_n;
	double M0;

	double Cuc;
	double e;
	double Cus;
	double A;

	double Toe;
	double Cic;
	double Omega0;
	double Cis;

	double i0;
	double Crc;
	double w0;
	double Omega_p;

	double IDOT;

	double Tgd;

	bool satVisible;
} NFileParams;//структура содержимого n-файла конретоного КА

typedef struct {
	Date time;
	NFileParams G[GNUM + 1];
} NFile;



VECTOR(NFile, NFile)

void parseNFile(StrVector* nFilenameList, int fileNumber, NFileVector* n);

#endif