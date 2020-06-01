#ifndef _C_PARSERI_H_
#define _C_PARSERI_H_

#include <stdio.h>
#include <stdbool.h>
#include "Date.h"
#include "Vector.h"
#include <string.h>

typedef struct {
	double lat[71];
	double lon[73];
	int TEC[71][73];
	double h;
	Date time;

} IonexFile;//значений ПЭС на конкретный момент времени, долготу и широту

VECTOR(IonexFile, IonexFile)

void parseIonex(StrVector* ionexFilenameList, int fileNumber, IonexFileVector* i);



#endif