#ifndef _C_DATE_H_
#define _C_DATE_H_
#include <stdbool.h>
#include <string.h>

typedef struct {
	int year, month, day, hours, minutes;
	double seconds;
} Date;
bool cmpDate(Date* d1, Date* d2);
 Date getFileDate(const char* name);

 double epoch2time(Date ep);
 double time2gpst(double t, int* week);
#endif