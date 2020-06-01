#ifndef _C_UTILS_H_
#define _C_UTILS_H_
#include <stdlib.h>
#include <math.h>
#include "Matrix.h"
#include "Vector.h"
#include "parsers.h"
#include "Date.h"
#include "errexit.h"

#include <windows.h>
#define Ra 6378137.0
#define Rb 6356752.3142
#define LSQITERMAX 10

#define SQRT(a) sqrt(a)
#define SQR(a) ((a)*(a))

#define max4(a,b,c,d) max(max((a),(b)),max((c),(d)))

void RK4(void(*der)(double, double*, double*), double* X1, double* X2, double t1, double t2, double h, const int  size);
Matrix* lsq(Matrix* A, Matrix* B);
void DecToSpher(double X, double Y, double Z, double* lat, double* lon);
void findAllFilename(StrVector* filenameList, const char* pattern);
double norm(int L, double* vec, int size);
double l2norm(double* vec, int size);
double linfnorm(double* vec, int size);

double linInter(double X[2], double F[2], double x);
double biLinInter(double X[2], double Y[2], double F[2][2], double x, double y);
double triLinInter(double X1[2], double X2[2], double X3[2], double F[2][2][2], double x1, double x2, double x3);
#endif