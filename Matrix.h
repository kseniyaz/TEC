#ifndef _C_MATRIX_H_
#define _C_MATRIX_H_
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

typedef struct {
	double** _el;
	int rows, columns;
} Matrix;

Matrix* MatrixInit(Matrix* m, int rows, int columns);
Matrix* MatrixCreate(int rows, int columns);
Matrix* MatrixCreateSq(int size);

double* MatrixGetEl(Matrix* m, int row, int column);
double  MatrixGetV(Matrix* m, int row, int column);
double* MatrixGetRow(Matrix* m, int row);
double* MatrixGetColumn(Matrix* m, int column);

Matrix* MatrixSet(Matrix* m, int row, int column, double value);
Matrix* MatrixSetRow(Matrix* m, int row, double* rowA);
Matrix* MatrixSetColumn(Matrix* m, int column, double* columnA);
Matrix* MatrixSetMatrix(Matrix* m1, Matrix* m2);
Matrix* MatrixSetArrOfArr(Matrix* m, double** rowAA);

Matrix* MatrixSum(Matrix* m1, Matrix* m2);
Matrix* MatrixSumRow(Matrix* m1, int rowN, double* row);
Matrix* MatrixSumColumn(Matrix* m1, int columnN, double* column);

Matrix* MatrixSub(Matrix* m1, Matrix* m2);
Matrix* MatrixSubRow(Matrix* m1, int rowN, double* row);
Matrix* MatrixSubColumn(Matrix* m1, int columnN, double* column);

bool MatrixEq(Matrix* m1, Matrix* m2);

Matrix* MatrixMul(Matrix* m1, Matrix* m2);
Matrix* MatrixMulSc(Matrix* m1, double val);

Matrix* MatrixTr(Matrix* m);

Matrix* MatrixInv(Matrix* m);

void MatrixFree(Matrix* m);

#endif