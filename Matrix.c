#include "Matrix.h"
#include "errexit.h"

Matrix* MatrixInit(Matrix* m, int rows, int columns) {
	m->rows = rows;
	m->columns = columns;
	m->_el = (double**)malloc(rows*sizeof(double*));
	for (int i = 0; i < rows; i++) {
		m->_el[i]= (double*)calloc(columns, sizeof(double));
	}
	return m;
}

Matrix* MatrixCreate(int rows, int columns) {
	Matrix* m= (Matrix*)malloc(sizeof(Matrix));
	return MatrixInit(m,rows,columns);
}

Matrix* MatrixCreateSq(int size) {
	return MatrixCreate( size, size);
}

double* MatrixGetEl(Matrix* m, int row, int column) {
	errexit(row >= m->rows || column >= m->columns, "ERROR(MatrixGetEl):index out of bounds");
	return m->_el[row]+column;
}

double MatrixGetV(Matrix* m, int row, int column) {
	errexit(row >= m->rows || column >= m->columns, "ERROR(MatrixGetV):index out of bounds");
	return *MatrixGetEl(m, row, column);
}

double* MatrixGetRow(Matrix* m, int row) {
	errexit(row >= m->rows , "ERROR(MatrixGetRow):row out of bounds");
	double* _row = (double*)malloc(m->columns*sizeof(double));
	memcpy(_row, m->_el[row], m->columns * sizeof(double));
	return _row;
}

double* MatrixGetColumn(Matrix* m, int column) {
	errexit(column >= m->rows, "ERROR(MatrixGetColumn):column out of bounds");
	double* _column = (double*)malloc(m->rows*sizeof(double));
	for (int row = 0; row < m->rows; row++) {
		_column[row] = m->_el[row][column];
	}
	return _column;
}

Matrix* MatrixSet(Matrix* m, int row, int column,double value) {
	errexit(row >= m->rows || column >= m->columns, "ERROR(MatrixSet):index out of bounds");
	m->_el[row][column] = value;
	return m;
}

Matrix* MatrixSetRow(Matrix* m, int row, double *rowA) {
	errexit(row >= m->rows , "ERROR(MatrixSetRow):row out of bounds");
	memcpy(m->_el[row], rowA, m->columns * sizeof(double));
	return m;
}

Matrix* MatrixSetColumn(Matrix* m, int column, double* columnA) {
	errexit(column >= m->columns, "ERROR(MatrixSetColumn):column out of bounds");
	for (int row = 0; row < m->rows; row++) {
		m->_el[row][column] = columnA[row];
	}
	return m;
}

Matrix* MatrixSetMatrix(Matrix* m1, Matrix* m2) {
	errexit(m1->rows != m2->rows|| m1->columns != m2->columns, "ERROR(MatrixSetMatrix):size of matrix no equal");
	for (int row = 0; row < m1->rows; row++) {
		memcpy(m1->_el[row], m2->_el[row], m1->columns * sizeof(double));
	}
	return m1;
}

Matrix* MatrixSetArrOfArr(Matrix* m, double** rowAA) {
	for (int row = 0; row < m->rows; row++) {
		memcpy(m->_el[row], rowAA[row], m->columns * sizeof(double));
	}
	return m;
}

Matrix* MatrixSum(Matrix* m1, Matrix* m2) {
	errexit(m1->rows != m2->rows || m1->columns != m2->columns, "ERROR(MatrixSum):size of matrix no equal");
	Matrix* m = MatrixCreate(m1->rows, m1->columns);
	for (int row = 0; row < m->rows; row++) {
		for (int column = 0; column < m->columns; column++) {
			MatrixSet(m, row, column, MatrixGetV(m1, row, column) + MatrixGetV(m2, row, column));
		}
	}
	return m;
}

Matrix* MatrixSumRow(Matrix* m1,int rowN, double* row) {
	Matrix* m = MatrixCreate(m1->rows, m1->columns);
	for (int column = 0; column < m->columns; column++) {
		MatrixSet(m, rowN, column, MatrixGetV(m1, rowN, column) + row[rowN]);
	}
	return m;
}

Matrix* MatrixSumColumn(Matrix* m1, int columnN, double* column) {
	Matrix* m = MatrixCreate(m1->rows, m1->columns);
	for (int row = 0; row < m->rows; row++) {
		MatrixSet(m, row, columnN, MatrixGetV(m1, row, columnN) + column[columnN]);
	}
	return m;
}

Matrix* MatrixSub(Matrix* m1, Matrix* m2) {
	Matrix* m = MatrixCreate(m1->rows, m1->columns);
	for (int row = 0; row < m->rows; row++) {
		for (int column = 0; column < m->columns; column++) {
			MatrixSet(m, row, column, MatrixGetV(m1, row, column) - MatrixGetV(m2, row, column));
		}
	}
	return m;
}

Matrix* MatrixSubRow(Matrix* m1, int rowN, double* row) {
	Matrix* m = MatrixCreate(m1->rows, m1->columns);
	for (int column = 0; column < m->columns; column++) {
		MatrixSet(m, rowN, column, MatrixGetV(m1, rowN, column) - row[rowN]);
	}
	return m;
}

Matrix* MatrixSubColumn(Matrix* m1, int columnN, double* column) {
	Matrix* m = MatrixCreate(m1->rows, m1->columns);
	for (int row = 0; row < m->rows; row++) {
		MatrixSet(m, row, columnN, MatrixGetV(m1, row, columnN) - column[columnN]);
	}
	return m;
}

bool MatrixEq(Matrix* m1, Matrix* m2) {
	for (int row = 0; row < m1->rows; row++) {
		for (int column = 0; column < m1->columns; column++) {
			if (MatrixGetV(m1, row, column) != MatrixGetV(m2, row, column)) return false;
		}
	}
	return true;
}


Matrix* MatrixMul (Matrix* m1, Matrix* m2) {
	Matrix* m = MatrixCreate(m1->rows, m2->columns);
	for (int row = 0; row < m->rows; row++) {
		for (int column = 0; column < m->columns; column++) {			
			double mv = 0;
			for (int colrow = 0; colrow < m1->columns; colrow++) {
				double m1v = MatrixGetV(m1, row, colrow);
				double m2v = MatrixGetV(m2, colrow, column);
				mv += m1v * m2v;
			}
			MatrixSet(m, row, column, mv);
		}
	}
	return m;
}

Matrix* MatrixMulSc(Matrix* m1, double val) {
	Matrix* m = MatrixCreate(m1->rows, m1->columns);
	for (int row = 0; row < m->rows; row++) {
		for (int column = 0; column < m->columns; column++) {
			MatrixSet(m, row, column, val * MatrixGetV(m1, row, column));
		}
	}
	return m;
}

Matrix* MatrixTr(Matrix* m) {
	Matrix* m1 = MatrixCreate( m->columns, m->rows);
	for (int row = 0; row < m->rows; row++) {
		for (int column = 0; column < m->columns; column++) {
			MatrixSet(m1, column, row, MatrixGetV(m, row, column));
		}
	}
	return m1;
}


Matrix* MatrixInv(Matrix* m) {
	Matrix* mExt = MatrixCreate(m->rows, m->columns * 2);
	Matrix* m1 = MatrixCreateSq(m->rows);

	for (int row = 0; row < m->rows; row++) {
		for (int column = 0; column < m->columns; column++) {
			MatrixSet(mExt, row, column, MatrixGetV(m, row, column));
			MatrixSet(mExt, row, column + m->columns, (row == column) ? 1 : 0);
		}
	}

	for (int N = 0; N < mExt->rows-1; N++) {
		for (int row = N; row < mExt->rows; row++) {
			double elem = MatrixGetV(mExt, row, N);
			for (int col = mExt->columns - 1; col >= 0; col--) {
				MatrixSet(mExt, row, col, MatrixGetV(mExt, row, col) / elem);
			}
		}
		for (int row = N + 1; row < mExt->rows; row++) {
			for (int col = N; col < mExt->columns; col++) {
				MatrixSet(mExt, row, col, MatrixGetV(mExt, row, col)-MatrixGetV(mExt, N, col)  );
			}
		}
	}
	for (int N = mExt->rows - 1; N > 0; N--) {
		for (int row = N; row >=0; row--) {
			double elem = MatrixGetV(mExt, row, N);
			for (int col = mExt->columns - 1; col >= 0; col--) {
				MatrixSet(mExt, row, col, MatrixGetV(mExt, row, col) / elem);
			}
		}
		for (int row = N-1; row >= 0; row--) {
			for (int col = N; col < mExt->columns; col++) {
				MatrixSet(mExt, row, col, MatrixGetV(mExt, row, col) - MatrixGetV(mExt, N, col));
			}
		}
	}

	for (int col = mExt->columns - 1; col >= 0; col--) {
		MatrixSet(mExt, 0, col, MatrixGetV(mExt, 0, col) / MatrixGetV(mExt, 0, 0));
	}

	for (int row = 0; row < m->rows; row++) {
		for (int column = 0; column < m->columns; column++) {
			MatrixSet(m1, row, column, MatrixGetV(mExt, row, column+m->columns));
		}
	}
	MatrixFree(mExt);
	free(mExt);
	return m1;
}

void MatrixFree(Matrix* m) {	
	for (int row = 0; row < m->rows; row++) {
		free(m->_el[row]);
	}
	free(m->_el);
}