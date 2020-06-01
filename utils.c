#include "utils.h"

/**********
метод Рунге_Кутты 4 порядка
void(*der)(double, double *, double*)
контракт системы ОДУ где первый аргумент независимый параметр, воторой значения переменных, третий значения производных
X1 начальные условия
X2 результат интегрирования
t1 начальный момент
t2 конечный момент
h шаг интегрирования
size количество уравнений в системе
***********/
void RK4(void(*der)(double, double*, double*), double* X1, double* X2, double t1, double t2, double h, const int  size) {
	double* k1, * k2, * k3, * k4, * X;
	k1 = (double*)calloc(size, sizeof(double));
	k2 = (double*)calloc(size, sizeof(double));
	k3 = (double*)calloc(size, sizeof(double));
	k4 = (double*)calloc(size, sizeof(double));
	X = (double*)calloc(size, sizeof(double));

	memcpy(X2, X1, size * sizeof(double));

	while (t1 < t2) {
		memcpy(X, X2, size * sizeof(double));

		der(t1, X2, k1);

		for (int i = 0; i < size; i++) {
			X[i] = X2[i] + h * k1[i] / 2;
		}
		der(t1 + h / 2, X, k2);

		for (int i = 0; i < size; i++) {
			X[i] = X2[i] + h * k2[i] / 2;
		}
		der(t1 + h / 2, X, k3);

		for (int i = 0; i < size; i++) {
			X[i] = X2[i] + h * k3[i];
		}
		der(t1 + h, X, k4);

		for (int i = 0; i < size; i++) {
			X2[i] = X2[i] + h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
		}
		t1 += h;

		if (t2 - t1 < h) h = t2 - t1;//последний шаг если он меньше h равен оставшемуся  времени до конца
	}

	free(k1); free(k2); free(k3); free(k4); free(X);
}


Matrix* lsq(Matrix* A, Matrix* B) {
	Matrix* AT = MatrixTr(A);
	Matrix* ATA = MatrixMul(AT, A);
	Matrix* ATAInv = MatrixInv(ATA);
	Matrix* AtB = MatrixMul(AT, B);
	Matrix* x = MatrixMul(ATAInv, AtB);

	MatrixFree(AT); MatrixFree(ATA); MatrixFree(ATAInv); MatrixFree(AtB);
	free(AT); free(ATA); free(ATAInv); free(AtB);
	return x;
}


//перевод из прямоугольных координат в геодезические
void DecToSpher(double X, double Y, double Z, double* lat, double* lon) {
	double e = (Ra * Ra - Rb * Rb) / (Ra * Ra);
	double e_ = (Ra * Ra - Rb * Rb) / (Rb * Rb);

	double tet = atan2(Ra * Z, Rb * sqrt(X * X + Y * Y));
	*lat = atan2(Z + e_ * e_ * Rb * pow(sin(tet), 3), sqrt(X * X + Y * Y) - e * e * Ra * pow(cos(tet), 3));
	*lon = atan2(Y, X);
}


/*********************
findAllFilename поиск всех файлов по заданному шаблону(pattern) в текущей папке и добавление в filenameList
*************************/
void findAllFilename(StrVector* filenameList, const char* pattern) {
	WIN32_FIND_DATA findFileData;//переменная содержащая свойства найденного файла
	WIN32_FIND_DATA findFile;

	HANDLE oHandle = FindFirstFile((LPCSTR)pattern, &findFileData);//поиск необходимого файла
	HANDLE oHandleFile = FindFirstFile((LPCSTR)"*.o", &findFile);//поиск объектного файла для исключения его из результатов

	//если не найдены файлы то завершается функция
	if (oHandle != INVALID_HANDLE_VALUE) {
		//если найден объектный файл, то проверяется не равен ли он найденному по шаблону, если нет то найденный по шаблону добавляется в массив
		if (oHandleFile == INVALID_HANDLE_VALUE) StrPush(filenameList, findFileData.cFileName);
		else {
			if (strcmp(findFileData.cFileName, findFile.cFileName)) StrPush(filenameList, findFileData.cFileName);
			else {
				FindNextFile(oHandle, &findFileData);
				StrPush(filenameList, findFileData.cFileName);
			}
		}
	}
	else return;

	int i = 1;

	while (true) {
		FindNextFile(oHandle, &findFileData);
		// сравнивается с предыдущим найденным файлом, если они равны, то это последний файл и можно выйти из функции
		if (!strcmp(findFileData.cFileName, filenameList->container[i - 1])) break;
		//если не найден объектный файл, то файлы просто добавляются, иначе проверяются на каждой итерации 
		if (oHandleFile == INVALID_HANDLE_VALUE) {
			StrPush(filenameList, findFileData.cFileName);
		}
		else {
			if (strcmp(findFileData.cFileName, findFile.cFileName)) {
				StrPush(filenameList, findFileData.cFileName);
			}
			else {
				FindNextFile(oHandle, &findFileData);
				oHandleFile = INVALID_HANDLE_VALUE;
				StrPush(filenameList, findFileData.cFileName);
			}
		}
		i++;
	}
};


double norm(double L, double* vec, int size) {//L==-1 --> L inf
	double n = 0;
	if (L == -1) {
		for (int i = 0; i < size; i++) {
			n = vec[i] > n ? vec[i] : n;
		}
	} else {
		for (int i = 0; i < size; i++) {
			n += pow(fabs(vec[i]), L);
		}
		n = pow(n, 1 / L);
	}
	return n;
}

double l2norm(double* vec, int size) {
	return norm(2, vec, size);
}

double linfnorm(double* vec, int size) {
	return norm(-1, vec, size);
}

double linInter(double X[2], double F[2], double x) {
	return F[0] + (F[1] - F[0]) * (x - X[0]) / (X[1] - X[0]);
}

double biLinInter(double X1[2], double X2[2], double F[2][2], double x1, double x2) {
	double F_[2] = { linInter(X1, F[0], x1), linInter(X1, F[1], x1) };
	return linInter(X2, F_, x2);
}

double triLinInter(double X1[2], double X2[2], double X3[2], double F[2][2][2], double x1, double x2, double x3) {
	double F_[2] = { biLinInter(X1, X2, F[0], x1, x2),biLinInter(X1, X2, F[1], x1, x2) };
	return linInter(X3, F_, x3);
}