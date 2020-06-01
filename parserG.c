#include "parserG.h"
//чтение maxCount символов файла в переменную типа double
static void fgetd(double* var, int maxCount, FILE* stream) {
	char* buffer;
	buffer = (char*)calloc(maxCount, sizeof(char));

	fgets(buffer, maxCount, stream);

	for (int j = 0; j < maxCount - 1; j++) {
		if (buffer[j] == 'D')   buffer[j] = 'e';
	}

	*(var) = atof(buffer);

	free(buffer);
}
//чтение maxCount символов файла в переменную типа int
static void fgeti(int* var, int maxCount, FILE* stream) {
	char* buffer;
	buffer = (char*)calloc(maxCount, sizeof(char));

	fgets(buffer, maxCount, stream);

	for (int j = 0; j < maxCount - 1; j++) {
		if (buffer[j] == 'D')   buffer[j] = 'e';
	}

	*(var) = atoi(buffer);

	free(buffer);
}

//переход на следующую строку файла
static void nextLine(FILE* file) {
	char buffer[256];
	fgets(buffer, 256, file);
}

//парсинг блока g файла
static void getGFileParams(FILE* nav, GFileParams* _N) {

	char buffer[82];
	fgets(buffer, 3, nav);
	fgetd(&_N->TauN, 20, nav);
	fgetd(&_N->GammaN, 20, nav);
	nextLine(nav);
	_N->satVisible = true;
	fgets(buffer, 4, nav);

	fgetd(&_N->X, 20, nav);
	fgetd(&_N->Vx, 20, nav);
	fgetd(&_N->jx, 20, nav);

	nextLine(nav);
	fgets(buffer, 4, nav);

	fgetd(&_N->Y, 20, nav);
	fgetd(&_N->Vy, 20, nav);
	fgetd(&_N->jy, 20, nav);
	fgeti(&_N->K, 20, nav);
	nextLine(nav);
	fgets(buffer, 4, nav);

	fgetd(&_N->Z, 20, nav);
	fgetd(&_N->Vz, 20, nav);
	fgetd(&_N->jz, 20, nav);
	_N->X *= 1000;
	_N->Vx *= 1000;
	_N->jx *= 1000;
	_N->Y *= 1000;
	_N->Vy *= 1000;
	_N->jy *= 1000;
	_N->Z *= 1000;
	_N->Vz *= 1000;
	_N->jz *= 1000;

	nextLine(nav);
}

 
//парсинг g файла аналогично n
void parseGFile(StrVector* gFilenameList, int fileNumber, GFileVector* g) {
	Date fileDate = getFileDate(gFilenameList->container[fileNumber]);

	FILE* gFile = fopen(gFilenameList->container[fileNumber], "r");

	char buffer[82];
	double TauC,leap;
	do {
		fgets(buffer, 81, gFile);
		if (strstr(buffer, "CORR TO SYSTEM TIME")) {
			char TauCS[20];
			memcpy(TauCS,buffer+21,20);
			TauCS[15]='e';
			TauC = atof(TauCS);
		}		
		if (strstr(buffer, "LEAP SECONDS")) {
			char leapS[3];
			memcpy(leapS, buffer +3,3);
			leap = atof(leapS);
		}
	} while (!strstr(buffer, "END OF HEADER"));

	int month, day, hrs, min, sec;
	int satNumber;

	
	while (true) {
		fgets(buffer, 21, gFile);
		if (feof(gFile)) break;

		sscanf(buffer, "%d %*d %*d %d %d %d %d", &satNumber, &day, &hrs, &min, &sec);

		if (day == fileDate.day) {
			break;
		}
		for (int i = 0; i < 4; i++) nextLine(gFile);
	}
	fileDate.hours = hrs;
	fileDate.minutes = min;
	fileDate.seconds = sec;

	GFile _g;
	for (int i = 1; i <= RNUM; i++) {
		_g.R[i].satVisible = false;
	}

	GFileParams _G;
	getGFileParams(gFile, &_G);
	_G.TauC = TauC;
	_G.leap = leap;
	_g.R[satNumber] = _G;

	while (true) {
		fgets(buffer, 21, gFile);
		if (feof(gFile)) break;

		sscanf(buffer, "%d %*d %d %d %d %d %d", &satNumber, &month, &day, &hrs, &min, &sec);
		getGFileParams(gFile, &_G);

		Date _time = { fileDate.year,month, day,hrs,min,sec };
		_G.TauC = TauC;
		_G.leap = leap;
		if (cmpDate(&fileDate, &_time)) {
			_g.R[satNumber] = _G;
		}else {
			_g.time = fileDate;
			
			GFilePush(g, &_g);
			fileDate = _time;
			for (int i = 1; i <= RNUM; i++) {
				_g.R[i].satVisible = false;
			}
			_g.R[satNumber] = _G;
		}
	}

	fclose(gFile);

}