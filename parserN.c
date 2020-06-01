#include "parserN.h"

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

//переход на следующую строку файла
static void nextLine(FILE* file) {
	char buffer[256];
	fgets(buffer, 256, file);
}

void getNFileParams(FILE* nav, NFileParams* _G) {
	char buffer[82];
	fgets(buffer, 3, nav);
	fgetd(&_G->af0, 20, nav);
	fgetd(&_G->af1, 20, nav);
	fgetd(&_G->af2, 20, nav);
	
	nextLine(nav);
	_G->satVisible = true;
	fgets(buffer, 4, nav);

	fgets(buffer, 20, nav);

	fgetd(&_G->Crs, 20, nav);
	fgetd(&_G->delta_n, 20, nav);
	fgetd(&_G->M0, 20, nav);

	nextLine(nav);
	fgets(buffer, 4, nav);

	fgetd(&_G->Cuc, 20, nav);
	fgetd(&_G->e, 20, nav);
	fgetd(&_G->Cus, 20, nav);
	fgetd(&_G->A, 20, nav);

	nextLine(nav);
	fgets(buffer, 4, nav);

	fgetd(&_G->Toe, 20, nav);
	fgetd(&_G->Cic, 20, nav);
	fgetd(&_G->Omega0, 20, nav);
	fgetd(&_G->Cis, 20, nav);

	nextLine(nav);
	fgets(buffer, 4, nav);

	fgetd(&_G->i0, 20, nav);
	fgetd(&_G->Crc, 20, nav);
	fgetd(&_G->w0, 20, nav);
	fgetd(&_G->Omega_p, 20, nav);

	nextLine(nav);
	fgets(buffer, 4, nav);

	fgetd(&_G->IDOT, 20, nav);

	nextLine(nav);
	fgets(buffer, 4, nav);

	fgets(buffer, 20, nav);
	fgets(buffer, 20, nav);
	fgetd(&_G->Tgd, 20, nav);

	nextLine(nav);
	nextLine(nav);
}




void parseNFile(StrVector* nFilenameList, int fileNumber, NFileVector* n) {
	Date fileDate = getFileDate(nFilenameList->container[fileNumber]);

	FILE* nFile = fopen(nFilenameList->container[fileNumber], "r");

	char buffer[82];

	do {
		fgets(buffer, 81, nFile);
	} while (!strstr(buffer, "END OF HEADER"));

	int month, day, hrs, min, sec;
	int satNumber;

	while (true) {
		fgets(buffer, 21, nFile);
		if (feof(nFile)) break;

		sscanf(buffer, "%d %*d %*d %d %d %d %d", &satNumber, &day, &hrs, &min, &sec);

		if (day == fileDate.day && (sec == 0 || sec == 30)) {
			break;
		}
		for (int i = 0; i < 8; i++) nextLine(nFile);
	}

	fileDate.hours = hrs;
	fileDate.minutes = min;
	fileDate.seconds = sec;

	NFile _n;
	for (int i = 1; i <= GNUM; i++) {
		_n.G[i].satVisible = false;
	}

	NFileParams _N;
	getNFileParams(nFile, &_N);
	_n.G[satNumber] = _N;

	while (true) {
		fgets(buffer, 21, nFile);
		if (feof(nFile)) break;

		sscanf(buffer, "%d %*d %d %d %d %d %d", &satNumber, &month, &day, &hrs, &min, &sec);
		if ((sec != 0 && sec != 30)) {
			for (int i = 0; i < 8; i++) nextLine(nFile);
			continue;
		}

		getNFileParams(nFile, &_N);

		Date _time = { fileDate.year,month, day,hrs,min,sec };

		if (cmpDate(&fileDate, &_time)) {
			_n.G[satNumber] = _N;
		}
		else {
			_n.time = fileDate;
			NFilePush(n, &_n);
			fileDate = _time;
			for (int i = 1; i <= GNUM; i++) {
				_n.G[i].satVisible = false;
			}
			_n.G[satNumber] = _N;
		}
	}

	fclose(nFile);
}