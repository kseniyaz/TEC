#include "parserO.h"


//переход на следующую строку файла
static void nextLine(FILE* file) {
	char buffer[256];
	fgets(buffer, 256, file);
}

static bool getOFileParamsNav(FILE* obs, OFileParams* _O) {
	char buffer[82];
	char strVal[14];
	fgets(buffer, 82, obs);
	strlen(buffer);
	if (strlen(buffer) != 79) {
		nextLine(obs);
		return false;
	}
	memcpy(strVal, buffer + 0, 14);
	_O->c1 = atof(strVal);
	memcpy(strVal, buffer + 16, 14);
	_O->p1 = atof(strVal);
	memcpy(strVal, buffer + 32, 14);
	_O->l1 = atof(strVal);
	memcpy(strVal, buffer + 48, 14);
	_O->c2 = atof(strVal);
	memcpy(strVal, buffer + 65, 14);
	_O->p2 = atof(strVal);
	fgets(buffer, 82, obs);
	if (strlen(buffer) != 17) {
		return false;
	}
	memcpy(strVal, buffer + 0, 14);
	_O->l2 = atof(strVal);
	if (_O->c1 == 0 || _O->p1 == 0 || _O->l1 == 0 || _O->c2 == 0 || _O->p2 == 0 || _O->l2 == 0) {
		return false;
	}
	_O->satVisible = true;
	return true;
}


//выделение из буффера названий КА в общий массив
static void bufferToSatNameArray(char* buffer, char satNameArray[36][4], int size, int string) {
	for (int i = 0; i < size; i++) {
		satNameArray[string * 12 + i][3] = '\0';
		memcpy(satNameArray[string * 12 + i], buffer + 32 + 3 * i, 3);
	}
}




void parseOFile(StrVector* oFilenameList, int fileNumber, OFileVector* o) {
	Date fileDate = getFileDate(oFilenameList->container[fileNumber]);

	FILE* oFile = fopen(oFilenameList->container[fileNumber], "rt");

	char buffer[82];

	do {
		fgets(buffer, 81, oFile);
	} while (!strstr(buffer, "END OF HEADER"));

	fgets(buffer, 81, oFile);

	while (!feof(oFile)) {
		int sec;
		sscanf(buffer, "%*d %*d %*d %d %d %d", &fileDate.hours, &fileDate.minutes, &sec);
		fileDate.seconds = sec;
		char numberOfSatS[2] = { buffer[30], buffer[31] };
		const int numberOfSat = atoi(numberOfSatS);

		char satNameArray[36][4];

		bufferToSatNameArray(buffer, satNameArray, 12, 0);

		nextLine(oFile);
		fgets(buffer, 81, oFile);

		if (numberOfSat < 25) {
			bufferToSatNameArray(buffer, satNameArray, numberOfSat - 12, 1);
		}

		if (numberOfSat >= 25) {
			bufferToSatNameArray(buffer, satNameArray, 12, 1);
			fgets(buffer, 81, oFile);
			bufferToSatNameArray(buffer, satNameArray, numberOfSat - 24, 2);
		}

		OFileList _satList;
		OFileParams sat = { .c1 = 0, .p1 = 0, .l1 = 0, .c2 = 0, .p2 = 0, .l2 = 0, .satVisible = false };
		for (int i = 1; i <= GNUM; i++) {
			_satList.G[i] = sat;
		}
		for (int i = 0; i <= RNUM; i++) {
			_satList.R[i] = sat;
		}

		for (int i = 0; i < numberOfSat; i++) {


			if (!getOFileParamsNav(oFile, &sat)) continue;

			char  satNumberS[2] = { satNameArray[i][1], satNameArray[i][2] };
			if (satNameArray[i][0] == 'G') {
				_satList.G[atoi(satNumberS)] = sat;
			}
			else if (satNameArray[i][0] == 'R') {
				_satList.R[atoi(satNumberS)] = sat;
			}
		}
		OFile _o = { fileDate, _satList };
		OFilePush(o, &_o);
		fgets(buffer, 81, oFile);
	}
	fclose(oFile);
}