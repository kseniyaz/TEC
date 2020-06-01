#include "parserI.h"
#include <stdio.h>
#include <stdbool.h>
#include "Date.h"
#include "Vector.h"

//переход на следующую строку файла
static void nextLine(FILE* file) {
	char buffer[256];
	fgets(buffer, 256, file);
}

//парсинг ionex файлов
void parseIonex(StrVector* ionexFilenameList, int fileNumber, IonexFileVector* i) {
	Date fileDate;
	char* dot = strstr(ionexFilenameList->container[fileNumber], ".");

	const char fileDayS[3] = { *(dot - 4), *(dot - 3), *(dot - 2) };
	const char fileYearS[2] = { *(dot + 1), *(dot + 2) };

	const int _day = atoi(fileDayS);

	fileDate.year = atoi(fileYearS);
	/**
	int day = 0, month;
	//перевод из дн€ с начала года в мес€ц и день (в названии день с начала года)
	for (month = 1; month < 13; month++) {
		char monthDay;
		TIMECONV_GetNumberOfDaysInMonth(2000 + fileDate.year, month, &monthDay);
		if (day + monthDay >= _day) break;
		day += monthDay;
	}

	fileDate.month = month;
	/**/
	FILE* ionexFile = fopen(ionexFilenameList->container[fileNumber], "r");

	char buffer[82];

	do {
		fgets(buffer, 82, ionexFile);
	} while (!strstr(buffer, "END OF HEADER"));

	nextLine(ionexFile);

	while (true) {
		fgets(buffer, 82, ionexFile);
		if (strstr(buffer, "END OF FILE")) break;
		int sec;
		sscanf(buffer, "%*i %i %i %i %i %i", &fileDate.month, &fileDate.day, &fileDate.hours, &fileDate.minutes, &sec);
		fileDate.seconds = sec;
		nextLine(ionexFile);
		fgets(buffer, 82, ionexFile);
		IonexFile _i;
		// парсинг блоков ѕЁ— по широтам 
		for (int k = 0; k < 71; k++) {
			//парсинг значений ѕЁ— на конкретной широте по долготам
			for (int n = 0; n < 4; n++) {
				sscanf(buffer, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", &_i.TEC[k][n * 16 + 0], &_i.TEC[k][n * 16 + 1], &_i.TEC[k][n * 16 + 2], &_i.TEC[k][n * 16 + 3], &_i.TEC[k][n * 16 + 4], &_i.TEC[k][n * 16 + 5], &_i.TEC[k][n * 16 + 6], &_i.TEC[k][n * 16 + 7], &_i.TEC[k][n * 16 + 8], &_i.TEC[k][n * 16 + 9], &_i.TEC[k][n * 16 + 10], &_i.TEC[k][n * 16 + 11], &_i.TEC[k][n * 16 + 12], &_i.TEC[k][n * 16 + 13], &_i.TEC[k][n * 16 + 14], &_i.TEC[k][n * 16 + 15]);
				fgets(buffer, 82, ionexFile);
			}
			sscanf(buffer, "%d %d %d %d %d %d %d %d %d", &_i.TEC[k][64], &_i.TEC[k][65], &_i.TEC[k][66], &_i.TEC[k][67], &_i.TEC[k][68], &_i.TEC[k][69], &_i.TEC[k][70], &_i.TEC[k][71], &_i.TEC[k][72]);
			nextLine(ionexFile);
			fgets(buffer, 82, ionexFile);
		}
		_i.time = fileDate;
		_i.h = 450;
		for (int k = 0; k < 73; k++) _i.lon[k] = -180 + k * 5.0;
		for (int k = 0; k < 71; k++) _i.lat[k] = 87.5 - k * 2.5;
		IonexFilePush(i, &_i);
		int N;
		sscanf(buffer, "%i", &N);
		if (!strstr(buffer, "START OF TEC MAP")) break;
	}
	fclose(ionexFile);
}