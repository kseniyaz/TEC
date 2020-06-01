#include "Date.h"

bool cmpDate(Date* d1, Date* d2) {
	if (d1->seconds == d2->seconds &&
		d1->minutes == d2->minutes &&
		d1->hours == d2->hours &&
		d1->day == d2->day &&
		d1->month == d2->month &&
		d1->year == d2->year) 
        return true;
	return false;
}

 Date getFileDate(const char* name) {
	char* dot = strstr(name, ".");

	const char fileDayS[2] = { *(dot - 3), *(dot - 2) };
	const char fileMonthS[2] = { *(dot - 5), *(dot - 4) };
	const char fileYearS[2] = { *(dot + 1), *(dot + 2) };

	Date date = { .year = atoi(fileYearS),.month = atoi(fileMonthS), .day = atoi(fileDayS),.hours = 0,.minutes = 0,.seconds = 0 };
	return date;
}


double epoch2time(Date ep)
 {
     const int doy[] = { 1,32,60,91,121,152,182,213,244,274,305,335 };
     double time = 0;
     int days, sec, year = ep.year+2000, mon = ep.month, day = ep.day;

     if (year < 1970 || 2099 < year || mon < 1 || 12 < mon) return time;

     days = (year - 1970) * 365 + (year - 1969) / 4 + doy[mon - 1] + day - 2 + (year % 4 == 0 && mon >= 3 ? 1 : 0);
     time = days * 86400. + ep.hours * 3600. + ep.minutes * 60. + ep.seconds;
     return time;
 }


 double time2gpst(double t, int* week)
{
    double t0 = epoch2time((Date) { -20, 1, 6, 0, 0, 0 });
    double sec = t - t0;
    int w = (int)(sec / (86400 * 7));

    if (week) *week = w;
    return (double)(sec - (double)w * 86400 * 7);
}
