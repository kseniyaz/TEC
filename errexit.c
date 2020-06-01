#include "errexit.h"

void errexit(bool condition,const char * message) {
	if (condition) {
		printf(message);
		printf("\n");
		system("pause");
		_Exit(1);
	}
}