#include "Vector.h"





 void StrInit(StrVector* strList, int capacity, int strLength) {
	strList->size = 0;
	strList->capacity = capacity;
	strList->sizeCapacity = capacity;
	strList->strLength = strLength;

	strList->container = (char**)calloc(strList->capacity, sizeof(char*));

	for (int i = 0; i < strList->capacity; i++) strList->container[i] = (char*)calloc(strLength, sizeof(char));
}

 void StrPush(StrVector* strList, char* elem) {
	strcpy(strList->container[strList->size], elem);

	strList->size++;

	if (strList->size >= strList->capacity) {
		strList->capacity += strList->sizeCapacity;
		strList->container = (char**)realloc(strList->container, sizeof(char*) * strList->capacity);
		for (int i = strList->size; i < strList->capacity; i++) strList->container[i] = (char*)calloc(strList->strLength, sizeof(char));
	}
}

 void StrFree(StrVector* strList) {
	for (int i = 0; i < strList->size; i++) free(strList->container[i]);
	free(strList->container);
}