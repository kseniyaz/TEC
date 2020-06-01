#ifndef _C_VECTOR_H_
#define _C_VECTOR_H_
#include <stdlib.h>
/******************
VECTOR динамический массив
	VECTORSTRUCT каркас вектора типа X
		capacity размер выделенной памяти
		sizeCapacity размер выделяемой памяти
		size количество элементво в векторе
		container массив элементов типа Х
	VECTORINIT
	void NAME ## Init (NAME ## Vector *vector, int capacity)
		выделение первоначальной памяти capacity под вектор
	VECTORPUSH
	void NAME ## Push(NAME ## Vector *vector, X *value)
		добавление элемента в вектор и расширение выделенной памяти на sizeCapacity
	VECTORCOPY(X, NAME)
	void NAME ## Copy (NAME ## Vector *vector,NAME ## Vector *sourceVector)
		освобождение памяти vector и инициализация его парамерами и значениями sourceVector
******************/
#define VECTORSTRUCT(X, NAME) \
	typedef struct {\
		X *container; \
		int capacity;\
		int sizeCapacity;\
		int size; } \
	NAME ## Vector;

#define VECTORINIT(X, NAME) \
	inline void NAME ## Init (NAME ## Vector *vector, int capacity){\
		vector->size = 0;\
		vector->capacity = capacity;\
		vector->sizeCapacity = capacity;\
		vector->container=(X*)calloc(vector->capacity, sizeof(X));\
	}

#define VECTORPUSH(X, NAME) \
	inline void NAME ## Push(NAME ## Vector *vector, X *value) {\
        vector->container[vector->size]=*value;\
		vector->size++;\
		if (vector->size >= vector->capacity) {\
			vector->capacity += vector->sizeCapacity;\
			vector->container = (X*)realloc(vector->container, sizeof(X)*vector->capacity);\
		}\
	}

#define VECTORFREE(X, NAME)\
	inline void NAME ## Free (NAME ## Vector *vector){\
		vector->size = 0;\
		vector->capacity = 0;\
		vector->sizeCapacity = 0;\
		free(vector->container);\
	}

#define VECTORCOPY(X, NAME)\
	inline void NAME ## Copy (NAME ## Vector *vector,NAME ## Vector *sourceVector){\
		NAME ## Free(vector);\
		NAME ## Init(vector,sourceVector->sizeCapacity);\
		for(int i=0;i<sourceVector->size;i++){\
			NAME ## Push(vector,sourceVector->container+i);\
		}\
	}

#define VECTOR(X, NAME) VECTORSTRUCT(X, NAME) VECTORINIT(X, NAME) VECTORPUSH(X, NAME) VECTORFREE(X, NAME) VECTORCOPY(X, NAME)


typedef struct {
	char** container;
	int capacity;
	int sizeCapacity;
	int size;
	int strLength;
} StrVector;


 void StrInit(StrVector* strList, int capacity, int strLength);
 void StrPush(StrVector* strList, char* elem);
 void StrFree(StrVector* strList);

#endif