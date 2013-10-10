/* 
 * File:   readPoints.h
 * Author: xy
 *
 * Created on September 22, 2013, 9:06 PM
 */

#ifndef READPOINTS_H
#define	READPOINTS_H

typedef struct DataStruct {
    float xin;
    float yin;
    float xref;
    float yref;
} DataStruct;

#define MAX_LINE_LENGTH 1024

int readPoints(const char *fName, DataStruct **points, int *pointNum);
int writePoints(const char *fName, DataStruct *points, int pointNum);

#endif	/* READPOINTS_H */

