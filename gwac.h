/* 
 * File:   gwac.h
 * Author: xy
 *
 * Created on September 9, 2013, 9:10 AM
 */

#ifndef GWAC_H
#define	GWAC_H


#include <iostream>
#include <vector>
#include <string>

using namespace std;

typedef struct {
    float x;
    float y;
    float ra;
    float dec;
    float flux;
    float mag;
    float fwhm;
    unsigned int pixnum;
    int flags;
    float thresh;
} ST_STAR;

typedef struct {
    ST_STAR obj;
    ST_STAR ref;
} ST_STARPEER;

#define RESIDUALS
#define GWAC_TEST

#define MAX_LINE_LENGTH 1024

#define GWAC_SUCCESS 0
#define GWAC_ERROR 3001
#define GWAC_OPEN_FILE_ERROR 3002
#define GWAC_MALLOC_ERROR 3003

int Gwac_geomap(vector<ST_STARPEER> matchpeervec,
        unsigned int order,
        unsigned int iter,
        float rejsigma,
        float &xrms,
        float &yrms,
        const char outfilename[],
        char statusstr[]);
int testGeomap();

#endif	/* GWAC_H */

