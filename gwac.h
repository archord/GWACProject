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

#define RESIDUALS1
#define GWAC_TEST1

#define MAX_LINE_LENGTH 1024

#define GWAC_SUCCESS 0
#define GWAC_ERROR 3001
#define GWAC_OPEN_FILE_ERROR 3002
#define GWAC_MALLOC_ERROR 3003
#define GWAC_FUNCTION_INPUT_NULL 3004
#define GWAC_FUNCTION_INPUT_EMPTY 3005
#define GWAC_STATUS_STR_NULL 3006

#define CHECK_RETURN_SATUS(var) {if(var!=GWAC_SUCCESS)return var;}
#define CHECK_STATUS_STR_IS_NULL(statusstr) \
        {if(statusstr==NULL)return GWAC_STATUS_STR_NULL;}
#define CHECK_INPUT_IS_NULL(fun,var,varname) \
        {if(var==NULL){\
            sprintf(statusstr, "Error Code: %d\n"\
                "In function %s, the input parameter \"%s\" is NULL!\n", \
                GWAC_FUNCTION_INPUT_NULL, fun, varname);\
            return GWAC_FUNCTION_INPUT_NULL;}}
#define CHECK_MALLOC_IS_NULL(fun,var,varname) \
        {if(var==NULL){\
            sprintf(statusstr, "Error Code: %d\n"\
                "In function %s, melloc memory for \"%s\" error!\n", \
                GWAC_MALLOC_ERROR, fun, varname);\
            return GWAC_MALLOC_ERROR;}}
#define CHECK_OPEN_FILE(fun,fp,fname) \
        {if(fp==NULL){\
            sprintf(statusstr, "Error Code: %d\n"\
                "In function %s, open file \"%s\" error!\n", \
                GWAC_OPEN_FILE_ERROR, fun, fname);\
            return GWAC_OPEN_FILE_ERROR;}}

int Gwac_geomap(vector<ST_STARPEER> matchpeervec,
        unsigned int order,
        unsigned int iter,
        float rejsigma,
        float &xrms,
        float &yrms,
        const char outfilename[],
        char statusstr[]);
int cofun(double x1, 
        double x2, 
        double *afunc, 
        int cofNum);
int GetShift(const char transfilename[], 
        float &xshift, 
        float &yshift,
        char statusstr[]);
int Gwac_geoxytran(vector<ST_STAR> &objvec,
        const char transfilename[],
        int flag,
        char statusstr[]);
int testGeomap();
void testGeoxytran();
void test2To5();

#endif	/* GWAC_H */

