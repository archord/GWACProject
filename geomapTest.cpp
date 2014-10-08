
#include <iostream>
#include <vector>
#include <string>

using namespace std;

#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gwac.h"
#include "readPoints.h"

float testFunction1(float x, float y) {
    float rst = 3.0 + 20.0 * x + 2.0 * x * x + 1.4 * y + 0.06 * x * y + 0.8 * y*y;
    int tmp = rand() % 20;
    if (tmp > 10) {
        tmp -= 10;
        rst = rst * (1.0 + 1.0 * tmp / 1000);
    } else {
        rst = rst * (1.0 - 1.0 * tmp / 1000);
    }
    return rst;
}

float testFunction2(float x, float y) {
    float rst = 2.0 + 30.5 * x + 8.0 * x * x + 0.008 * y + 0.8 * x * y + 10 * y*y;
    int tmp = rand() % 40;
    if (tmp > 20) {
        tmp -= 20;
        rst = rst * (1.0 + 1.0 * tmp / 1000);
    } else {
        rst = rst * (1.0 - 1.0 * tmp / 1000);
    }
    return rst;
}

int testGeomap() {

    unsigned int order = 5;
    unsigned int iter = 1;
    float rejsigma = 2.5;
    float xrms, yrms;
    float xrms2, yrms2;
    const char outfilename[MAX_LINE_LENGTH] = "gwac_geomap_result2.txt";
    char statusstr[MAX_LINE_LENGTH];
    memset(statusstr, 0, MAX_LINE_LENGTH);

    char *fName = "data";   //data  test.txt
    DataStruct *points;
    int pointNum;
    readPoints(fName, &points, &pointNum);

    vector<ST_STARPEER> matchpeervec;
    int i;
    for (i = 0; i < pointNum; i++) {
        ST_STARPEER peer;
        peer.ref.x = points[i].xref;
        peer.ref.y = points[i].yref;
        peer.obj.x = points[i].xin;
        peer.obj.y = points[i].yin;
        matchpeervec.push_back(peer);
    }

    //    int pointNum = 600;
    //    vector<ST_STARPEER> matchpeervec;
    //    int i;
    //    for (i = 0; i < pointNum; i++) {
    //        ST_STARPEER peer;
    //        int x = i + 10.32;
    //        int y = 0.2 * i + 1;
    //        peer.ref.x = x;
    //        peer.ref.y = y;
    //        peer.obj.x = testFunction1(x, y);
    //        peer.obj.y = testFunction2(x, y);
    //        matchpeervec.push_back(peer);
    //    }

    int rstStatus = 0;
    rstStatus = Gwac_geomap(matchpeervec, order, iter, rejsigma, xrms, yrms, xrms2, yrms2, outfilename, statusstr);
    printf("xrms=%f\n", xrms);
    printf("yrms=%f\n", yrms);
    printf("xrms2=%f\n", xrms2);
    printf("yrms2=%f\n", yrms2);
    printf("%s\n", statusstr);

    float xshift = 0, yshift = 0;
    rstStatus = GetShift(outfilename, xshift, yshift, statusstr);
    printf("xshift=%f\nyshift=%f\n", xshift, yshift);
    printf("%s\n", statusstr);
    
    free(points);
}

void testGeoxytran() {

    const char *dataName = "20130111_177d550752_60d851283-420.fit.sex";
    const char cofName[MAX_LINE_LENGTH] = "gwac_geomap_result.txt";
    const char outName[MAX_LINE_LENGTH] = "fits_result.txt";
    char statusstr[MAX_LINE_LENGTH];
    memset(statusstr, 0, MAX_LINE_LENGTH);

    DataStruct *points;
    int pointNum;
    readPoints(dataName, &points, &pointNum);

    vector<ST_STAR> stars;
    int i;
    for (i = 0; i < pointNum; i++) {
        ST_STAR star;
        star.x = points[i].xref;
        star.y = points[i].yref;
        stars.push_back(star);
    }

    int flag = -1;
    Gwac_geoxytran(stars, cofName, flag, statusstr);
    
    for (i = 0; i < pointNum; i++) {
        points[i].xin = stars.at(i).x;
        points[i].yin = stars.at(i).y;
    }
    
    writePoints(outName, points, pointNum);
}

void test2To5(){
    
    const char order2[MAX_LINE_LENGTH] = "fits_result_2order.txt";
    const char order5[MAX_LINE_LENGTH] = "fits_result_5order.txt";
    const char outfile[MAX_LINE_LENGTH] = "residual_2_5.txt";
    char statusstr[MAX_LINE_LENGTH];
    memset(statusstr, 0, MAX_LINE_LENGTH);

    DataStruct *points2, *points5;
    int pointNum2, pointNum5;
    readPoints(order2, &points2, &pointNum2);
    readPoints(order5, &points5, &pointNum5);
    
    float *residualsx = (float*)malloc(pointNum2*sizeof(float));
    float *residualsy = (float*)malloc(pointNum2*sizeof(float));
    
    FILE *fp = fopen(outfile, "w");
    
    int i;
    for(i=0; i<pointNum2; i++){
        residualsx[i] = points2[i].xin - points5[i].xin;
        residualsy[i] = points2[i].yin - points5[i].yin;
        float rms = sqrt(residualsx[i]*residualsx[i]+residualsy[i]*residualsy[i]);
        fprintf(fp, "%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", 
                points2[i].xin, points2[i].yin, 
                points5[i].xin, points5[i].yin, 
                residualsx[i], residualsy[i], rms);
    }
    
    fclose(fp);
    free(residualsx);
    free(residualsy);
}