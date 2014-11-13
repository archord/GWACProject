
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

#include "mathFunction.h"

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
    const char outfilename[MAX_LINE_LENGTH] = "cctran/map.txt";
    char statusstr[MAX_LINE_LENGTH];
    memset(statusstr, 0, MAX_LINE_LENGTH);

    char *fName = "cctran/data"; //data  test.txt
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
    const char outName[MAX_LINE_LENGTH] = "fits_result2.txt";
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

void test2To5() {

    const char order2[MAX_LINE_LENGTH] = "fits_result_2order.txt";
    const char order5[MAX_LINE_LENGTH] = "fits_result_5order.txt";
    const char outfile[MAX_LINE_LENGTH] = "residual_2_5.txt";
    char statusstr[MAX_LINE_LENGTH];
    memset(statusstr, 0, MAX_LINE_LENGTH);

    DataStruct *points2, *points5;
    int pointNum2, pointNum5;
    readPoints(order2, &points2, &pointNum2);
    readPoints(order5, &points5, &pointNum5);

    float *residualsx = (float*) malloc(pointNum2 * sizeof (float));
    float *residualsy = (float*) malloc(pointNum2 * sizeof (float));

    FILE *fp = fopen(outfile, "w");

    int i;
    for (i = 0; i < pointNum2; i++) {
        residualsx[i] = points2[i].xin - points5[i].xin;
        residualsy[i] = points2[i].yin - points5[i].yin;
        float rms = sqrt(residualsx[i] * residualsx[i] + residualsy[i] * residualsy[i]);
        fprintf(fp, "%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
                points2[i].xin, points2[i].yin,
                points5[i].xin, points5[i].yin,
                residualsx[i], residualsy[i], rms);
    }

    fclose(fp);
    free(residualsx);
    free(residualsy);
}

int testMacro1(char *statusstr) {
    GWAC_REPORT_ERROR(GWAC_ERROR, "error test1");
}

int testMacro2(char *statusstr) {
    CHECK_RETURN_SATUS(GWAC_OPEN_FILE_ERROR);
}

int testMacro3(char *statusstr) {
    CHECK_STATUS_STR_IS_NULL(NULL);
}

int testMacro4(char *statusstr) {
    CHECK_INPUT_IS_NULL(NULL, "var test");
}

int testMacro5(char *statusstr) {
    MALLOC_IS_NULL();
}

int testMacro6(char *statusstr) {
    CHECK_OPEN_FILE(NULL, "file test");
}

void testMacro() {
    char statusstr[255] = "";

    int rst = testMacro1(statusstr);
    printf("1######### %s", statusstr);
    rst = testMacro2(statusstr);
    printf("2######### %s", statusstr);
    rst = testMacro3(statusstr);
    printf("3######### %s", statusstr);
    rst = testMacro4(statusstr);
    printf("4######### %s", statusstr);
    rst = testMacro5(statusstr);
    printf("5######### %s", statusstr);
    rst = testMacro6(statusstr);
    printf("6######### %s", statusstr);
}

void testProjection() {
    double ra0 = 120.4;
    double dec0 = 71.8;
    double ra = 119.6;
    double dec = 67.3;
    double xi;
    double eta;

    int rst = tanSphereToPlane(ra0, dec0, ra, dec, xi, eta);
    printf("xi=%f eta=%f\n", xi, eta);
    printf("rst=%d\n", rst);
    tanSphereToPlane2(ra0, dec0, ra, dec, xi, eta);
    printf("xi=%f eta=%f\n", xi, eta);

    xi = -0.005405;
    eta = -0.078667;


    tanPlaneToSphere(ra0, dec0, xi, eta, ra, dec);
    printf("ra=%f dec=%f\n", ra, dec);
    tanPlaneToSphere2(ra0, dec0, xi, eta, ra, dec);
    printf("ra=%f dec=%f\n", ra, dec);

}

void testCctran() {

    const char *dataName = "cctran/a.axycc";
    const char cofName[MAX_LINE_LENGTH] = "cctran/refcom.acc";
    const char outName1[MAX_LINE_LENGTH] = "cctran/xytoradec.txt";
    const char outName2[MAX_LINE_LENGTH] = "cctran/radectoxy.txt";
    char statusstr[MAX_LINE_LENGTH];
    memset(statusstr, 0, MAX_LINE_LENGTH);

    DataStruct *points;
    int pointNum;
    readPoints(dataName, &points, &pointNum);

    int flag = 1;
    if (flag == 1) {
        vector<ST_STAR> stars;
        int i;
        for (i = 0; i < pointNum; i++) {
            ST_STAR star;
            star.x = points[i].xref;
            star.y = points[i].yref;
            stars.push_back(star);
        }
        Gwac_cctran(stars, cofName, flag, statusstr);
        printf("%s\n", statusstr);

        for (i = 0; i < pointNum; i++) {
            points[i].xref = stars.at(i).x;
            points[i].yref = stars.at(i).y;
            points[i].xin = stars.at(i).ra;
            points[i].yin = stars.at(i).dec;
        }
        writePoints(outName1, points, pointNum);
    } else {
        vector<ST_STAR> stars;
        int i;
        for (i = 0; i < pointNum; i++) {
            ST_STAR star;
            star.ra = points[i].xin;
            star.dec = points[i].yin;
            stars.push_back(star);
        }
        Gwac_cctran(stars, cofName, flag, statusstr);
        printf("%s\n", statusstr);

        for (i = 0; i < pointNum; i++) {
            points[i].xref = stars.at(i).x;
            points[i].yref = stars.at(i).y;
        }
        writePoints(outName2, points, pointNum);
    }
}
