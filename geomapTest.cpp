
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

void testGeomap(char *pairCatalog, char *mapParm, unsigned int order, unsigned int iter, float rejsigma) {

  printf("start geomap\n");
//  unsigned int order = 5;
//  unsigned int iter = 4;
//  float rejsigma = 2.5;
  float xrms, yrms;
  float xrms2, yrms2;
//  const char outfilename[MAX_LINE_LENGTH] = "test/map-config.txt";
  char statusstr[MAX_LINE_LENGTH];
  memset(statusstr, 0, MAX_LINE_LENGTH);

//  char *fName = "test/M1_01_151101_1_010020_0019_Mtest.dat"; //data  test.txt
  DataStruct *points;
  int pointNum;
  
  readPoints(pairCatalog, &points, &pointNum);
  printf("read %d pos pairs\n", pointNum);

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

  printf("run Gwac_geomap\n");
  int rstStatus = 0;
  rstStatus = Gwac_geomap(matchpeervec, order, iter, rejsigma, xrms, yrms, xrms2, yrms2, mapParm, statusstr);
  printf("xrms=%f\n", xrms);
  printf("yrms=%f\n", yrms);
  printf("xrms2=%f\n", xrms2);
  printf("yrms2=%f\n", yrms2);
  printf("%s\n", statusstr);

  float xshift = 0, yshift = 0;
  rstStatus = GetShift(mapParm, xshift, yshift, statusstr);
  printf("xshift=%f\nyshift=%f\n", xshift, yshift);
  printf("%s\n", statusstr);

  printf("geomap done\n");
  
  free(points);
}

void testGeoxytran(char *dataCatalog, char *mapParm, char *outName, int direction) {
  //M1_01_151101_1_010020_0019_Mtest.dat M1_01_151101_1_010020_0019.fit.sexini
//  const char *dataName = "test/M1_01_151101_1_010020_0019.fit.sexini"; //20130111_177d550752_60d851283-420.fit.sex geomap/data
  //    const char *dataName = "test/M1_01_151101_1_010020_0099.fit.tran3"; //20130111_177d550752_60d851283-420.fit.sex geomap/data
//  const char cofName[MAX_LINE_LENGTH] = "test/map-config.txt";
//  const char outName[MAX_LINE_LENGTH] = "test/fits_result.txt";
  char statusstr[MAX_LINE_LENGTH];
  memset(statusstr, 0, MAX_LINE_LENGTH);

  DataStruct *points;
  int pointNum;
  readPoints(dataCatalog, &points, &pointNum);

  vector<ST_STAR> stars;
  int i;
  for (i = 0; i < pointNum; i++) {
    ST_STAR star;
    star.x = points[i].xin;
    star.y = points[i].yin;
    stars.push_back(star);
  }

//  int flag = -1;
  Gwac_geoxytran(stars, mapParm, direction, statusstr);

  for (i = 0; i < pointNum; i++) {
    points[i].xin = stars.at(i).x;
    points[i].yin = stars.at(i).y;
  }

  writePoints(outName, points, pointNum);
}

void test2To5() {

  const char order2[MAX_LINE_LENGTH] = "test/M1_01_151101_1_010020_0099.fit.tran3";
  //    const char order2[MAX_LINE_LENGTH] = "test/M1_01_151101_1_010020_0019_Mtest.dat";
  const char order5[MAX_LINE_LENGTH] = "test/fits_result.txt";
  const char outfile[MAX_LINE_LENGTH] = "test/residual.txt";
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
    residualsx[i] = points2[i].xref - points5[i].xin;
    residualsy[i] = points2[i].yref - points5[i].yin;
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

  //    xi = -0.005405;
  //    eta = -0.078667;


  tanPlaneToSphere(ra0, dec0, xi, eta, ra, dec);
  printf("ra=%f dec=%f\n", ra, dec);
  tanPlaneToSphere2(ra0, dec0, xi, eta, ra, dec);
  printf("ra=%f dec=%f\n", ra, dec);

}

void testProjection2() {
  //1039.780        1599.120    119.61845398        67.31501007
  double ra0 = 120.4146562756105;
  double dec0 = 71.85829263128276;
  double ra = 119.61845398;
  double dec = 67.31501007;
  double x = 1039.780;
  double y = 1599.120;
  double xi;
  double eta;

  double a = 1402.868558469206;
  double b = -972.282292734403;
  double c = -29871.65438087219;
  double d = -2550.406313606504;
  double e = 29880.52202304072;
  double f = -973.3084798633404;

  printf("x=%f y=%f\n", x, y);
  x = (x - 1500) / 1000;
  y = (y - 1500) / 1000;
  printf("x=%f y=%f\n", x, y);
  xi = a + b * x + c*y;
  eta = d + e * x + f*y;
  printf("xi=%f eta=%f\n", xi, eta);

  xi /= SECOND_TO_RADIANS;
  eta /= SECOND_TO_RADIANS;
  double ra2 = 0, dec2 = 0;
  tanPlaneToSphere(ra0, dec0, xi, eta, ra2, dec2);
  printf("ra=%f dec=%f\n", ra2, dec2);

  int rst = tanSphereToPlane(ra0, dec0, ra, dec, xi, eta);
  printf("xi=%f eta=%f\n", xi, eta);
  printf("rst=%d\n", rst);
  xi *= SECOND_TO_RADIANS;
  eta *= SECOND_TO_RADIANS;
  printf("xi=%f eta=%f\n", xi, eta);

  y = (e * xi - b * eta + b * d - a * e) / (e * c - b * f);
  x = (xi - c * y - a) / b;
  double x2 = (eta - f * y - d) / e;
  printf("x=%f\n", x);
  printf("x=%f\n", x2);
  printf("y=%f\n", y);
  x = x * 1000 + 1500;
  y = y * 1000 + 1500;
  printf("x=%f y=%f\n", x, y);

  //    xi = -0.005405;
  //    eta = -0.078667;


  tanPlaneToSphere(ra0, dec0, xi, eta, ra, dec);
  printf("ra=%f dec=%f\n", ra, dec);
  tanPlaneToSphere2(ra0, dec0, xi, eta, ra, dec);
  printf("ra=%f dec=%f\n", ra, dec);

}

void testCctran() {
  printf("start");

      const char *dataName = "testdata/data11.txt";
      const char cofName[MAX_LINE_LENGTH] = "testdata/refcom.acc";
      const char outName1[MAX_LINE_LENGTH] = "testdata/xy2radec.txt";
      const char outName2[MAX_LINE_LENGTH] = "testdata/radec2xy.txt";
//  const char *dataName = "cctran/a.axycc";
//  const char cofName[MAX_LINE_LENGTH] = "cctran/refcom.acc";
//  const char outName1[MAX_LINE_LENGTH] = "cctran/xytoradec.txt";
//  const char outName2[MAX_LINE_LENGTH] = "cctran/radectoxy.txt";
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
      star.ra = points[i].xin;
      star.dec = points[i].yin;
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
      star.x = points[i].xref;
      star.y = points[i].yref;
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
  printf("done");
}
