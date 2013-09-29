
#include <iostream>
#include <vector>
#include <string>

using namespace std;

#include <cstdlib>
#include <math.h>
#include "gwac.h"
#include "readPoints.h"

float testFunction1(float x, float y){
    float rst = 3.0 + 20.0*x + 2.0*x*x + 1.4*y + 0.06*x*y + 0.8*y*y;
    int tmp = rand()%20;
    if(tmp>10){
        tmp -= 10;
        rst = rst*(1.0 + 1.0*tmp/1000);
    }else{
        rst = rst*(1.0 - 1.0*tmp/1000);
    }
    return rst;
}

float testFunction2(float x, float y){
    float rst = 2.0 + 30.5*x + 8.0*x*x + 0.008*y + 0.8*x*y + 10*y*y;
    int tmp = rand()%40;
    if(tmp>20){
        tmp -= 20;
        rst = rst*(1.0 + 1.0*tmp/1000);
    }else{
        rst = rst*(1.0 - 1.0*tmp/1000);
    }
    return rst;
}

int testGeomap() {
    
    unsigned int order = 5;
    unsigned int iter = 3;
    float rejsigma = 2.5;
    float xrms, yrms;
    const char outfilename[MAX_LINE_LENGTH] = "gwac_geomap_result.txt";
    char statusstr[MAX_LINE_LENGTH];

    char *fName = "data";
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
    
    Gwac_geomap(matchpeervec, order, iter, rejsigma, xrms, yrms, outfilename, statusstr) ;
}
