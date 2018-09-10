/* 
 * File:   main.cpp
 * Author: xy
 *
 * Created on September 9, 2013, 9:07 AM
 */

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include "gwac.h"
#include "testFunction.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

  if (argc > 2) {
    char *method = argv[1];
    if (strcmp(method, "geomap")==0) {
      if (argc == 7) {
        //      testGeomap(pairCatalog, mapParm, order, iter, rejsigma);
        testGeomap(argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atof(argv[6]));
      } else {
        printf("exeprog geomap pairs.cat map_parm.txt order=5 iterNum=4 rejSigma=2.5\n");
      }
    } else if (strcmp(method, "geoxytran")==0) {

      if (argc == 6) {
        //      testGeoxytran(char *dataCatalog, char *mapParm, char *outName, int direction);
        testGeoxytran(argv[2], argv[3], argv[4], atoi(argv[5]));
      } else {
        printf("exeprog geoxytran data.cat mapParm.txt data.out direction=-1\n");
      }
    } else {
      printf("unknown method %s\n", method);
    }
  } else {
    printf("exeprog geomap pairs.cat map_parm.txt order=5 iterNum=4 rejSigma=2.5\n");
    printf("exeprog geoxytran data.cat mapParm.txt data.out direction=-1\n");
  }

  return 0;
}

