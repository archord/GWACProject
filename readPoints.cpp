
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gwac.h"
#include "readPoints.h"

int readPoints(char *fName, DataStruct **points, int *pointNum) {

    if (fName == NULL)
        return GWAC_ERROR;
    
    FILE *fp = fopen(fName, "r");
    if (fp == NULL)
        return GWAC_OPEN_FILE_ERROR;

    char *line = (char*) malloc(MAX_LINE_LENGTH);
    int lineNum = 0;
    //count file line number
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) {
        lineNum++;
    }
    *pointNum = lineNum;
    fseek(fp, 0, SEEK_SET);

    DataStruct *tPoints = (DataStruct *) malloc(lineNum * sizeof (DataStruct));
    
    if(tPoints == NULL)
        return GWAC_MALLOC_ERROR;

    int i;
    float xref, yref, xin, yin;
    for (i = 0; i < lineNum; i++) {
        fgets(line, MAX_LINE_LENGTH, fp);
        sscanf(line, "%f %f %f %f",
                &xref, &yref, &xin, &yin);
        tPoints[i].xref = xref;
        tPoints[i].yref = yref;
        tPoints[i].xin = xin;
        tPoints[i].yin = yin;
    }
    
    *points = tPoints;
    free(line);
    fclose(fp);

    return GWAC_SUCCESS;
}
