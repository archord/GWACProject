

#include "gwac.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*******************************************************************************
 * 功能：星表坐标转换，根据指定的参数文件中的多项式拟合系数，对给定的星表列表进行坐标转换
 * 
 **输入：
 *      objvec星表列表，输入时为待转换坐标，输出时为转换后坐标
 *      transfilename 拟合参数文件名，其格式为geomap软件输出格式
 *      flag 控制转换方向
 *          当flag= 1时，objvec输入为（xref,yref），输出为（xin,yin），直接读取参赛文件
 *          当flag=-1时，objvec输入为（xin,yin），输出为（xref,yref），读取参数文件
 *          transfilename.reverse
 * 
 **输出
 *      objvec星表列表，输入时为待转换坐标，输出时为转换后坐标
 *          
 **返回值:
 *      0表示正确，其它值为错误码，
 *      蔡使用2001～2999，徐使用3001～3999，苑使用4001～4999，李使用5001～5999
 ******************************************************************************/
int Gwac_geoxytran(vector<ST_STAR> &objvec,
        const char transfilename[],
        int flag,
        char statusstr[]) {

    /*检查错误结果输出参数statusstr是否为空*/
    CHECK_STATUS_STR_IS_NULL(statusstr);
    /*检测输入参数是否为空*/
    CHECK_INPUT_IS_NULL("Gwac_geoxytran", transfilename, "transfilename");

    /*检测星表数据数组是否为空*/
    if (objvec.empty()) {
        sprintf(statusstr, "Error Code: %d\n"
                "In Gwac_geoxytran, the input parameter objvec is empty!\n",
                GWAC_FUNCTION_INPUT_EMPTY);
        return GWAC_FUNCTION_INPUT_EMPTY;
    }

    char cofFileName[MAX_LINE_LENGTH];
    if (flag == 1) {
        sprintf(cofFileName, "%s", transfilename);
    } else {
        sprintf(cofFileName, "%s.reverse", transfilename);
    }

    FILE *fp = fopen(cofFileName, "r");
    CHECK_OPEN_FILE("Gwac_geoxytran", fp, transfilename);

    int surface2 = 0;
    int cofNum = 0, cofIdx = 0;
    double *xcof, *ycof;

    int startFlag = 0;
    int cofStartLine = 0;
    char line[MAX_LINE_LENGTH];
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) {
        if (strstr(line, "surface2") != NULL) {
            sscanf(line, "%*s%d", &surface2);
            /*前8行为其他信息*/
            cofNum = surface2 - 8;
            xcof = (double*) malloc(cofNum * sizeof (double));
            CHECK_MALLOC_IS_NULL("Gwac_geoxytran", xcof, "xcof");
            ycof = (double*) malloc(cofNum * sizeof (double));
            CHECK_MALLOC_IS_NULL("Gwac_geoxytran", ycof, "ycof");
            startFlag = 1;
        }
        if (startFlag == 1) {
            cofStartLine++;
            /*从surface2开始的第9行开始读取拟合参数*/
            if (cofStartLine > 9) {
                sscanf(line, "%lf %lf", &xcof[cofIdx], &ycof[cofIdx]);
                cofIdx++;
            }
        }
    }
    fclose(fp);
    
    double *afunc = (double*) malloc((cofNum + 1) * sizeof (double));
    CHECK_MALLOC_IS_NULL("Gwac_geoxytran", afunc, "afunc");

    int i, j;
    for (i = 0; i < objvec.size(); i++) {
        ST_STAR &star = objvec.at(i);
        cofun(star.x, star.y, afunc, cofNum);
        double tmpx = 0.0, tmpy = 0.0;
        for (j = 1; j <= cofNum; j++) {
            tmpx += xcof[j - 1] * afunc[j];
            tmpy += ycof[j - 1] * afunc[j];
        }
        star.x = tmpx;
        star.y = tmpy;
    }

    if (startFlag == 1) {
        free(xcof);
        free(ycof);
    }
    free(afunc);

    return GWAC_SUCCESS;
}
