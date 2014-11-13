

#include "gwac.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*******************************************************************************
 * 功能：星表坐标转换，平面坐标系到平面坐标系的转换，使用多项式拟合方式进行转换
 *      根据指定的参数文件中的多项式拟合系数，对给定的星表列表进行坐标转换
 * 
 **输入：
 *      objvec星表列表，输入时为待转换坐标，输出时为转换后坐标
 *      transfilename 拟合参数文件名，其格式为geomap软件输出格式
 *      flag 控制转换方向
 *          当flag= 1时，反向拟合，objvec输入为（xref,yref），输出为（xin,yin），读
 *          取参赛文件中第二组拟合系数中surface2后的系数
 *          当flag=-1时，正向拟合，objvec输入为（xin,yin），输出为（xref,yref），读
 *          取参赛文件中第一组拟合系数中surface2后的系数
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
    CHECK_INPUT_IS_NULL(transfilename, "transfilename");

    /*检测星表数据数组是否为空*/
    if (objvec.empty()) {
        sprintf(statusstr, "Error Code: %d\n"
                "In Gwac_geoxytran, the input parameter objvec is empty!\n",
                GWAC_FUNCTION_INPUT_EMPTY);
        return GWAC_FUNCTION_INPUT_EMPTY;
    }

    FILE *fp = fopen(transfilename, "r");
    CHECK_OPEN_FILE(fp, transfilename);

    int surface2 = 0;
    int cofNum = 0, cofIdx = 0;

    /*startFlag=1时，控制读取第二组参数，第二个surface2后面的参数*/
    int startFlag = 0;
    int reverseFlag = 0;
    char line[MAX_LINE_LENGTH];

    /*读取拟合参数的行数，也就是surface2后面的数据*/
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) {
        if (strstr(line, "surface2") != NULL) {
            /*如果flag==-1，则进行正向拟合，读取第一组多项式系数*/
            /*如果flag==1，则进行反向拟合，读取第二组多项式系数*/
            if ((flag == -1) || (flag == 1 && reverseFlag == 1)) {
                sscanf(line, "%*s%d", &surface2);
                startFlag = 1;
                break;
            }
            reverseFlag = 1;
        }
    }

    /*前8行为其他信息*/
    cofNum = surface2 - 8;
    double *xcof = (double*) malloc(cofNum * sizeof (double));
    double *ycof = (double*) malloc(cofNum * sizeof (double));
    double *afunc = (double*) malloc((cofNum + 1) * sizeof (double));
    if (xcof == NULL || ycof == NULL || afunc == NULL) {
        if (xcof != NULL) free(xcof);
        if (ycof != NULL) free(ycof);
        if (afunc != NULL) free(afunc);
        MALLOC_IS_NULL();
    }

    int cofStartLine = 0;
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) {
        if (startFlag == 1) {
            cofStartLine++;
            /*从surface2开始的第9行开始读取拟合参数*/
            if (cofStartLine >= 9) {
                sscanf(line, "%lf %lf", &xcof[cofIdx], &ycof[cofIdx]);
                cofIdx++;
            }
        }
        /*读完一组多项式系数，则跳出循环*/
        if (cofStartLine >= surface2) {
            break;
        }
    }
    fclose(fp);

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

    if (xcof != NULL) free(xcof);
    if (ycof != NULL) free(ycof);
    if (afunc != NULL) free(afunc);

    return GWAC_SUCCESS;
}
