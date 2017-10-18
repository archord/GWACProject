

#include "gwac.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mathFunction.h"

/*******************************************************************************
 * 功能：星表坐标转换，平面坐标系与天球坐标系的转换，平面坐标系与天球坐标系之间有一个标准
 *      坐标系作为过渡，平面坐标系与标准坐标系之间使用多项式拟合方式进行转换，标准坐标系
 *      与天球坐标系之间使用切平面投影方式进行转换。
 *      根据指定的参数文件中的多项式拟合系数，对给定的星表列表进行坐标转换
 * 
 **输入：
 *      objvec星表列表，输入时为待转换坐标，输出时为转换后坐标
 *      transfilename 拟合参数文件名，其格式为ccmap软件输出格式
 *      flag 控制转换方向
 *          当flag= 1时，正向拟合，objvec输入为（x,y），输出为（ra,dec）
 *          当flag=-1时，反向拟合，objvec输入为（ra,dec），输出为（x,y），读
 *          取参赛文件中surface2后的系数
 * 
 **输出
 *      objvec星表列表，输入时为待转换坐标，输出时为转换后坐标
 *          
 **返回值:
 *      0表示正确，其它值为错误码，
 *      蔡使用2001～2999，徐使用3001～3999，苑使用4001～4999，李使用5001～5999
 ******************************************************************************/
int Gwac_cctran(vector<ST_STAR> &objvec,
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

    int surface1 = 0;
    int surface2 = 0;
    int cofNum = 0;
    double lngref = 0.0, latref = 0.0;
    float minbnd = 0.0, maxbnd = 0.0;

    char line[MAX_LINE_LENGTH];

    /*读取拟合参数的行数，也就是surface2后面的数据*/
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) {
        if (strstr(line, "lngref") != NULL) {
            sscanf(line, "%*s%lf", &lngref);    //ra中值
        } else if (strstr(line, "latref") != NULL) {
            sscanf(line, "%*s%lf", &latref);    //dec中值
        } else if (strstr(line, "surface1") != NULL) {
            sscanf(line, "%*s%d", &surface1);  //线性参数项8+3
            break;
        }
    }
    printf("lngref=%lf\n", lngref);
    printf("latref=%lf\n", latref);

    double xcofl[3] = {0.0};
    double ycofl[3] = {0.0};
    int cofIdx = 0;
    int cofStartLine = 0;
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) {
        cofStartLine++;
        if(cofStartLine==5){ /*legendre拟合边界最小值*/
            sscanf(line, "%f %*s", &minbnd);
        }else if(cofStartLine==6){ /*legendre拟合边界最大值*/
            sscanf(line, "%f %*s", &maxbnd);
        }else if (cofStartLine >= 9) { /*从surface2开始的第9行开始读取拟合参数*/
            printf("%s\n", line);
            sscanf(line, "%lf %lf", &xcofl[cofIdx], &ycofl[cofIdx]);
            cofIdx++;
        }
        /*读完一组多项式系数，则跳出循环*/
        if (cofStartLine >= surface1) {
            break;
        }
    }
    
    fgets(line, MAX_LINE_LENGTH, fp);
    sscanf(line, "%*s%d", &surface2);  //高阶参数项

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

    cofIdx = 0;
    cofStartLine = 0;
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) {
        cofStartLine++;
        /*从surface2开始的第9行开始读取拟合参数*/
        if (cofStartLine >= 9) {
            printf("%s\n", line);
            sscanf(line, "%lf %lf", &xcof[cofIdx], &ycof[cofIdx]);
            cofIdx++;
        }
        /*读完一组多项式系数，则跳出循环*/
        if (cofStartLine >= surface2) {
            break;
        }
    }
    fclose(fp);

    float xcenter = (minbnd+maxbnd)/2;
    float hafbnd = (maxbnd - minbnd)/2;
    int i, j;
    if (flag == 1) {
        for (i = 0; i < objvec.size(); i++) {
            ST_STAR &star = objvec.at(i);
            
            double normalX = (star.x - xcenter)/hafbnd;
            double normalY = (star.y - xcenter)/hafbnd;
            double xi = xcofl[0] + xcofl[1]*normalX+ xcofl[2]*normalY;
            double eta = ycofl[0] + ycofl[1]*normalX + ycofl[2]*normalY;

            cofun_Legendre(normalX, normalY, afunc, cofNum);
            double xires = 0.0, etares = 0.0;
            for (j = 1; j <= cofNum; j++) {
                xires += xcof[j - 1] * afunc[j];
                etares += ycof[j - 1] * afunc[j];
            }
            
            xi = xi + xires;
            eta = eta + etares;
            xi/=SECOND_TO_RADIANS;
            eta/=SECOND_TO_RADIANS;
            
            double ra = 0.0, dec = 0.0;
            tanPlaneToSphere(lngref, latref, xi, eta, ra, dec);
            star.ra = ra;
            star.dec = dec;
        }
    } else if (flag == -1) {
        for (i = 0; i < objvec.size(); i++) {
            ST_STAR &star = objvec.at(i);
            
            double xi = 0.0, eta = 0.0;
            tanSphereToPlane(lngref, latref, star.ra, star.dec, xi, eta);
            xi*=SECOND_TO_RADIANS;
            eta*=SECOND_TO_RADIANS;
                        
            double x = xcofl[0] + xcofl[1]*xi + xcofl[2]*eta;
            double y = ycofl[0] + ycofl[1]*xi+ ycofl[2]*eta;
                        
            cofun_Legendre(xi, eta, afunc, cofNum);
            double xres = 0.0;
            double yres = 0.0;
            for (j = 1; j <= cofNum; j++) {
                xres += xcof[j - 1] * afunc[j];
                yres += ycof[j - 1] * afunc[j];
            }
            star.x = (x + xres)*hafbnd+xcenter;
            star.y = (y + yres)*hafbnd+xcenter;
        }
    }

    if (xcof != NULL) free(xcof);
    if (ycof != NULL) free(ycof);
    if (afunc != NULL) free(afunc);

    return GWAC_SUCCESS;
}
