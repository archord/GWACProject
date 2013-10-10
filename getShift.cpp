
#include "gwac.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*******************************************************************************
 * 功能：从指定文件获取线性位移项
 * 
 **输入：
 *      transfilename 拟合参数文件名，其格式为geomap软件输出格式
 * 
 **输出
 *      xshift x方向线性位移
 *      yshift y方向线性位移
 *          
 **返回值:
 *      0表示正确，其它值为错误码，
 *      蔡使用2001～2999，徐使用3001～3999，苑使用4001～4999，李使用5001～5999
 ******************************************************************************/
int GetShift(const char transfilename[], 
        float &xshift, 
        float &yshift,
        char statusstr[]) {

    /*检查错误结果输出参数statusstr是否为空*/
    CHECK_STATUS_STR_IS_NULL(statusstr);
    /*检测输入参数是否为空*/
    CHECK_INPUT_IS_NULL("GetShift", transfilename, "transfilename");

    FILE *fp = fopen(transfilename, "r");
    CHECK_OPEN_FILE("GetShift", fp, transfilename);

    int flag = 0;
    char line[MAX_LINE_LENGTH];
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL && flag < 2) {
        if (strstr(line, "xshift") != NULL) {
            sscanf(line, "%*s %f", &xshift);
            flag++;
        } else if (strstr(line, "yshift") != NULL) {
            sscanf(line, "%*s %f", &yshift);
            flag++;
        }
    }
    fclose(fp);

    return GWAC_SUCCESS;
}
