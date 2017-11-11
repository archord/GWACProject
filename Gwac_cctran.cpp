

#include "gwac.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mathFunction.h"
#include "WCSTNX.h"

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

  WCSTNX tnx;
  if (!tnx.LoadText(transfilename)) {
    sprintf(statusstr, "Error Code: %d\n"\
                "File %s line %d, open file \"%s\" error!\n",
            GWAC_OPEN_FILE_ERROR, __FILE__, __LINE__, transfilename);
    return GWAC_OPEN_FILE_ERROR;
  }

  if (flag == 1) {
    for (int i = 0; i < objvec.size(); i++) {
      ST_STAR &star = objvec.at(i);
      double ra = 0.0, dec = 0.0;
      tnx.XY2WCS(star.x, star.y, ra, dec);
      star.ra = ra;
      star.dec = dec;
    }
  } else if (flag == -1) {
    for (int i = 0; i < objvec.size(); i++) {
      ST_STAR &star = objvec.at(i);
      double x = 0.0, y = 0.0;
      tnx.WCS2XY(star.ra, star.dec, x, y);
      star.x = x;
      star.y = y;
    }
  }

  return GWAC_SUCCESS;
}
