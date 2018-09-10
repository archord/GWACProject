/*
 * @file WCSTNX.h 声明文件, 基于非标准WCS格式TNX, 计算图像坐标与WCS坐标之间的对应关系
 * @version 0.1
 * @date 2017年11月9日
 * - 从FITS文件头加载WCS TNX参数项
 * - 从文本文件加载WCS TNX参数项
 * - 计算(x,y)对应的WCS坐标(ra, dec)
 *
 * @note
 * 畸变改正参数项数量(order>=1， 二元)
 * none:   2*order-1
 * half:   order*(order+1)/2
 * full:   order*order
 *
 * @note
 * 畸变改正项排列顺序(order>=1, 二元)
 * X0Y0
 * X1Y0
 * ...
 * X(order-1)Y0
 * X1Y1
 * X2Y1
 * ...
 * X(order-1)Y1
 * ...
 * X0Y(order-1)
 */

#ifndef WCSTNX_H_
#define WCSTNX_H_

#include <string>
#include "ADefine.h"

using std::string;
using AstroUtil::PT2F;

class WCSTNX {
public:
	WCSTNX();
	virtual ~WCSTNX();

protected:
	/* 数据结构 */
	enum {// 平面畸变拟合多项式类型
		TNX_CHEB = 1,	//< 契比雪夫多项式
		TNX_LEG,			//< 勒让德多项式
		TNX_POLY			//< 线性多项式
	};

	enum {// 多项式交叉系数类型
		TNX_XNONE,		//< 无交叉项
		TNX_XFULL,		//< 全交叉
		TNX_XHALF		//< 半交叉
	};

	struct param_surface {// 畸变修正系数
		int xsurface, ysurface;		//< 多项式类型
		int xxorder, xyorder;		//< 阶次
		int yxorder, yyorder;		//< 阶次
		int xxterm, yxterm;			//< 交叉项类型
		double xxmin, xxmax, xymin, xymax;	//< x范围
		double yxmin, yxmax, yymin, yymax;	//< y范围
		int xncoef, yncoef;	//< 系数数量
		double *xcoef, *ycoef;		//< 系数
		double *xx, *xy, *xxy;		//< x多项式单项变量
		double *yx, *yy, *yxy;		//< y多项式单项变量

	public:
		param_surface() {
			xsurface = ysurface = -1;
			xxorder = xyorder = -1;
			yxorder = yyorder = -1;
			xxterm = yxterm = -1;
			xxmin = xxmax = xymin = xymax = 0.0;
			yxmin = yxmax = yymin = yymax = 0.0;
			xncoef = 0, yncoef = 0;
			xcoef = ycoef = NULL;
			xx = xy = xxy = NULL;
			yx = yy = yxy = NULL;
		}

		void free_array(double **array) {
			if ((*array) != NULL) {
				delete [](*array);
				(*array) = NULL;
			}
		}

		virtual ~param_surface() {
			free_array(&xcoef);
			free_array(&ycoef);
			free_array(&xx);
			free_array(&xy);
			free_array(&xxy);
			free_array(&yx);
			free_array(&yy);
			free_array(&yxy);
		}

		int item_count(int type, int xorder, int yorder) {
			int order = xorder < yorder ? xorder : yorder;
			int n;

			if      (type == TNX_XNONE) n = xorder + yorder - 1;
			else if (type == TNX_XFULL) n = xorder * yorder;
			else if (type == TNX_XHALF) 	n = xorder * yorder - order * (order - 1) / 2;

			return n;
		}

		void set_orderx(int x, int y) {
			if (xxorder != x) {
				xxorder = x;
				free_array(&xx);
			}
			if (xyorder != y) {
				xyorder = y;
				free_array(&xy);
			}

			if (!xx) xx = new double[x];
			if (!xy) xy = new double[y];
		}

		void set_ordery(int x, int y) {
			if (yxorder != x) {
				yxorder = x;
				free_array(&yx);
			}
			if (yyorder != y) {
				yyorder = y;
				free_array(&yy);
			}

			if (!yx) yx = new double[x];
			if (!yy) yy = new double[y];
		}

		void set_xtermx(int xterm) {
			int n = item_count(xterm, xxorder, xyorder);
			if (xxterm != xterm) xxterm = xterm;
			if (n != xncoef) {
				free_array(&xxy);
				free_array(&xcoef);
				xncoef = n;
			}
			if (!xxy) xxy = new double[n];
			if (!xcoef) xcoef = new double[n];
		}

		void set_xtermy(int xterm) {
			int n = item_count(xterm, yxorder, yyorder);
			if (yxterm != xterm) yxterm = xterm;
			if (n != yncoef) {
				free_array(&yxy);
				free_array(&ycoef);
				yncoef = n;
			}
			if (!yxy) yxy = new double[n];
			if (!ycoef) ycoef = new double[n];
		}
	};

	struct param_tnx {// TNX参数
		bool valid1, valid2;	//< 参数有效性标志
		PT2F ref_xymean;		//< 参考点: 平均XY坐标
		PT2F ref_wcsmean;	//< 参考点: 平均WCS坐标, 量纲: 弧度
		string pixsystem;	//< 图像像素坐标系名称
		string coosystem;	//< WCS坐标系名称
		string projection;	//< 投影模式
		PT2F ref_xy;			//< 参考点: XY坐标
		PT2F ref_wcs;		//< 参考点: WCS坐标, 量纲: 弧度
		string function;		//< 畸变改正函数类型
		PT2F shift;			//< 投影坐标偏移量
		PT2F mag;			//< 比例尺, 量纲: 弧度/像素
		PT2F rotation;		//< 旋转角, 量纲: 弧度
//		double cd[2][2];		//< 旋转矩阵. 由平均关系获得. 量纲: 弧度/像素
		param_surface surface1;	//< 一阶投影面拟合系数
		param_surface surface2;	//< 残差投影面拟合系数
	};

protected:
	/* 成员变量 */
	param_tnx param_;	//< TNX参数

public:
	/* 接口 */
	/*!
	 * @brief 从FITS文件头加载WCS参数
	 * @param filepath FITS文件路径
	 * @return
	 * 参数加载结果
	 */
	bool LoadImage(const char* filepath);
	/*!
	 * @brief 从文本文件加载WCS参数
	 * @param filepath 文本文件路径
	 * @return
	 * 参数加载结果
	 */
	bool LoadText(const char* filepath);
	/*!
	 * @brief 计算与图像坐标(x,y)对应的WCS坐标(ra,dec)
	 * @param x   X轴坐标
	 * @param y   Y轴坐标
	 * @param ra  赤经, 量纲: 弧度
	 * @param dec 赤纬, 量纲: 弧度
	 * @return
	 * 计算结果
	 *  0: 正确
	 * -1: 错误(未加载参数项)
	 */
	int XY2WCS(double x, double y, double& ra, double& dec);
        int WCS2XY(double ra, double dec, double& x, double& y);

protected:
	/* 功能 */
	/*!
	 * @brief 计算一元多阶线性数组
	 * @param x      自变量
	 * @param order  阶次
	 * @param ptr    输出数组
	 * @note
	 * 函数创建数组存储空间
	 */
	void linear_array(double x, int order, double* ptr);
	/*!
	 * @brief 计算一元多阶勒让德数组
	 * @param x      自变量
	 * @param xmin   自变量有效范围最小值
	 * @param xmax   自变量有效范围最大值
	 * @param order  阶次
	 * @param ptr    输出数组
	 * @note
	 * 函数创建数组存储空间
	 */
	void legendre_array(double x, double xmin, double xmax, int order, double* ptr);
	/*!
	 * @brief 计算一元多阶契比雪夫数组
	 * @param x      自变量
	 * @param xmin   自变量有效范围最小值
	 * @param xmax   自变量有效范围最大值
	 * @param order  阶次
	 * @param ptr    输出数组
	 * @note
	 * 函数创建数组存储空间
	 */
	void chebyshev_array(double x, double xmin, double xmax, int order, double* ptr);
	/*!
	 * @brief 计算二元多项式变量数组
	 * @param type   交叉项类型
	 * @param xorder x阶次
	 * @param yorder y阶次
	 * @param x      第一个变量数组, 其长度为order
	 * @param y      第二个变量数组, 其长度为order
	 * @param n      输出数组长度
	 * @param array  输出数组, 其长度与order和type有关
	 */
	void polyval_item(int type, int xorder, int yorder, double* x, double* y, int n, double* array);
	/*!
	 * @brief 计算多项式和
	 * @param n    数组长度
	 * @param coef 系数
	 * @param item 单项数值
	 * @return
	 * 多项式和
	 */
	double polysum(int n, double* coef, double* item);
	/*!
	 * @brief 计算残差修正项
	 * @param x    输入xi坐标
	 * @param y   输入eta坐标
	 * @param dx   xi修正量
	 * @param dy  eta修正量
	 */
	void correct(param_surface& surface, double x, double y, double& dx, double& dy);
	/*!
	 * @brief 图像坐标转换为投影平面平坐标
	 * @param x    图像X坐标
	 * @param y    图像Y坐标
	 * @param xi   投影xi坐标
	 * @param eta  投影eta坐标
	 */
	void image_to_plane(double x, double y, double& xi, double& eta);
	/*!
	 * @brief TAN投影逆变换, 将平面坐标转换为球面坐标
	 * @param x    平面坐标1
	 * @param y    平面坐标2
	 * @param ra   赤经, 量纲: 弧度
	 * @param dec  赤纬, 量纲: 弧度
	 */
	void plane_to_wcs(double x, double y, double& ra, double& dec);
        int wcs_to_plane(double ra, double dec, double &xi, double &eta);
};

#endif /* WCSTNX_H_ */
