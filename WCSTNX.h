/*
 * @file WCSTNX.h 声明文件, 基于非标准WCS格式TNX, 计算图像坐标与WCS坐标之间的对应关系
 * @version 0.1
 * @date 2017年11月9日
 * - 从FITS文件头加载TNX参数项
 * - 从文本文件加载TNX参数项
 * - 将TNX参数项写入FITS文件
 * - 计算(x,y)对应的WCS坐标(ra, dec)
 */

#ifndef WCSTNX_H_
#define WCSTNX_H_

#include <string>
#include "ADefine.h"

using std::string;
using namespace AstroUtil;

class WCSTNX {
public:
	WCSTNX();
	virtual ~WCSTNX();

public:
	/* 数据结构 */
	enum {// 畸变修正函数类型
		TNX_CHEBYSHEV = 1,	//< 契比雪夫
		TNX_LEGENDRE,		//< 勒让德
		TNX_LINEAR			//< 线性
	};

	enum {// 多项式交叉系数类型
		TNX_XNONE,		//< 无交叉项
		TNX_XFULL,		//< 全交叉
		TNX_XHALF		//< 半交叉
	};

	struct wcs_tnx {// TNX投影修正模式, 单轴wcs坐标拟合参数
		int function;		//< 函数类型
		int xorder, yorder;	//< 阶次
		int xterm;			//< 交叉项类型
		double xmin, xmax, ymin, ymax;	//< 归一化范围
		int ncoef;			//< 系数数量
		double *coef;		//< 系数
		double *x, *y;		//< 多项式单项变量存储区

	private:
		/*!
		 * @brief 计算交叉项数量
		 * @return
		 * 交叉项数量
		 */
		int xterm_count() {
			int order = xorder < yorder ? xorder : yorder;
			int n;

			if      (xterm == TNX_XNONE) n = xorder + yorder - 1;
			else if (xterm == TNX_XFULL) n = xorder * yorder;
			else if (xterm == TNX_XHALF) n = xorder * yorder - order * (order - 1) / 2;

			return n;
		}

		/*!
		 * @brief 生成一元线性数组
		 * @param value  自变量
		 * @param order  阶次
		 * @param ptr    输出数组
		 */
		void linear_array(double value, int order, double *ptr) {
			int i;

			ptr[0] = 1.0;
			for (i = 1; i < order; ++i) ptr[i] = value * ptr[i - 1];
		}

		/*!
		 * @brief 生成一元勒让德数组
		 * @param value  自变量
		 * @param min    最小值
		 * @param max    最大值
		 * @param order  阶次
		 * @param ptr    输出数组
		 */
		void legendre_array(double value, double min, double max, int order, double *ptr) {
			int i;
			double norm = (2 * value - (max + min)) / (max - min);

			ptr[0] = 1.0;
			if (order > 1) ptr[1] = norm;
			for (i = 2; i < order; ++i) {
				ptr[i] = ((2 * i - 1) * norm * ptr[i - 1] - (i - 1) * ptr[i - 2]) / i;
			}
		}

		/*!
		 * @brief 生成一元契比雪夫数组
		 * @param value  自变量
		 * @param min    最小值
		 * @param max    最大值
		 * @param order  阶次
		 * @param ptr    输出数组
		 */
		void chebyshev_array(double value, double min, double max, int order, double *ptr) {
			int i;
			double norm = (2 * value - (max + min)) / (max - min);

			ptr[0] = 1.0;
			if (order > 1) ptr[1] = norm;
			for (i = 2; i < order; ++i) 	ptr[i] = 2 * norm * ptr[i - 1] - ptr[i - 2];
		}

		/*!
		 * @brief 指定自变量后, 计算多项式对应因变量
		 * @return
		 * 因变量
		 */
		double polyval() {
			double sum, sum1;
			int maxorder = xorder > yorder ? xorder : yorder;
			int i, j, k, imax(xorder);

			for (j = k = 0, sum = 0.0; j < yorder; ++j) {
				if (j) {
					if (xterm == TNX_XNONE && imax != 1) imax = 1;
					else if (xterm == TNX_XHALF && (j + xorder) > maxorder) --imax;
				}

				for (i = 0, sum1 = 0.0; i < imax; ++i, ++k) sum1 += (coef[k] * x[i]);
				sum += (sum1 * y[j]);
			}
			return sum;
		}

	public:
		wcs_tnx() {
			function = -1;
			xorder = yorder = -1;
			xterm = -1;
			xmin = xmax = ymin = ymax = 0.0;
			ncoef = 0;
			coef = NULL;
			x = y = NULL;
		}

		void free_array(double **array) {
			if ((*array) != NULL) {
				delete [](*array);
				(*array) = NULL;
			}
		}

		~wcs_tnx() {
			free_array(&coef);
			free_array(&x);
			free_array(&y);
		}

		/*
		 * @note set_order()应在set_xterm()之前执行
		 */
		void set_orderx(int order) {
			if (order <= 0) return;

			if (xorder != order) {
				xorder = order;
				free_array(&x);
			}
			if (!x) x = new double[order];
		}

		void set_ordery(int order) {
			if (order <= 0) return;

			if (yorder != order) {
				yorder = order;
				free_array(&y);
			}
			if (!y) y = new double[order];
		}

		/*
		 * @note set_xterm()应在set_order()之后执行
		 */
		void set_xterm(int xt) {
			if (xt < TNX_XNONE || xt > TNX_XHALF) return;

			if (xterm != xt) xterm = xt;
			int n = xterm_count();
			if (n != ncoef) {
				free_array(&coef);
				ncoef = n;
			}
			if (!coef) coef = new double[n];
		}

		/*!
		 * @brief 计算图像坐标(x,y)对应的投影位置
		 * @param vx X轴坐标
		 * @param vy Y轴坐标
		 * @return
		 * 投影位置, 量纲: 弧度
		 */
		double project_reverse(double vx, double vy) {
			if (!(x && y && coef)) return 0.0;

			if (function == TNX_CHEBYSHEV) {
				chebyshev_array(vx, xmin, xmax, xorder, x);
				chebyshev_array(vy, ymin, ymax, yorder, y);
			}
			else if (function == TNX_LEGENDRE) {
				legendre_array(vx, xmin, xmax, xorder, x);
				legendre_array(vy, ymin, ymax, yorder, y);
			}
			else if (function == TNX_LINEAR) {
				linear_array(vx, xorder, x);
				linear_array(vy, xorder, y);
			}

			return (polyval() * AS2R);
		}
	};

	struct param_tnx {// TNX参数
		PT2F ref_xy;			//< 参考点: XY坐标
		PT2F ref_wcs;		//< 参考点: WCS坐标, 量纲: 弧度
		double cd[2][2];		//< 旋转矩阵. 由平均关系获得. 量纲: 弧度/像素
		bool valid[2];		//< 修正模型有效性
		wcs_tnx tnx1[2];		//< 一阶投影面修正模型
		wcs_tnx tnx2[2];		//< 残差投影面修正模型
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
	 *  0: 正确
	 * -1: 不能打开FITS文件
	 * -2: 非TNX标准WCS信息
	 * -3: 缺少一阶修正模型
	 * -4: 残差修正模型格式错误
	 */
	int LoadImage(const char* filepath);
	/*!
	 * @brief 从文本文件加载WCS参数
	 * @param filepath 文本文件路径
	 * @return
	 * 参数加载结果
	 */
	bool LoadText(const char* filepath);
	/*!
	 * @brief 将LoadText加载的TNX参数写入filepath指代的FITS文件
	 * @param filepath FITS文件路径
	 * @return
	 * 操作结果.
	 * 2017-11-15
	 *   0: 成功
	 *  -1: 未加载位置定标文件, 位置定标必须为TNX格式
	 *  -2: 不能打开FITS文件
	 * 其它: fitsio错误
	 */
	int WriteImage(const char* filepath);
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
	/*!
	 * @brief 输出参数指针
	 * @return
	 * 参数指针
	 */
	param_tnx *GetParam() const;

protected:
	/* 功能 */
	/*!
	 * @brief 图像坐标转换为投影平面平坐标
	 * @param x    图像X坐标
	 * @param y    图像Y坐标
	 * @param xi   投影xi坐标
	 * @param eta  投影eta坐标
	 */
	void image_to_plane(double x, double y, double& xi, double& eta);
	/*!
	 * @brief 由投影平面坐标转换为赤道坐标
	 * @param xi   投影xi坐标, 量纲: 弧度
	 * @param eta  投影eta坐标, 量纲: 弧度
	 * @param ra   赤经, 量纲: 弧度
	 * @param dec  赤纬, 量纲: 弧度
	 */
	void plane_to_wcs(double xi, double eta, double &ra, double &dec);
	/*!
	 * @brief 最大可能保留浮点数精度, 将浮点数转换为字符串
	 * @param output 输出字符串
	 * @param value  浮点数
	 * @return
	 * 转换后字符串长度
	 */
	int output_precision_double(char *output, double value);
	/*!
	 * @brief 解析FITS头中以字符串记录的TNX修正模型
	 * @param strcor 字符串
	 * @param tnx    模型参数存储区
	 * @return
	 * 解析结果
	 */
	int resolve_tnxaxis(char *strcor, wcs_tnx *tnx);
};

#endif /* WCSTNX_H_ */
