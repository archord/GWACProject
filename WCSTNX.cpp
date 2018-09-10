/*
 * @file WCSTNX.cpp 定义文件, 基于非标准WCS格式TNX, 计算图像坐标与WCS坐标之间的对应关系
 * @version 0.1
 * @date 2017年11月9日
 * - 从FITS文件头加载WCS TNX参数项
 * - 从文本文件加载WCS TNX参数项
 * - 计算(x,y)对应的WCS坐标(ra, dec)
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include "ADefine.h"
#include "WCSTNX.h"

using namespace AstroUtil;

WCSTNX::WCSTNX() {

}

WCSTNX::~WCSTNX() {

}

bool WCSTNX::LoadImage(const char* filepath) {
	return true;
}

bool WCSTNX::LoadText(const char* filepath) {
	FILE *fp = fopen(filepath, "r");
	if (!fp) return false;
	const int size(100);
	char line[size];
	char seps[] = " \t\r\n";
	char *token, *token1;

	while(!feof(fp)) {
		if (fgets(line, size, fp) == NULL) continue;
		token = strtok(line, seps);

		if (token[0] == 'x') {
			if (!strcmp(token, "xrefmean")) param_.ref_xymean.x = atof(strtok(NULL, seps));
			else if (!strcmp(token, "xpixref")) param_.ref_xy.x = atof(strtok(NULL, seps));
			else if (!strcmp(token, "xishift")) param_.shift.x = atof(strtok(NULL, seps));
			else if (!strcmp(token, "xmag")) param_.mag.x = atof(strtok(NULL, seps));
			else if (!strcmp(token, "xrotation")) param_.rotation.x = atof(strtok(NULL, seps));
		}
		else if (token[0] == 'y') {
			if (!strcmp(token, "yrefmean")) param_.ref_xymean.y = atof(strtok(NULL, seps));
			else if (!strcmp(token, "ypixref")) param_.ref_xy.y = atof(strtok(NULL, seps));
			else if (!strcmp(token, "ymag")) param_.mag.y = atof(strtok(NULL, seps));
			else if (!strcmp(token, "yrotation")) param_.rotation.y = atof(strtok(NULL, seps));
		}
		else if (token[0] == 'l') {
			if (!strcmp(token, "lngmean")) param_.ref_wcsmean.x = atof(strtok(NULL, seps));
			else if (!strcmp(token, "latmean")) param_.ref_wcsmean.y = atof(strtok(NULL, seps));
			else if (!strcmp(token, "lngref"))  param_.ref_wcs.x = atof(strtok(NULL, seps));
			else if (!strcmp(token, "latref"))  param_.ref_wcs.y = atof(strtok(NULL, seps));
		}
		else if (token[0] == 'e') {
			if (!strcmp(token, "etashift")) param_.shift.y = atof(strtok(NULL, seps));
		}
		else if (token[0] == 'p') {
			if (!strcmp(token, "pixsystem"))  param_.pixsystem  = strtok(NULL, seps);
			else if (!strcmp(token, "projection")) param_.projection = strtok(NULL, seps);
		}
		else if (token[0] == 's') {
			int n;
			int xxorder, yxorder, xyorder, yyorder;
			int xterm, yterm;

			if (!strcmp(token, "surface1")) {
				param_surface& srfc = param_.surface1;
				int i, j(0);

				n = atoi(strtok(NULL, seps));
				for (i = 0; i < n && !feof(fp); ++i) {
					fgets(line, size, fp);
					token = strtok(line, seps);
					token1= strtok(NULL, seps);

					if (i == 0) { srfc.xsurface = atoi(token); srfc.ysurface = atoi(token1); }
					else if (i == 1) { xxorder = atoi(token);  yxorder = atoi(token1); }
					else if (i == 2) {
						xyorder = atoi(token);
						yyorder = atoi(token1);
						srfc.set_orderx(xxorder, xyorder);
						srfc.set_ordery(yxorder, yyorder);
					}
					else if (i == 3) {
						xterm = atoi(token);
						yterm = atoi(token1);
						srfc.set_xtermx(xterm);
						srfc.set_xtermy(yterm);
					}
					else if (i == 4) { srfc.xxmin = atof(token); srfc.yxmin = atof(token1); }
					else if (i == 5) { srfc.xxmax = atof(token); srfc.yxmax = atof(token1); }
					else if (i == 6) { srfc.xymin = atof(token); srfc.yymin = atof(token1); }
					else if (i == 7) { srfc.xymax = atof(token); srfc.yymax = atof(token1); }
					else {
						srfc.xcoef[j] = atof(token);
						srfc.ycoef[j] = atof(token1);
						++j;
					}
				}
				param_.valid1 = (n && i == n);
				/*----- 验证 -----*/
//				printf("surface1:\n");
//				printf("surface = %d\t%d\n", srfc.xsurface, srfc.ysurface);
//				printf("x order = %d\t%d\n", srfc.xxorder, srfc.xyorder);
//				printf("y order = %d\t%d\n", srfc.yxorder, srfc.yyorder);
//				printf("xterm   = %d\t%d\n", srfc.xxterm, srfc.yxterm);
//				printf("x (%4.0f %4.0f) <==> (%4.0f %4.0f)\n", srfc.xxmin, srfc.xymin, srfc.xxmax, srfc.xymax);
//				printf("y (%4.0f %4.0f) <==> (%4.0f %4.0f)\n", srfc.yxmin, srfc.yymin, srfc.yxmax, srfc.yymax);
//				printf("coef:\n");
//				for (i = 0; i < j; ++i) {
//					printf("\t%f\t%f\n", srfc.xcoef[i], srfc.ycoef[i]);
//				}
			}
			else if (!strcmp(token, "surface2")) {
				param_surface& srfc = param_.surface2;
				int i(0), j(0);

				n = atoi(strtok(NULL, seps));
				for (i = 0; i < n && !feof(fp); ++i) {
					fgets(line, size, fp);
					token = strtok(line, seps);
					token1= strtok(NULL, seps);

					if (i == 0) { srfc.xsurface = atoi(token); srfc.ysurface = atoi(token1); }
					else if (i == 1) { xxorder = atoi(token);  yxorder = atoi(token1); }
					else if (i == 2) {
						xyorder = atoi(token);
						yyorder = atoi(token1);
						srfc.set_orderx(xxorder, xyorder);
						srfc.set_ordery(yxorder, yyorder);
					}
					else if (i == 3) {
						xterm = atoi(token);
						yterm = atoi(token1);
						srfc.set_xtermx(xterm);
						srfc.set_xtermy(yterm);
					}
					else if (i == 4) { srfc.xxmin = atof(token); srfc.yxmin = atof(token1); }
					else if (i == 5) { srfc.xxmax = atof(token); srfc.yxmax = atof(token1); }
					else if (i == 6) { srfc.xymin = atof(token); srfc.yymin = atof(token1); }
					else if (i == 7) { srfc.xymax = atof(token); srfc.yymax = atof(token1); }
					else {
						srfc.xcoef[j] = atof(token);
						srfc.ycoef[j] = atof(token1);
						++j;
					}
				}
				param_.valid2 = (n && i == n);
				/*----- 验证 -----*/
//				printf("surface2:\n");
//				printf("surface = %d\t%d\n", srfc.xsurface, srfc.ysurface);
//				printf("x order = %d\t%d\n", srfc.xxorder, srfc.xyorder);
//				printf("y order = %d\t%d\n", srfc.yxorder, srfc.yyorder);
//				printf("xterm   = %d\t%d\n", srfc.xxterm, srfc.yxterm);
//				printf("x (%4.0f %4.0f) <==> (%4.0f %4.0f)\n", srfc.xxmin, srfc.xymin, srfc.xxmax, srfc.xymax);
//				printf("y (%4.0f %4.0f) <==> (%4.0f %4.0f)\n", srfc.yxmin, srfc.yymin, srfc.yxmax, srfc.yymax);
//				printf("coef:\n");
//				for (i = 0; i < j; ++i) {
//					printf("\t%f\t%f\n", srfc.xcoef[i], srfc.ycoef[i]);
//				}
			}
		}
		else if (!strcmp(token, "coosystem"))  param_.coosystem  = strtok(NULL, seps);
		else if (!strcmp(token, "function"))   param_.function = strtok(NULL, seps);
	}

	if (param_.valid1) {// 计算旋转矩阵
//		double xrot = param_.rotation.x * D2R - API;
//		double yrot = param_.rotation.y * D2R - API;
//		double xs = param_.mag.x * AS2D;
//		double ys = param_.mag.y * AS2D;
//		param_.cd[0][0] =  xs * cos(xrot);
//		param_.cd[0][1] = -ys * sin(yrot);
//		param_.cd[1][0] =  xs * sin(xrot);
//		param_.cd[1][1] =  ys * cos(yrot);
	}

	fclose(fp);
	return true;
}

int WCSTNX::XY2WCS(double x, double y, double& ra, double& dec) {
	if (!param_.valid1) return -1;
	double xi, eta, xi1, eta1;

	image_to_plane(x, y, xi, eta);
	if (param_.valid2) {
		correct(param_.surface2, x, y, xi1, eta1);
		xi += xi1;
		eta += eta1;
	}
	plane_to_wcs(xi, eta, ra, dec);

	return 0;
}

int WCSTNX::WCS2XY(double ra, double dec, double& x, double& y) {
	if (!param_.valid1) return -1;
	double xi, eta, xi1, eta1;

	wcs_to_plane(ra, dec, xi, eta);
        correct(param_.surface1, xi, eta, x, y);
	if (param_.valid2) {
		correct(param_.surface2, xi, eta, xi1, eta1);
		x += xi1;
		y += eta1;
	}

	return 0;
}

void WCSTNX::linear_array(double x, int order, double* ptr) {
	int i;

	ptr[0] = 1.0;
	for (i = 1; i < order; ++i) {
		ptr[i] = x * ptr[i - 1];
	}
}

void WCSTNX::legendre_array(double x, double xmin, double xmax, int order, double* ptr) {
	int i;
	double xnorm = (2 * x - (xmax + xmin)) / (xmax - xmin);

	ptr[0] = 1.0;
	if (order > 1) ptr[1] = xnorm;
	for (i = 2; i < order; ++i) {
		ptr[i] = ((2 * i - 1) * xnorm * ptr[i - 1] - (i - 1) * ptr[i - 2]) / i;
	}
}

void WCSTNX::chebyshev_array(double x, double xmin, double xmax, int order, double* ptr) {
	int i;
	double xnorm = (2 * x - (xmax + xmin)) / (xmax - xmin);

	ptr[0] = 1.0;
	if (order > 1) ptr[1] = xnorm;
	for (i = 2; i < order; ++i) {
		ptr[i] = 2 * xnorm * ptr[i - 1] - ptr[i - 2];
	}
}

void WCSTNX::polyval_item(int type, int xorder, int yorder, double* x, double* y, int n, double* array) {
	double *ptr, t;
	int maxorder = xorder > yorder ? xorder : yorder;
	int i, j, imin(0), imax(xorder);

	for (j = 0, ptr = array; j < yorder; ++j) {
		if (j) {
			if (type == TNX_XNONE && imax != 1) imax = 1;
			else if (type == TNX_XHALF && (j + xorder) > maxorder) --imax;
		}

		for (i = imin, t = y[j]; i < imax; ++i, ++ptr) *ptr = x[i] * t;
	}
}

double WCSTNX::polysum(int n, double* coef, double* item) {
	double sum(0.0);

	for (int i = 0; i < n; ++i) sum += (coef[i] * item[i]);
	return sum;
}

void WCSTNX::correct(param_surface& surface, double x, double y, double& dx, double& dy) {
	// xi残差
	int xsurface(surface.xsurface), xxterm(surface.xxterm);
	int xxorder(surface.xxorder), xyorder(surface.xyorder);
	double xxmin(surface.xxmin), xxmax(surface.xxmax);
	double xymin(surface.xymin), xymax(surface.xymax);
	double *xx = surface.xx;
	double *xy = surface.xy;
	double *xxy = surface.xxy;

	if (xsurface == TNX_CHEB) {
		chebyshev_array(x, xxmin, xxmax, xxorder, xx);
		chebyshev_array(y, xymin, xymax, xyorder, xy);
	}
	else if (xsurface == TNX_LEG) {
		legendre_array(x, xxmin, xxmax, xxorder, xx);
		legendre_array(y, xymin, xymax, xyorder, xy);
	}
	else if (xsurface == TNX_POLY) {
		linear_array(x, xxorder, xx);
		linear_array(y, xyorder, xy);
	}
	polyval_item(xxterm, xxorder, xyorder, xx, xy, surface.xncoef, xxy);
	dx = polysum(surface.xncoef, surface.xcoef, xxy);
	// eta残差
	int ysurface(surface.ysurface), yxterm(surface.yxterm);
	int yxorder(surface.yxorder), yyorder(surface.yyorder);
	double yxmin(surface.yxmin), yxmax(surface.yxmax);
	double yymin(surface.yymin), yymax(surface.yymax);
	double *yx = surface.yx;
	double *yy = surface.yy;
	double *yxy = surface.yxy;

	if (ysurface == TNX_CHEB) {
		chebyshev_array(x, yxmin, yxmax, yxorder, yx);
		chebyshev_array(y, yymin, yymax, yyorder, yy);
	}
	else if (ysurface == TNX_LEG) {
		legendre_array(x, yxmin, yxmax, yxorder, yx);
		legendre_array(y, yymin, yymax, yyorder, yy);
	}
	else if (ysurface == TNX_POLY) {
		linear_array(x, yxorder, yx);
		linear_array(y, yyorder, yy);
	}
	polyval_item(yxterm, yxorder, yyorder, yx, yy, surface.yncoef, yxy);
	dy = polysum(surface.yncoef, surface.ycoef, yxy);

	dx *= AS2D;
	dy *= AS2D;
}

/*
 * @note 图像坐标转换为投影(中间)坐标
 */
void WCSTNX::image_to_plane(double x, double y, double& xi, double& eta) {
	if (param_.valid1) correct(param_.surface1, x, y, xi, eta);
	else {
//		double dx(x - param_.ref_xy.x), dy(y - param_.ref_xy.y);
//
//		xi  = param_.cd[0][0] * dx + param_.cd[0][1] * dy;
//		eta = param_.cd[1][0] * dx + param_.cd[1][1] * dy;
	}
}

/*
 * @note TAN投影逆变换, 平面坐标==>球面坐标
 */
void WCSTNX::plane_to_wcs(double x, double y, double& ra, double& dec) {
	double ra0(param_.ref_wcs.x * D2R), dec0(param_.ref_wcs.y * D2R);
	double x1(x * D2R), y1(y * D2R);
	double dra = atan2(x1 / cos(dec0), 1 - y1 * tan(dec0));

	ra  = reduce(ra0 + dra, A2PI) * R2D;
	dec = atan2((y1 + tan(dec0)) * cos(dra), 1 - y1 * tan(dec0)) * R2D;
}


int WCSTNX::wcs_to_plane(double ra, double dec, double &xi, double &eta) {

    double ra0(param_.ref_wcs.x * D2R), dec0(param_.ref_wcs.y * D2R);
    ra *= D2R;
    dec *= D2R;

    int j = 0;
    double tiny = 1e-6;

    double sdec0 = sin(dec0);
    double sdec = sin(dec);
    double cdec0 = cos(dec0);
    double cdec = cos(dec);
    double radif = ra - ra0;
    double sradif = sin(radif);
    double cradif = cos(radif);

    double denom = sdec * sdec0 + cdec * cdec0*cradif;

    if (denom > tiny) {
        j = 0;
    } else if (denom >= 0.0) {
        j = 1;
        denom = tiny;
    } else if (denom > -tiny) {
        j = 2;
        denom = -tiny;
    } else {
        j = 3;
    }

    xi = cdec * sradif / denom;
    eta = (sdec * cdec0 - cdec * sdec0 * cradif) / denom;

    return j;
}
