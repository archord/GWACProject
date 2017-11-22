/*
 * @file WCSTNX.cpp 定义文件, 基于非标准WCS格式TNX, 计算图像坐标与WCS坐标之间的对应关系
 * @version 0.1
 * @date 2017年11月9日
 * - 从FITS文件头加载WCS TNX参数项
 * - 从文本文件加载WCS TNX参数项
 * - 计算(x,y)对应的WCS坐标(ra, dec)
 */

#include <stdio.h>
#include <string.h>
#include <strings.h>
#include "fits_handler.h"
#include "WCSTNX.h"

WCSTNX::WCSTNX() {

}

WCSTNX::~WCSTNX() {

}

int WCSTNX::LoadImage(const char* filepath) {
	fits_handler hfits;
	if (!hfits(filepath)) return -1; // FITS文件打开失败

	char keyword[10], value[70], CTYPE1[10], CTYPE2[10];
	char strcor[2][2048], (*ptr)[2048];
	int status(0), i, j, n;

	param_.valid[0] = param_.valid[1] = false;

	// 检查TNX关键字
	fits_read_key(hfits(), TSTRING, "CTYPE1", CTYPE1, NULL, &status);
	fits_read_key(hfits(), TSTRING, "CTYPE2", CTYPE2, NULL, &status);
	if (status || strcasecmp(CTYPE1, "RA---TNX") || strcasecmp(CTYPE2, "DEC--TNX"))
		return -2;

	// 检查一阶修正模型
	fits_read_key(hfits(), TDOUBLE, "CRVAL1", &param_.ref_wcs.x, NULL, &status);
	fits_read_key(hfits(), TDOUBLE, "CRVAL2", &param_.ref_wcs.y, NULL, &status);
	fits_read_key(hfits(), TDOUBLE, "CRPIX1", &param_.ref_xy.x,  NULL, &status);
	fits_read_key(hfits(), TDOUBLE, "CRPIX2", &param_.ref_xy.y,  NULL, &status);
	fits_read_key(hfits(), TDOUBLE, "CD1_1",  &param_.cd[0][0],  NULL, &status);
	fits_read_key(hfits(), TDOUBLE, "CD1_2",  &param_.cd[0][1],  NULL, &status);
	fits_read_key(hfits(), TDOUBLE, "CD2_1",  &param_.cd[1][0],  NULL, &status);
	fits_read_key(hfits(), TDOUBLE, "CD2_2",  &param_.cd[1][1],  NULL, &status);
	if (status)
		return -3;
	param_.ref_wcs.x *= D2R;
	param_.ref_wcs.y *= D2R;

	fits_read_key(hfits(), TSTRING, "WAT1_001", value, NULL, &status);
	fits_read_key(hfits(), TSTRING, "WAT2_001", value, NULL, &status);
	if (status) return 0; // 无残差修正模型

	// 解析残差修正模型
	for (j = 1; j <= 2; ++j) {
		ptr = &strcor[j - 1];
		i = n = 0;
		while(!status) {
			sprintf(keyword, "WAT%d_%03d", j, ++i);
			fits_read_key(hfits(), TSTRING, keyword, value, NULL, &status);
			if (!status) n += sprintf((*ptr) + n, "%s", value);
		}
		status = 0;
	}
	param_.valid[1] = !(resolve_tnxaxis(&strcor[0][0], &param_.tnx2[0])
			|| resolve_tnxaxis(&strcor[1][0], &param_.tnx2[1]));

	return param_.valid[1] ? 0 : -4;
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

		if      (!strcmp(token, "xpixref")) param_.ref_xy.x = atof(strtok(NULL, seps));
		else if (!strcmp(token, "ypixref")) param_.ref_xy.y = atof(strtok(NULL, seps));
		else if (!strcmp(token, "lngref"))  param_.ref_wcs.x = atof(strtok(NULL, seps)) * D2R;
		else if (!strcmp(token, "latref"))  param_.ref_wcs.y = atof(strtok(NULL, seps)) * D2R;
		else if (token[0] == 's') {
			int srfc = !strcmp(token, "surface1") ? 1 : (!strcmp(token, "surface2") ? 2 : 0);
			if (srfc) {
				wcs_tnx *tnx = srfc == 1 ? &param_.tnx1[0] : &param_.tnx2[0];
				int i, j, n;

				n = atoi(strtok(NULL, seps));
				for (i = j = 0; i < n && !feof(fp); ++i) {
					fgets(line, size, fp);
					token = strtok(line, seps);
					token1= strtok(NULL, seps);

					if (i > 7) {
						tnx[0].coef[j] = atof(token);
						tnx[1].coef[j] = atof(token1);
						++j;
					}
					else if (i == 0) {
						tnx[0].function = int(atof(token)  + 0.5);
						tnx[1].function = int(atof(token1) + 0.5);
					}
					else if (i == 1) {
						tnx[0].set_orderx(int(atof(token)  + 0.5));
						tnx[1].set_orderx(int(atof(token1) + 0.5));
					}
					else if (i == 2) {
						tnx[0].set_ordery(int(atof(token)  + 0.5));
						tnx[1].set_ordery(int(atof(token1) + 0.5));
					}
					else if (i == 3) {
						tnx[0].set_xterm(int(atof(token)  + 0.5));
						tnx[1].set_xterm(int(atof(token1) + 0.5));
					}
					else if (i == 4) { tnx[0].xmin = atof(token); tnx[1].xmin = atof(token1); }
					else if (i == 5) { tnx[0].xmax = atof(token); tnx[1].xmax = atof(token1); }
					else if (i == 6) { tnx[0].ymin = atof(token); tnx[1].ymin = atof(token1); }
					else if (i == 7) { tnx[0].ymax = atof(token); tnx[1].ymax = atof(token1); }
				}
				param_.valid[srfc - 1] = n && i == n
						&& TNX_CHEBYSHEV <= tnx[0].function && tnx[0].function <= TNX_LINEAR
						&& TNX_CHEBYSHEV <= tnx[1].function && tnx[1].function <= TNX_LINEAR
						&& tnx[0].xorder > 0 && tnx[0].yorder > 0
						&& tnx[1].xorder > 0 && tnx[1].yorder > 0
						&& TNX_XNONE <= tnx[0].xterm && tnx[0].xterm <= TNX_XHALF
						&& TNX_XNONE <= tnx[1].xterm && tnx[1].xterm <= TNX_XHALF
						&& tnx[0].xmax > tnx[0].xmin && tnx[0].ymax > tnx[0].ymin
						&& tnx[1].xmax > tnx[1].xmin && tnx[1].ymax > tnx[1].ymin;
			}
		}
	}
	fclose(fp);

	if (param_.valid[0]) {// 计算旋转矩阵
		wcs_tnx *tnx = &param_.tnx1[0];
		param_.cd[0][0] = 2 * tnx[0].coef[1] * AS2D / (tnx[0].xmax - tnx[0].xmin);
		param_.cd[0][1] = 2 * tnx[0].coef[2] * AS2D / (tnx[0].ymax - tnx[0].ymin);
		param_.cd[1][0] = 2 * tnx[1].coef[1] * AS2D / (tnx[1].xmax - tnx[1].xmin);
		param_.cd[1][1] = 2 * tnx[1].coef[2] * AS2D / (tnx[1].ymax - tnx[1].ymin);
	}

	return (param_.valid[0]);
}

int WCSTNX::WriteImage(const char* filepath) {
	if (!param_.valid[0]) return -1; // 至少需要线性项
	fits_handler hfits;
	if (!hfits(filepath, 1)) return -2; // FITS文件打开失败

	string WCSASTRM = "ct4m.19990714T012701 (USNO-K V) by F. Valdes 1999-08-02";
	string CTYPE1 = "RA---TNX";
	string CTYPE2 = "DEC--TNX";
	string WAT0_001 = "system=image";
	string WAT1 = "wtype=tnx axtype=ra lngcor = ";
	string WAT2 = "wtype=tnx axtype=dec latcor = ";
	int WCSDIM = 2;
	double CRVAL1 = param_.ref_wcs.x * R2D;
	double CRVAL2 = param_.ref_wcs.y * R2D;
	double CRPIX1 = param_.ref_xy.x;
	double CRPIX2 = param_.ref_xy.y;
	double CD1_1  = param_.cd[0][0];
	double CD1_2  = param_.cd[0][1];
	double CD2_1  = param_.cd[1][0];
	double CD2_2  = param_.cd[1][1];
	int status(0);

	if (param_.valid[1]) {
		wcs_tnx *tnx = &param_.tnx2[0];
		int n, nc, i, j;
		char strcor[2][2048];
		char txtdbl[20];
		char (*ptr)[2048];
		for (j = 0, tnx = &param_.tnx2[0]; j < 2; ++j, ++tnx) {
			ptr = &strcor[j];
			n = sprintf(*ptr, "\"%d %d %d %d ", tnx->function, tnx->xorder, tnx->yorder, tnx->xterm);
			output_precision_double(txtdbl, tnx->xmin);
			n += sprintf((*ptr) + n, "%s ", txtdbl);
			output_precision_double(txtdbl, tnx->xmax);
			n += sprintf((*ptr) + n, "%s ", txtdbl);
			output_precision_double(txtdbl, tnx->ymin);
			n += sprintf((*ptr) + n, "%s ", txtdbl);
			output_precision_double(txtdbl, tnx->ymax);
			n += sprintf((*ptr) + n, "%s ", txtdbl);

			nc = tnx->ncoef;
			for (i = 0; i < nc; ++i) {
				output_precision_double(txtdbl, tnx->coef[i]);
				n += sprintf((*ptr) + n, "%s ", txtdbl);
			}
			n += sprintf((*ptr) + n, "\"");
			(*ptr)[n] = 0;
		}

		// 合并修正项
		WAT1 += strcor[0];
		WAT2 += strcor[1];
	}

	// 写入FITS头
	fits_update_key(hfits(), TSTRING, "WCSASTRM", (void*) WCSASTRM.c_str(), "WCS Source",         &status);
	fits_update_key(hfits(), TINT,    "WCSDIM",   &WCSDIM,                  "WCS dimensionality", &status);
	fits_update_key(hfits(), TSTRING, "CTYPE1",   (void*) CTYPE1.c_str(),   "Coordinate type",    &status);
	fits_update_key(hfits(), TSTRING, "CTYPE2",   (void*) CTYPE2.c_str(),   "Coordinate type",    &status);
	fits_update_key(hfits(), TDOUBLE, "CRVAL1",   &CRVAL1, "Coordinate reference value", &status);
	fits_update_key(hfits(), TDOUBLE, "CRVAL2",   &CRVAL2, "Coordinate reference value", &status);
	fits_update_key(hfits(), TDOUBLE, "CRPIX1",   &CRPIX1, "Coordinate reference pixel", &status);
	fits_update_key(hfits(), TDOUBLE, "CRPIX2",   &CRPIX2, "Coordinate reference pixel", &status);
	fits_update_key(hfits(), TDOUBLE, "CD1_1",    &CD1_1,  "Coordinate matrix",          &status);
	fits_update_key(hfits(), TDOUBLE, "CD1_2",    &CD1_2,  "Coordinate matrix",          &status);
	fits_update_key(hfits(), TDOUBLE, "CD2_1",    &CD2_1,  "Coordinate matrix",          &status);
	fits_update_key(hfits(), TDOUBLE, "CD2_2",    &CD2_2,  "Coordinate matrix",          &status);
	if (param_.valid[1]) {// 畸变改正项
		char item[70]; // 每行实际可存储数据68字节
		char keyword[10]; // 关键字
		int i, j, byteleft, byteitem(68);
		int len[] = { WAT1.size(), WAT2.size() };
		const char *head[] = { WAT1.c_str(), WAT2.c_str() };
		const char *ptr;

		fits_update_key(hfits(), TSTRING, "WAT0_001", (void*) WAT0_001.c_str(), "Coordinate system", &status);
		for (j = 0; j < 2; ++j) {// 0: R.A.; 1: DEC.
			for (byteleft = len[j], i = 0, ptr = head[j]; byteleft; ptr += byteitem) {
				sprintf(keyword, "WAT%d_%03d", j + 1, ++i);
				if (byteleft > byteitem) {
					strncpy(item, ptr, byteitem);
					item[byteitem] = 0;
					fits_update_key(hfits(), TSTRING, keyword, (void*) item, "", &status);
				}
				else {
					fits_update_key(hfits(), TSTRING, keyword, (void*) ptr, "", &status);
				}
				byteleft = byteleft > byteitem ? byteleft - byteitem : 0;
			}
		}
	}

	return status;
}

int WCSTNX::XY2WCS(double x, double y, double& ra, double& dec) {
	double xi, eta;

	image_to_plane(x, y, xi, eta);
	if (param_.valid[1]) {
		xi  += param_.tnx2[0].project_reverse(x, y);
		eta += param_.tnx2[1].project_reverse(x, y);
	}
	plane_to_wcs(xi, eta, ra, dec);

	return 0;
}

WCSTNX::param_tnx *WCSTNX::GetParam() const {
	return (param_tnx * const) &param_;
}

/*
 * @note 图像坐标转换为投影(中间)坐标
 */
void WCSTNX::image_to_plane(double x, double y, double& xi, double& eta) {
	if (param_.valid[0]) {
		xi  = param_.tnx1[0].project_reverse(x, y);
		eta = param_.tnx1[1].project_reverse(x, y);
	}
	else {
		double dx(x - param_.ref_xy.x), dy(y - param_.ref_xy.y);
		xi  = (param_.cd[0][0] * dx + param_.cd[0][1] * dy) * D2R;
		eta = (param_.cd[1][0] * dx + param_.cd[1][1] * dy) * D2R;
	}
}

/*
 * @note 投影坐标转换为WCS(赤道)坐标
 */
void WCSTNX::plane_to_wcs(double xi, double eta, double &ra, double &dec) {
	double A0(param_.ref_wcs.x), D0(param_.ref_wcs.y);
	double fract = cos(D0) - eta * sin(D0);
	ra  = reduce(A0 + atan2(xi, fract), A2PI);
	dec = atan2(((eta * cos(D0) + sin(D0)) * cos(ra - A0)), fract);
}

int WCSTNX::output_precision_double(char *output, double value) {
	int n(0), pos(0), valid(0), intv;
	if (value < 0) {
		output[0] = '-';
		++pos;
		value = -value;
	}

	// 整数部分
	intv = int(value + AEPS);
	n = sprintf(output + pos, "%d.", intv);
	pos += n;
	valid = n - 1;
	// 小数部分
	while(valid < 17 && (value = (value - intv) * 10) > AEPS) {
		intv = int(value);
		sprintf(output + pos, "%d", intv);
		++pos;
		++valid;
	}
	output[pos] = 0;

	return n;
}

int WCSTNX::resolve_tnxaxis(char *strcor, wcs_tnx *tnx)  {
	char	 *pstr, *ptr;
	const char seps[] = " ";
	const char flags[] = "1234567890-+.";
	int xorder, yorder, xterm, i;

	if (NULL == (pstr = strpbrk(strcor, flags)))
		return 1;

	tnx->function = int(atof(strtok_r(pstr, seps, &ptr)) + 0.5);
	xorder        = (pstr = strtok_r(NULL, seps, &ptr)) != NULL ? int(atof(pstr) + 0.5) : 0;
	yorder        = (pstr = strtok_r(NULL, seps, &ptr)) != NULL ? int(atof(pstr) + 0.5) : 0;
	xterm         = (pstr = strtok_r(NULL, seps, &ptr)) != NULL ? int(atof(pstr) + 0.5) : 0;
	tnx->xmin     = (pstr = strtok_r(NULL, seps, &ptr)) != NULL ? atof(pstr) : 0.0;
	tnx->xmax     = (pstr = strtok_r(NULL, seps, &ptr)) != NULL ? atof(pstr) : 0.0;
	tnx->ymin     = (pstr = strtok_r(NULL, seps, &ptr)) != NULL ? atof(pstr) : 0.0;
	tnx->ymax     = (pstr = strtok_r(NULL, seps, &ptr)) != NULL ? atof(pstr) : 0.0;
	if (TNX_CHEBYSHEV > tnx->function || TNX_LINEAR < tnx->function
			|| xorder <= 0 || yorder <= 0
			|| TNX_XNONE > xterm || TNX_XHALF < xterm
			|| tnx->xmax <= tnx->xmin || tnx->ymax < tnx->ymin)
		return 2;

	tnx->set_orderx(xorder);
	tnx->set_ordery(yorder);
	tnx->set_xterm(xterm);

	for (i=0; i < tnx->ncoef && (pstr = strtok_r(NULL, seps, &ptr)) != NULL; ++i)
		tnx->coef[i] = atof(pstr);

	return (i == tnx->ncoef) ? 0 : 3;
}
