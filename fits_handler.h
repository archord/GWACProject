/*
 * @file fits_handler.h 基于cfitsio封装FITS文件基本操作
 */

#ifndef FITS_HANDLER_H_
#define FITS_HANDLER_H_

#include <vector>
#include <string>
#include <longnam.h>
#include <fitsio.h>

struct fits_handler {// FITS图像操作接口
	fitsfile *fitsptr;	//< 基于cfitsio接口的文件操作接口
	int width, height;	//< 图像宽度和高度
	float exptime;		//< 曝光时间, 量纲: 秒

public:
	fits_handler() {
		fitsptr = NULL;
		width = height = 0;
		exptime = 0.0;
	}

	bool close() {
		int status(0);
		if (fitsptr) fits_close_file(fitsptr, &status);
		if (!status) fitsptr = NULL;
		return !status;
	}

	virtual ~fits_handler() {
		close();
	}

	/*!
	 * @brief 获得基于cfitsio的FITS文件访问接口
	 * @return
	 * FITS文件访问接口
	 */
	fitsfile* operator()() {
		return fitsptr;
	}

	/*!
	 * @brief 打开或创建FITS文件并获得基于cfitsio的FITS文件访问接口
	 * @param pathname    文件路径
	 * @param mode        0: 只读; 1: 读写; 其它: 新文件
	 * @return
	 * FITS文件访问接口
	 */
	fitsfile* operator()(const std::string &pathname, const int mode = 0) {
		if (!close()) return NULL;
		int status(0);
		if (mode == 0 || mode == 1) {
			int naxis;
			long naxes[2];
			fits_open_file(&fitsptr, pathname.c_str(), mode, &status);
			fits_get_img_dim(fitsptr, &naxis, &status);
			if (naxis != 2) {
				close();
				return NULL;
			}
			fits_get_img_size(fitsptr, naxis, naxes, &status);
			width = naxes[0];
			height= naxes[1];
			fits_read_key_flt(fitsptr, "EXPTIME", &exptime, NULL, &status);
		}
		else fits_create_file(&fitsptr, pathname.c_str(), &status);
		if (status) close();
		return fitsptr;
	}

	/*!
	 * @brief 从文件中读取一行数据到缓冲区
	 * @param row  行编号, 从0开始
	 * @param buff 数据缓冲区
	 * @return
	 * 读取结果
	 */
	bool loadrow(const int row, float* buff) {
		if (!fitsptr) return false;
		int status(0);
		fits_read_img(fitsptr, TFLOAT, row * width + 1, width, NULL, buff, NULL, &status);

		return !status;
	}

	/*!
	 * @brief 加载图像数据到缓冲区
	 * @param buff 数据缓冲区
	 * @return
	 * 读取结果
	 */
	bool loadimg(float* buff) {
		if (!fitsptr) return false;
		int status(0);
		fits_read_img(fitsptr, TFLOAT, 1, width * height, NULL, buff, NULL, &status);
		return !status;
	}
};
typedef fits_handler * HFITS;
typedef std::vector<HFITS> HFITSVEC;

#endif /* FITS_HANDLER_H_ */
