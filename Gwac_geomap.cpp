
#include "gwac.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mathFunction.h"

int cofun(double x1, double x2, double *afunc, int cofNum);
int printBasicInfo(const char *fName, char *otype, double xrefmean, double yrefmean, double xmean,
        double ymean, double xrms2, double yrms2, double xrmsn, double yrmsn,
        double *cof2x, double * cof2y, int cof2Num, double *cofnx, double * cofny,
        int cofnNum, char statusstr[]);
int readCof(char *fName, double *xcof, double *ycof, int cofNum);
int GWAC_fitting(double *x1, double *y1, double *x2, double *y2, int dataNum,
        double *ax, double *ay, int cofNum, double *xrms, double *yrms,
        int iter, float rejsigma, char statusstr[]);
int printFitDiff(const char *fName, double refx[], double refy[], double inx[],
        double iny[], int dataNum, double ax[], double ay[], int cofNum,
        char statusstr[]);

/*******************************************************************************
 * 功能：星表坐标拟合，包括正向拟合(in到ref)和反向拟合两部分(ref到in)
 *      1，实现iraf中geomap的功能，对输入的一组星表的坐标值in(x,y)进行坐标转换到另一组
 *      星表ref(x,y)。其中geomap的实现是先对这两组星表坐标进行一次线性拟合，然后再对ref
 *      与拟合结果的残差进行高阶多项式拟合。
 *      2，这里实现时，如果order大于1，首先对这两组星表坐标进行一次线性拟合，然后再直接
 *      对这两组星表坐标进行高阶多项式拟合，没有针对惨差的高阶多项式拟合。如果order等于1，
 *      则只进行一阶线性拟合。
 *      3，拟合完成后，按照geomap的输出文件格式对结果进行输出。
 *      4，该函数首先进行正向拟合，拟合结果输出到outfilename文件中；然后进行反向拟合，拟
 *      合结果附加到outfilename文件中。
 * 
 **输入：
 *      matchpeervec   亮星匹配对
 *      order        拟合阶数
 *      iter         迭代次数
 *      rejsigma     数据抛弃sigma值
 * 
 **输出
 *      xrms         正向拟合x方向拟合精度
 *      yrms         正向拟合y方向拟合精度
 *      xrms2        反向拟合x方向拟合精度
 *      yrms2        反向拟合y方向拟合精度
 *      outfilename  输出拟合文件，格式同iraf geomap输出，具体描述见
 *      statusstr    返回值说明（不能超过128个字符）
 * 
 **返回值:
 *      0表示正确，其它值为错误码，
 *      蔡使用2001～2999，徐使用3001～3999，苑使用4001～4999，李使用5001～5999
 ******************************************************************************/
int Gwac_geomap(vector<ST_STARPEER> matchpeervec,
        unsigned int order,
        unsigned int iter,
        float rejsigma,
        float &xrms,
        float &yrms,
        float &xrms2,
        float &yrms2,
        const char outfilename[],
        char statusstr[]) {

    int retStatus1, retStatus2, retStatus3, retStatus4, retStatus5, retStatus6;

    /*检查错误结果输出参数statusstr是否为空*/
    CHECK_STATUS_STR_IS_NULL(statusstr);
    /*检测星表数据数组是否为空*/
    if (matchpeervec.empty()) {
        sprintf(statusstr, "Error Code: %d\n"
                "File %s line %d, the input parameter matchpeervec is empty!\n",
                GWAC_FUNCTION_INPUT_EMPTY, __FILE__, __LINE__);
        return GWAC_FUNCTION_INPUT_EMPTY;
    }
    /*检测输出结果文件名是否为空*/
    CHECK_STRING_NULL_OR_EMPTY(outfilename, "outfilename");

    int pointNum = matchpeervec.size(); //数据点的个数
    int cofNum = (order + 1)*(order + 2) / 2; //高阶拟合系数个数
    int lineCof = 3; //一阶拟合系数个数

    /**
     * pointNum + 1， 这里面的多项式拟合调用的《C数值算法》一书的实现，但是该书的C语言函
     * 数是直接基于fortran转变过来的，而fortran 的数组都是从1开始，为了保持统一，这里多
     * 申请1个单位的内存空间，所以该实现的内存分配都有+1的现象，以及for循环是从1开始的。
     */
    double *refx = (double *) malloc((pointNum + 1) * sizeof (double));
    double *refy = (double *) malloc((pointNum + 1) * sizeof (double));
    double *objx = (double *) malloc((pointNum + 1) * sizeof (double));
    double *objy = (double *) malloc((pointNum + 1) * sizeof (double));

    double *ax = (double *) malloc((cofNum + 1) * sizeof (double));
    double *ay = (double *) malloc((cofNum + 1) * sizeof (double));
    double *lineax = (double *) malloc((lineCof + 1) * sizeof (double));
    double *lineay = (double *) malloc((lineCof + 1) * sizeof (double));

    double *afunc = (double *) malloc((cofNum + 1) * sizeof (double));

    if (refx == NULL || refy == NULL || objx == NULL || objy == NULL || ax == NULL ||
            ay == NULL || lineax == NULL || lineay == NULL || afunc == NULL) {
        if (refx != NULL) free(refx);
        if (refy != NULL) free(refy);
        if (objx != NULL) free(objx);
        if (objy != NULL) free(objy);
        if (ax != NULL) free(ax);
        if (ay != NULL) free(ay);
        if (lineax != NULL) free(lineax);
        if (lineay != NULL) free(lineay);
        if (afunc != NULL) free(afunc);
        MALLOC_IS_NULL();
    }

    double objxrms2 = 0.0, objyrms2 = 0.0, objxrmsn = 0.0, objyrmsn = 0.0;
    double xrefmean = 0.0, yrefmean = 0.0, xmean = 0.0, ymean = 0.0;

    memset(ax, 0, (cofNum + 1) * sizeof (double));
    memset(ay, 0, (cofNum + 1) * sizeof (double));
    memset(lineax, 0, (lineCof + 1) * sizeof (double));
    memset(lineay, 0, (lineCof + 1) * sizeof (double));

    /*正向拟合，使用参考星拟合目标星*/
    int i, j;
    for (i = 0; i < pointNum; i++) {
        ST_STAR tref = matchpeervec.at(i).ref;
        ST_STAR tobj = matchpeervec.at(i).obj;
        refx[i + 1] = tref.x;
        refy[i + 1] = tref.y;
        objx[i + 1] = tobj.x;
        objy[i + 1] = tobj.y;
    }

    xrefmean = meand(refx, pointNum);
    yrefmean = meand(refy, pointNum);
    xmean = meand(objx, pointNum);
    ymean = meand(objy, pointNum);

    /*首先进行一阶线性拟合*/
    retStatus1 = GWAC_fitting(objx, objy, refx, refy, pointNum, lineax,
            lineay, lineCof, &objxrms2, &objyrms2, iter, rejsigma, statusstr);

    xrms = objxrms2;
    yrms = objyrms2;

#ifdef GWAC_TEST
    printFitDiff("order2diff.txt", refx, refy, objx, objy, pointNum,
            lineax, lineay, lineCof, statusstr);
#endif

    /*如果order大于1，再进行高阶拟合*/
    if (order > 1) {
#ifdef RESIDUALS
        /*对目标值与拟合值求残差，该程序对残差的拟合效果不好*/
        for (i = 1; i <= pointNum; i++) {
            cofun(objx[i], objy[i], afunc, lineCof);
            double tmpx = 0.0;
            double tmpy = 0.0;
            for (j = 1; j <= lineCof; j++) {
                tmpx += lineax[j] * afunc[j];
                tmpy += lineay[j] * afunc[j];
            }
            refx[i] = fabs(refx[i] - tmpx);
            refy[i] = fabs(refy[i] - tmpy);
            //            refx[i] = (refx[i] - tmpx);
            //            refy[i] = (refy[i] - tmpy);
        }
#endif
        retStatus2 = GWAC_fitting(objx, objy, refx, refy, pointNum, ax, ay, cofNum,
                &objxrmsn, &objyrmsn, iter, rejsigma, statusstr);

        xrms = objxrmsn;
        yrms = objyrmsn;

#ifdef GWAC_TEST
        printFitDiff("orderndiff.txt", refx, refy, objx, objy, pointNum,
                ax, ay, cofNum, statusstr);
#endif
    }

    /*将拟合结果输出到文件*/
    retStatus3 = printBasicInfo(outfilename, "w", xrefmean, yrefmean, xmean, ymean, objxrms2,
            objyrms2, objxrmsn, objyrmsn, lineax, lineay, lineCof, ax, ay,
            cofNum, statusstr);

    /*反向拟合，使用目标星拟合参考星*/
    for (i = 0; i < pointNum; i++) {
        ST_STAR tref = matchpeervec.at(i).ref;
        ST_STAR tobj = matchpeervec.at(i).obj;
        refx[i + 1] = tobj.x;
        refy[i + 1] = tobj.y;
        objx[i + 1] = tref.x;
        objy[i + 1] = tref.y;
    }

    xrefmean = meand(refx, pointNum);
    yrefmean = meand(refy, pointNum);
    xmean = meand(objx, pointNum);
    ymean = meand(objy, pointNum);

    /*首先进行一阶线性拟合*/
    retStatus4 = GWAC_fitting(objx, objy, refx, refy, pointNum, lineax, lineay, lineCof,
            &objxrms2, &objyrms2, iter, rejsigma, statusstr);

    xrms2 = objxrms2;
    yrms2 = objyrms2;

#ifdef GWAC_TEST
    printFitDiff("reverse_order2diff.txt", refx, refy, objx, objy, pointNum,
            lineax, lineay, lineCof, statusstr);
#endif

    /*如果order大于1，再进行高阶拟合*/
    if (order > 1) {
#ifdef RESIDUALS
        /*对目标值与拟合值求残差，该程序对残差的拟合效果不好*/
        for (i = 1; i <= pointNum; i++) {
            cofun(objx[i], objy[i], afunc, lineCof);
            double tmpx = 0.0;
            double tmpy = 0.0;
            for (j = 1; j <= lineCof; j++) {
                tmpx += lineax[j] * afunc[j];
                tmpy += lineay[j] * afunc[j];
            }
            refx[i] = fabs(refx[i] - tmpx);
            refy[i] = fabs(refy[i] - tmpy);
            //            refx[i] = (refx[i] - tmpx);
            //            refy[i] = (refy[i] - tmpy);
        }
#endif
        retStatus5 = GWAC_fitting(objx, objy, refx, refy, pointNum, ax, ay, cofNum,
                &objxrmsn, &objyrmsn, iter, rejsigma, statusstr);

        xrms2 = objxrmsn;
        yrms2 = objyrmsn;

#ifdef GWAC_TEST
        printFitDiff("reverse_orderndiff.txt", refx, refy, objx, objy, pointNum,
                ax, ay, cofNum, statusstr);
#endif
    }

    /*将拟合结果输出到文件*/
    retStatus6 = printBasicInfo(outfilename, "a", xrefmean, yrefmean, xmean, ymean, objxrms2,
            objyrms2, objxrmsn, objyrmsn, lineax, lineay, lineCof, ax, ay,
            cofNum, statusstr);

    if (refx != NULL) free(refx);
    if (refy != NULL) free(refy);
    if (objx != NULL) free(objx);
    if (objy != NULL) free(objy);
    if (ax != NULL) free(ax);
    if (ay != NULL) free(ay);
    if (lineax != NULL) free(lineax);
    if (lineay != NULL) free(lineay);
    if (afunc != NULL) free(afunc);

    CHECK_RETURN_SATUS(retStatus1);
    CHECK_RETURN_SATUS(retStatus2);
    CHECK_RETURN_SATUS(retStatus3);
    CHECK_RETURN_SATUS(retStatus4);
    CHECK_RETURN_SATUS(retStatus5);
    CHECK_RETURN_SATUS(retStatus6);

    return GWAC_SUCCESS;
}

/*******************************************************************************
 * 功能：二元多次多项式拟合，分别对x2和y2使用x1，y1进行拟合，进行iter次迭代，每次去掉拟合
 *      误差大于rejsigma倍的均方差的数据点，得到拟合系数ax，ay，以及最终的均方差
 * 
 **输入：
 *      x1, y1 输入的二元数据
 *      x2, y2 待拟合的数据
 *      dataNum 数据点的个数
 *      cofNum 拟合系数的个数，假设阶数为order，则其值为order*(order+1)/2
 *      iter 迭代次数
 *      rejsigma 去除偏差过大的数据点的阈值
 * 
 **输出
 *      ax, ay 分别对应拟合x2和y2的系数
 *      xrms, yrms 取出偏差过大的点后，最终ax和ay的均方差
 *          
 **返回值:
 *      0表示正确，其它值为错误码，
 *      蔡使用2001～2999，徐使用3001～3999，苑使用4001～4999，李使用5001～5999
 ******************************************************************************/
int GWAC_fitting(double *x1,
        double *y1,
        double *x2,
        double *y2,
        int dataNum,
        double *ax,
        double *ay,
        int cofNum,
        double *xrms,
        double *yrms,
        int iter,
        float rejsigma,
        char statusstr[]) {

    /*检查错误结果输出参数statusstr是否为空*/
    CHECK_STATUS_STR_IS_NULL(statusstr);
    /*检测输入参数是否为空*/
    CHECK_INPUT_IS_NULL(x1, "x1");
    CHECK_INPUT_IS_NULL(y1, "y1");
    CHECK_INPUT_IS_NULL(x2, "x2");
    CHECK_INPUT_IS_NULL(y2, "y2");
    CHECK_INPUT_IS_NULL(ax, "ax");
    CHECK_INPUT_IS_NULL(ay, "ay");
    CHECK_INPUT_IS_NULL(xrms, "xrms");
    CHECK_INPUT_IS_NULL(yrms, "yrms");

    double chisq;
    int i, j;
    bool hasError = false; //是否出现错误

    bool *rmvflg = (bool *) malloc((dataNum + 1) * sizeof (double));
    if (rmvflg == NULL) hasError = true;
    double *tx1 = (double *) malloc((dataNum + 1) * sizeof (double));
    if (tx1 == NULL) hasError = true;
    double *tx2 = (double *) malloc((dataNum + 1) * sizeof (double));
    if (tx2 == NULL) hasError = true;
    double *ty1 = (double *) malloc((dataNum + 1) * sizeof (double));
    if (ty1 == NULL) hasError = true;
    double *ty2 = (double *) malloc((dataNum + 1) * sizeof (double));
    if (ty2 == NULL) hasError = true;
    double *sig = (double *) malloc((dataNum + 1) * sizeof (double));
    if (sig == NULL) hasError = true;
    double *afunc = (double *) malloc((cofNum + 1) * sizeof (double));
    if (afunc == NULL) hasError = true;
    int *ia = (int *) malloc((cofNum + 1) * sizeof (int));
    if (ia == NULL) hasError = true;
    double **covar = (double **) malloc((cofNum + 1) * sizeof (double*)); //matrix cofNum * cofNum
    if (covar == NULL) hasError = true;

    for (i = 0; i <= cofNum; i++) {
        covar[i] = (double *) malloc((cofNum + 1) * sizeof (double));
        if (covar[i] == NULL) hasError = true;
    }

    if (!hasError) {
        for (i = 1; i <= dataNum; i++) {
            tx1[i] = x1[i];
            ty1[i] = y1[i];
            tx2[i] = x2[i];
            ty2[i] = y2[i];
            rmvflg[i] = false;
            sig[i] = 1;
        }
        for (i = 1; i <= cofNum; i++) {
            ia[i] = 1;
        }

        int leftNum = dataNum;
        int m, n;
        for (j = 0; j < iter; j++) {

            int tmpNum = leftNum;
            polynomialFit(tx1, ty1, tx2, sig, leftNum, ax, ia, cofNum, covar, &chisq, cofun);
            *xrms = sqrt(chisq / leftNum);
            polynomialFit(tx1, ty1, ty2, sig, leftNum, ay, ia, cofNum, covar, &chisq, cofun);
            *yrms = sqrt(chisq / leftNum);

            float xlimit = *xrms*rejsigma;
            float ylimit = *yrms*rejsigma;
#ifdef GWAC_TEST
            printf("iter=%d chisq=%f xrms=%f yrms=%f 3xrms=%f 3yrms=%f\n",
                    j + 1, chisq, *xrms, *yrms, ylimit, xlimit);
#endif

            for (m = 1; m <= dataNum; m++) {
                if (rmvflg[m])
                    continue;
                cofun(x1[m], y1[m], afunc, cofNum);
                double tmpx = 0.0;
                double tmpy = 0.0;
                for (n = 1; n <= cofNum; n++) {
                    tmpx += ax[n] * afunc[n];
                    tmpy += ay[n] * afunc[n];
                }

                float xdiff = fabs(x2[m] - tmpx);
                float ydiff = fabs(y2[m] - tmpy);
                if (xdiff > xlimit || ydiff > ylimit) {
                    //printf("%d \t%8.3f %8.3f %8.3f%%\n", m, y[m], tmpy, ydiff * 100);
                    rmvflg[m] = true;
                    tmpNum--;
                }
            }
#ifdef GWAC_TEST
            printf("remove %d point\n", leftNum - tmpNum);
#endif

            if (tmpNum >= leftNum) {
                break;
            }
            leftNum = tmpNum;

            n = 1;
            for (m = 1; m <= dataNum; m++) {
                if (!rmvflg[m]) {
                    tx1[n] = x1[m];
                    tx2[n] = x2[m];
                    ty1[n] = y1[m];
                    ty2[n] = y2[m];
                    n++;
                }
            }

#ifdef GWAC_TEST
            if (leftNum != (n - 1)) {
                printf("error leftNum=%d n-1=%d\n", leftNum, n - 1);
            }
#endif

        }
#ifdef GWAC_TEST
        printf("***********************************\n");
        printf("left points: %d\n", leftNum);
        printf("***********************************\n");
#endif
    }

    if (NULL != rmvflg) free(rmvflg);
    if (NULL != tx1) free(tx1);
    if (NULL != tx2) free(tx2);
    if (NULL != ty1) free(ty1);
    if (NULL != ty2) free(ty2);
    if (NULL != sig) free(sig);
    if (NULL != ia) free(ia);
    for (i = 0; i <= cofNum; i++) {
        if (NULL != covar[i]) free(covar[i]);
    }
    if (NULL != covar) free(covar);
    if (NULL != afunc) free(afunc);
    
    if(hasError){
        MALLOC_IS_NULL();
    }

    return GWAC_SUCCESS;
}

/*******************************************************************************
 * 功能：输出拟合数据的相关信息到指定的文件fName， 包括参考星表和目标星表的x，y坐标的均值和
 *      拟合后的均方差，平移，旋转，缩放，一阶拟合和高阶拟合的系数。
 * 
 **输入：
 *      fName 输出文件名
 *      otype 文件打开方式，只应该为"w"或"a"，"w"为清空原有内容再写入数据，"a"直接追加数据
 *      xrefmean, yrefmean, xmean, ymean 参考星表和目标星表的x，y坐标的均值
 *      cof2x, cof2y 二阶拟合系数
 *      cof2Num 二阶拟合系数的个数
 *      cofnx, cofny N阶拟合系数
 *      cofnNum N阶拟合系数的个数
 * 
 **输出
 *      fName 结果输出文件
 *      statusstr 错误返回值
 *          
 **返回值:
 *      0表示正确，其它值为错误码，
 *      蔡使用2001～2999，徐使用3001～3999，苑使用4001～4999，李使用5001～5999
 ******************************************************************************/
int printBasicInfo(const char *fName, char *otype, double xrefmean, double yrefmean,
        double xmean, double ymean, double xrms2, double yrms2, double xrmsn,
        double yrmsn, double *cof2x, double * cof2y, int cof2Num, double *cofnx,
        double * cofny, int cofnNum, char statusstr[]) {

    /*检查错误结果输出参数statusstr是否为空*/
    CHECK_STATUS_STR_IS_NULL(statusstr);
    /*检测输入参数是否为空*/
    CHECK_STRING_NULL_OR_EMPTY(fName, "fName");
    CHECK_INPUT_IS_NULL(otype, "otype");
    CHECK_INPUT_IS_NULL(cof2x, "cof2x");
    CHECK_INPUT_IS_NULL(cof2y, "cof2y");
    CHECK_INPUT_IS_NULL(cofnx, "cofnx");
    CHECK_INPUT_IS_NULL(cofny, "cofny");

    FILE *fp = fopen(fName, otype);
    CHECK_OPEN_FILE(fp, fName);

    int order2 = sqrt(2 * cof2Num);
    int ordern = sqrt(2 * cofnNum);

    double xshift, yshift, xmag, ymag, xrotation, yrotation, xrms, yrms;
    double a, b, c, d, e, f;
    a = cof2x[1];
    b = cof2x[2];
    c = cof2x[3];
    d = cof2y[1];
    e = cof2y[2];
    f = cof2y[3];
    xshift = a;
    yshift = d;
    xrotation = atan(e / b);
    yrotation = atan(c / f);
    xmag = b / cos(xrotation);
    ymag = c / sin(yrotation);

    if (cofnNum == 3) {
        xrms = xrms2;
        yrms = yrms2;
    } else {
        xrms = xrmsn;
        yrms = yrmsn;
    }

    fprintf(fp, "begin    first\n");
    fprintf(fp, "    xrefmean  %f\n", xrefmean);
    fprintf(fp, "    yrefmean  %f\n", yrefmean);
    fprintf(fp, "    xmean     %f\n", xmean);
    fprintf(fp, "    ymean     %f\n", ymean);
    fprintf(fp, "    geometry  general\n");
    fprintf(fp, "    function  polynomial\n");
    fprintf(fp, "    xshift    %f\n", xshift);
    fprintf(fp, "    yshift    %f\n", yshift);
    fprintf(fp, "    xmag      %f\n", xmag);
    fprintf(fp, "    ymag      %f\n", ymag);
    fprintf(fp, "    xrotation %f\n", xrotation);
    fprintf(fp, "    yrotation %f\n", yrotation);
    fprintf(fp, "    xrms      %f\n", xrms);
    fprintf(fp, "    yrms      %f\n", yrms);

    fprintf(fp, "    surface1    %d\n", 8 + cof2Num);
    fprintf(fp, "            3.    3.\n");
    fprintf(fp, "            %d.    %d.\n", order2, order2);
    fprintf(fp, "            %d.    %d.\n", order2, order2);
    fprintf(fp, "            0.    0.\n");
    fprintf(fp, "            0.0   0.0\n");
    fprintf(fp, "            0.0   0.0\n");
    fprintf(fp, "            0.0   0.0\n");
    fprintf(fp, "            0.0   0.0\n");

    int i;
#ifdef GWAC_TEST
    printf("fit coefficient: %d\n", cof2Num);
#endif
    for (i = 1; i <= cof2Num; i++) {
        fprintf(fp, "            %e\t%e\n", cof2x[i], cof2y[i]);
#ifdef GWAC_TEST
        printf("%e\t%e\n", cof2x[i], cof2y[i]);
#endif
    }

    fprintf(fp, "    surface2    %d\n", 8 + cofnNum);
    fprintf(fp, "            3.    3.\n");
    fprintf(fp, "            %d.    %d.\n", ordern, ordern);
    fprintf(fp, "            %d.    %d.\n", ordern, ordern);
    fprintf(fp, "            2.    2.\n");
    fprintf(fp, "            0.0   0.0\n");
    fprintf(fp, "            0.0   0.0\n");
    fprintf(fp, "            0.0   0.0\n");
    fprintf(fp, "            0.0   0.0\n");

#ifdef GWAC_TEST
    printf("fit coefficient: %d\n", cofnNum);
#endif

    for (i = 1; i <= cofnNum; i++) {
        if (cofnNum != cof2Num)
            fprintf(fp, "            %e\t%e\n", cofnx[i], cofny[i]);
        else
            fprintf(fp, "            %e\t%e\n", cof2x[i], cof2y[i]);
#ifdef GWAC_TEST
        if (cofnNum != cof2Num)
            printf("%e\t%e\n", cofnx[i], cofny[i]);
        else
            printf("%e\t%e\n", cof2x[i], cof2y[i]);
#endif
    }

    fclose(fp);
    return GWAC_SUCCESS;
}


/**
 * 测试函数
 * 输出目标星的实际坐标值和拟合坐标值的原始值，差值，等信息到文件
 * @param fName
 * @param refx
 * @param refy
 * @param inx
 * @param iny
 * @param dataNum
 * @param ax
 * @param ay
 * @param cofNum
 * @param statusstr
 * @return 
 */
int printFitDiff(const char *fName, double refx[], double refy[], double inx[], double iny[],
        int dataNum, double ax[], double ay[], int cofNum, char statusstr[]) {

    double xrms2 = 0.0;
    double yrms2 = 0.0;
    double *inxfit = (double *) malloc((dataNum + 1) * sizeof (double));
    double *inyfit = (double *) malloc((dataNum + 1) * sizeof (double));
    double *afunc = (double *) malloc((cofNum + 1) * sizeof (double));

    int i, j;
    for (i = 1; i <= dataNum; i++) {
        cofun(refx[i], refy[i], afunc, cofNum);
        double tmpx = 0.0;
        double tmpy = 0.0;
        for (j = 1; j <= cofNum; j++) {
            tmpx += ax[j] * afunc[j];
            tmpy += ay[j] * afunc[j];
        }
        inxfit[i] = tmpx;
        inyfit[i] = tmpy;
    }

    if (fName == NULL)
        return GWAC_ERROR;

    FILE *fp = fopen(fName, "w");
    if (fp == NULL)
        return GWAC_OPEN_FILE_ERROR;

    for (i = 1; i <= dataNum; i++) {
        double xdiff = inx[i] - inxfit[i];
        double ydiff = iny[i] - inyfit[i];
        xrms2 += xdiff*xdiff;
        yrms2 += ydiff*ydiff;
    }
    fprintf(fp, "xrms=%f\t", sqrt(xrms2 / dataNum));
    fprintf(fp, "yrms=%f\n", sqrt(yrms2 / dataNum));
    fprintf(fp, "\n");

    fprintf(fp, "objx\tobjy\tobjxfit\tobjyfit\txdiff\tydiff\n");
    for (i = 1; i <= dataNum; i++) {
        double xdiff = inx[i] - inxfit[i];
        double ydiff = iny[i] - inyfit[i];
        fprintf(fp, "%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f\n",
                inx[i], iny[i], inxfit[i], inyfit[i], xdiff, ydiff);
    }
    fclose(fp);
    free(inxfit);
    free(inyfit);
    free(afunc);

    return GWAC_SUCCESS;
}

/**
 * 测试函数
 * 从参数文件中读取参数
 * @param fName
 * @param xcof
 * @param ycof
 * @param cofNum
 * @return 
 */
int readCof(char *fName, double *xcof, double *ycof, int cofNum) {

    if (fName == NULL)
        return GWAC_ERROR;

    FILE *fp = fopen(fName, "r");
    if (fp == NULL)
        return GWAC_OPEN_FILE_ERROR;

    int i;
    for (i = 1; i <= cofNum; i++) {
        fscanf(fp, "%lf %lf", &xcof[i], &ycof[i]);
    }
    fclose(fp);

    return GWAC_SUCCESS;
}

/*******************************************************************************
 * 功能：函数头注释示例
 * 
 **输入：
 *      a 输入数据a
 *      b 输入数据b
 * 
 **输出
 *      c 输出数据c
 *          
 **返回值:
 *      0表示正确，其它值为错误码，
 *      蔡使用2001～2999，徐使用3001～3999，苑使用4001～4999，李使用5001～5999
 ******************************************************************************/
int testFunction(int a, int b, int &c) {
    return GWAC_SUCCESS;
}