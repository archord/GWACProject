
#include "gwac.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mathFunction.h"

int cofun(double x1, double x2, double *afunc, int cofNum);
int printFitDiff(const char *fName, double refx[], double refy[], double inx[], double iny[],
        int dataNum, double ax[], double ay[], int cofNum, char statusstr[]);
int printCof(const char *fName, double *cofx, double * cofy, int cofNum, char statusstr[]);
int readCof(char *fName, double *xcof, double *ycof, int cofNum);

/*******************************************************************************
 **输入：
 *  matchpeervec   亮星匹配对
 *  order        拟合阶数
 *  iter         迭代次数
 *  rejsigma     数据抛弃sigma值
 **输出
 *  xrms         x方向拟合精度
 *  yrms         y方向拟合精度
 *  outfilename  输出拟合文件，格式同iraf geomap输出，具体描述见
 *  statusstr    返回值说明（不能超过128个字符）
 **返回值:         0表示正确，其它值为错误码，
 *  0蔡使用2001～2999，徐使用3001～3999，苑使用4001～4999，李使用5001～5999
 ******************************************************************************/
int Gwac_geomap(vector<ST_STARPEER> matchpeervec,
        unsigned int order,
        unsigned int iter,
        float rejsigma,
        float &xrms,
        float &yrms,
        const char outfilename[],
        char statusstr[]) {

    int pointNum = matchpeervec.size();
    int cofNum = (order + 1)*(order + 2) / 2;

    double *refx = (double *) malloc((pointNum + 1) * sizeof (double));
    double *refy = (double *) malloc((pointNum + 1) * sizeof (double));
    double *objx = (double *) malloc((pointNum + 1) * sizeof (double));
    double *objy = (double *) malloc((pointNum + 1) * sizeof (double));
    double *sig = (double *) malloc((pointNum + 1) * sizeof (double));
    double *ax = (double *) malloc((cofNum + 1) * sizeof (double));
    double *ay = (double *) malloc((cofNum + 1) * sizeof (double));
    int *ia = (int *) malloc((cofNum + 1) * sizeof (int));
    double **covar = (double **) malloc((cofNum + 1) * sizeof (double*)); //matrix cofNum * cofNum
    double chisq;

    int i, j;
    for (i = 0; i <= cofNum; i++) {
        covar[i] = (double *) malloc((cofNum + 1) * sizeof (double));
    }

    double *afunc = (double *) malloc((cofNum + 1) * sizeof (double));

    int validNum = pointNum;
    for (i = 0; i < validNum; i++) {
        ST_STAR tref = matchpeervec.at(i).ref;
        ST_STAR tobj = matchpeervec.at(i).obj;
        refx[i + 1] = tref.x;
        refy[i + 1] = tref.y;
        objx[i + 1] = tobj.x;
        objy[i + 1] = tobj.y;
        sig[i + 1] = 1;
    }
    for (i = 1; i <= cofNum; i++) {
        ia[i] = 1;
    }

    cofNum = 3;
    
    lfit(refx, refy, objx, sig, validNum, ax, ia, cofNum, covar, &chisq, cofun);
    xrms = sqrt(chisq / validNum);
    //printf("chisq=%f\txrms=%f\n", chisq, xrms);
    lfit(refx, refy, objy, sig, validNum, ay, ia, cofNum, covar, &chisq, cofun);
    yrms = sqrt(chisq / validNum);
    //printf("chisq=%f\tyrms=%f\n", chisq, yrms);
    printf("iter=%d\txrms=%f\tyrms=%f\n", 0, xrms, yrms);
    for(i=1; i<=cofNum; i++){
        printf("%d %f\t%f\n", i+1, ax[i], ay[i]);
    }

    int m, n;
    for (m = 1; m <= validNum; m++) {
        ST_STARPEER tmPeer = matchpeervec.at(m - 1);
        cofun(tmPeer.ref.x, tmPeer.ref.y, afunc, cofNum);
        double tmpx = 0.0;
        double tmpy = 0.0;
        for (n = 1; n <= cofNum; n++) {
            tmpx += ax[n] * afunc[n];
            tmpy += ay[n] * afunc[n];
        }
        //        tmPeer.obj.x = fabs(tmPeer.obj.x - tmpx);
        //        tmPeer.obj.y = fabs(tmPeer.obj.y - tmpy);
        matchpeervec.at(m - 1).obj.x = (tmPeer.obj.x - tmpx);
        matchpeervec.at(m - 1).obj.y = (tmPeer.obj.y - tmpy);
    }

    for (i = 0; i < validNum; i++) {
        ST_STAR tref = matchpeervec.at(i).ref;
        ST_STAR tobj = matchpeervec.at(i).obj;
        refx[i + 1] = tref.x;
        refy[i + 1] = tref.y;
        objx[i + 1] = tobj.x;
        objy[i + 1] = tobj.y;
        sig[i + 1] = 1;
    }


    cofNum = (order + 1)*(order + 2) / 2;
    for (j = 0; j < iter; j++) {
        int tmpNum = validNum;

        lfit(refx, refy, objx, sig, validNum, ax, ia, cofNum, covar, &chisq, cofun);
        xrms = sqrt(chisq / validNum);
        //printf("chisq=%f\txrms=%f\n", chisq, xrms);

        lfit(refx, refy, objy, sig, validNum, ay, ia, cofNum, covar, &chisq, cofun);
        yrms = sqrt(chisq / validNum);
        //printf("chisq=%f\tyrms=%f\n", chisq, yrms);

        printf("iter=%d\txrms=%f\tyrms=%f\n", j + 1, xrms, yrms);

        float xlimit = xrms*rejsigma;
        float ylimit = yrms*rejsigma;
        for (m = 1; m <= validNum; m++) {
            ST_STARPEER tmPeer = matchpeervec.at(m - 1);
            cofun(tmPeer.ref.x, tmPeer.ref.y, afunc, cofNum);
            double tmpx = 0.0;
            double tmpy = 0.0;
            for (n = 1; n <= cofNum; n++) {
                tmpx += ax[n] * afunc[n];
                tmpy += ay[n] * afunc[n];
            }
//            float xdiff = fabs(tmpx - tmPeer.obj.x);
//            float ydiff = fabs(tmpy - tmPeer.obj.y);
            float xdiff = fabs((tmpx - tmPeer.obj.x) / tmPeer.obj.x);
            float ydiff = fabs((tmpy - tmPeer.obj.y) / tmPeer.obj.y);
            if (xdiff > xlimit || ydiff > ylimit) {
                printf("%d \t%.3f \t%.3f \t%.3f \t%.3f %.3f \t%.3f\n", 
                        m, tmPeer.obj.x, tmpx, xdiff * 100 / tmPeer.obj.x,
                        tmPeer.obj.y, tmpy, ydiff * 100 / tmPeer.obj.y);
                matchpeervec.erase(matchpeervec.begin() + m - 1);
                validNum--;
                m--;
            }
        }

        if (tmpNum <= matchpeervec.size()) {
            break;
        }
        printf("remove %d point\n", tmpNum - validNum);
        validNum = matchpeervec.size();

        for (i = 0; i < validNum; i++) {
            ST_STAR tref = matchpeervec.at(i).ref;
            ST_STAR tobj = matchpeervec.at(i).obj;
            refx[i + 1] = tref.x;
            refy[i + 1] = tref.y;
            objx[i + 1] = tobj.x;
            objy[i + 1] = tobj.y;
            sig[i + 1] = 1;
        }
    }
    //    printf("left %d point\n", validNum);

    //    char *cofName = "coefficient/order6";
    //    readCof(cofName, ax, ay, cofNum);

    printCof(outfilename, ax, ay, cofNum, statusstr);
    char *fName = "fitDiff.txt";
    printFitDiff(fName, refx, refy, objx, objy, validNum, ax, ay, cofNum, statusstr);


    free(refx);
    free(refy);
    free(objx);
    free(objy);
    free(sig);
    free(ax);
    free(ay);
    free(ia);
    for (i = 0; i <= cofNum; i++) {
        free(covar[i]);
    }
    free(covar);
    free(afunc);
}

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

int printCof(const char *fName, double *cofx, double * cofy, int cofNum, char statusstr[]) {

    if (fName == NULL)
        return GWAC_ERROR;

    FILE *fp = fopen(fName, "w");
    if (fp == NULL)
        return GWAC_OPEN_FILE_ERROR;

    int order = sqrt(2 * cofNum);
    fprintf(fp, "	surface2	%d\n", 8 + cofNum);
    fprintf(fp, "			3.	3.\n");
    fprintf(fp, "			%d.	%d.\n", order, order);
    fprintf(fp, "			%d.	%d.\n", order, order);
    fprintf(fp, "			2.	2.\n");
    fprintf(fp, "\n");
    fprintf(fp, "\n");
    fprintf(fp, "\n");
    fprintf(fp, "\n");

    int i;
    printf("fit coefficient:\n");
    for (i = 1; i <= cofNum; i++) {
        fprintf(fp, "			%e\t%e\n", cofx[i], cofy[i]);
        printf("%e\t%e\n", cofx[i], cofy[i]);
    }

    fclose(fp);
    return GWAC_SUCCESS;
}

int printBasicInfo(FILE *fp, double x1[], double x2[], double y[], double sig[],
        int dataNum, double a[], int ia[], int cofNum, double **covar,
        double *chisq, void (*cofun)(double, double, double [], int),
        char statusstr[]) {
    double xrefmean, yrefmean, xmean, ymean, xshift, yshift;
    double xmag, ymag, xrotation, yrotation, xrms, yrms;

    return GWAC_SUCCESS;
}

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

    fprintf(fp, "inx\tinxfit\txdiff\txpercent(100%)\tiny\tinyfit\tydiff\typercent(100%)\n");
    for (i = 1; i <= dataNum; i++) {
        double xdiff = inx[i] - inxfit[i];
        double ydiff = iny[i] - inyfit[i];
        fprintf(fp, "%9.3f%9.3f%9.3f%9.3f%%%9.3f%9.3f%9.3f%9.3f%%\n",
                inx[i], inxfit[i], xdiff, xdiff * 100 / inx[i],
                iny[i], inyfit[i], ydiff, ydiff * 100 / iny[i]);
    }
    fclose(fp);
    free(inxfit);
    free(inyfit);
    free(afunc);

    return GWAC_SUCCESS;
}

int cofun(double x1, double x2, double *afunc, int cofNum) {

    int order = sqrt(2 * cofNum);
    int j, k;
    int i = 1;
    for (j = 0; j < order; j++) {
        for (k = 0; k < order - j; k++) {
            double x = pow(x1, k);
            double y = pow(x2, j);
            afunc[i++] = x * y;
        }
    }

    return GWAC_SUCCESS;
}

int cofun1(double x1, double x2, double *afunc, int cofNum) {

    //    int order = sqrt(2 * cofNum);
    //    int j, k;
    //    int i = 1;
    //    for (j = 0; j < order; j++) {
    //        for (k = 0; k < order - j; k++) {
    //            double x = pow(x1, k);
    //            double y = pow(x2, j);
    //            afunc[i++] = x * y;
    //        }
    //    }

    afunc[1] = 1.0;
    afunc[2] = x1;
    afunc[3] = x1*x1;
    afunc[4] = x2;
    afunc[5] = x1*x2;
    afunc[6] = x2*x2;

    return GWAC_SUCCESS;
}
