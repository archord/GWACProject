
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gwac.h"
#include "nrutil.h"
#include "mathFunction.h"

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void covsrt(double **covar, int ma, int ia[], int mfit);
void gaussj(double **a, int n, double **b, int m);

//come from book of cai, has large residual

void tanPlaneToSphere2(double ra0, double dec0, double xi, double eta, double &ra, double &dec) {

    ra0 /= DEGREE_TO_RADIANS;
    dec0 /= DEGREE_TO_RADIANS;

    double a = cos(dec0) - eta * sin(dec0);
    double b = cos(ra - ra0);
    double c = eta * cos(dec0) + sin(dec0);
    double d = cos(dec0) - eta * sin(dec0);

    ra = (atan2(xi, a) + ra0) * DEGREE_TO_RADIANS;
    dec = (atan2(b * c, d)) * DEGREE_TO_RADIANS;
}

//the calculate formula come form slalib

void tanPlaneToSphere(double ra0, double dec0, double xi, double eta, double &ra, double &dec) {

    ra0 /= DEGREE_TO_RADIANS;
    dec0 /= DEGREE_TO_RADIANS;

    double sdec0 = sin(dec0);
    double cdec0 = cos(dec0);
    double denom = cdec0 - eta*sdec0;

    ra = (atan2(xi, denom) + ra0) * DEGREE_TO_RADIANS;
    dec = (atan2(sdec0 + eta*cdec0, sqrt(xi * xi + denom * denom))) * DEGREE_TO_RADIANS;
}

//come from book of cai

void tanSphereToPlane2(double ra0, double dec0, double ra, double dec, double &xi, double &eta) {

    ra0 /= DEGREE_TO_RADIANS;
    dec0 /= DEGREE_TO_RADIANS;
    ra /= DEGREE_TO_RADIANS;
    dec /= DEGREE_TO_RADIANS;

    double a = sin(ra - ra0) / tan(dec);
    double b = cos(dec0) * cos(ra - ra0) / tan(dec);
    double c = sin(dec0) * cos(ra - ra0) / tan(dec);

    xi = a / (sin(dec0) + b);
    eta = (cos(dec0) - c) / (sin(dec0) + b);
}

//the calculate formula come form slalib

int tanSphereToPlane(double ra0, double dec0, double ra, double dec, double &xi, double &eta) {

    ra0 /= DEGREE_TO_RADIANS;
    dec0 /= DEGREE_TO_RADIANS;
    ra /= DEGREE_TO_RADIANS;
    dec /= DEGREE_TO_RADIANS;

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

/*******************************************************************************
 * 
 * 
 ******************************************************************************/
void polynomialFit(double x1[], double x2[], double y[], double sig[], int dataNum, double a[], int ia[],
        int cofNum, double **covar, double *chisq, int (*cofun)(double, double, double [], int)) {

    int i, j, k, l, m, mfit = 0;
    double ym, wt, sum, sig2i, **beta;
    double *afunc;

    beta = dmatrix(1, cofNum, 1, 1);
    afunc = dvector(1, cofNum);
    for (j = 1; j <= cofNum; j++)
        if (ia[j]) mfit++;
    if (mfit == 0) nrerror("lfit: no parameters to be fitted");
    for (j = 1; j <= mfit; j++) {
        for (k = 1; k <= mfit; k++) covar[j][k] = 0.0;
        beta[j][1] = 0.0;
    }
    for (i = 1; i <= dataNum; i++) {
        (*cofun)(x1[i], x2[i], afunc, cofNum);
        ym = y[i];
        if (mfit < cofNum) {
            for (j = 1; j <= cofNum; j++)
                if (!ia[j]) ym -= a[j] * afunc[j];
        }
        sig2i = 1.0 / SQR(sig[i]);
        for (j = 0, l = 1; l <= cofNum; l++) {
            if (ia[l]) {
                wt = afunc[l] * sig2i;
                for (j++, k = 0, m = 1; m <= l; m++)
                    if (ia[m]) covar[j][++k] += wt * afunc[m];
                beta[j][1] += ym*wt;
            }
        }
    }
    for (j = 2; j <= mfit; j++)
        for (k = 1; k < j; k++)
            covar[k][j] = covar[j][k];
    gaussj(covar, mfit, beta, 1);
    for (j = 0, l = 1; l <= cofNum; l++)
        if (ia[l]) a[l] = beta[++j][1];
    *chisq = 0.0;
    for (i = 1; i <= dataNum; i++) {
        (*cofun)(x1[i], x2[i], afunc, cofNum);
        for (sum = 0.0, j = 1; j <= cofNum; j++) sum += a[j] * afunc[j];
        *chisq += SQR((y[i] - sum) / sig[i]);
        //        double tdiff = (y[i] - sum)/y[i];
        //        *chisq += ((tdiff*tdiff) /sig[i]);
    }
    covsrt(covar, cofNum, ia, mfit);
    free_dvector(afunc, 1, cofNum);
    free_dmatrix(beta, 1, cofNum, 1, 1);
}

/*******************************************************************************
 * 
 * 
 ******************************************************************************/
void gaussj(double **a, int n, double **b, int m) {
    int *indxc, *indxr, *ipiv;
    int i, icol, irow, j, k, l, ll;
    double big, dum, pivinv, temp;

    indxc = ivector(1, n);
    indxr = ivector(1, n);
    ipiv = ivector(1, n);
    for (j = 1; j <= n; j++) ipiv[j] = 0;
    for (i = 1; i <= n; i++) {
        big = 0.0;
        for (j = 1; j <= n; j++)
            if (ipiv[j] != 1)
                for (k = 1; k <= n; k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(a[j][k]) >= big) {
                            big = fabs(a[j][k]);
                            irow = j;
                            icol = k;
                        }
                    } else if (ipiv[k] > 1) nrerror("gaussj: Singular Matrix-1");
                }
        ++(ipiv[icol]);
        if (irow != icol) {
            for (l = 1; l <= n; l++) SWAP(a[irow][l], a[icol][l])
                for (l = 1; l <= m; l++) SWAP(b[irow][l], b[icol][l])
                }
        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix-2");
        pivinv = 1.0 / a[icol][icol];
        a[icol][icol] = 1.0;
        for (l = 1; l <= n; l++) a[icol][l] *= pivinv;
        for (l = 1; l <= m; l++) b[icol][l] *= pivinv;
        for (ll = 1; ll <= n; ll++)
            if (ll != icol) {
                dum = a[ll][icol];
                a[ll][icol] = 0.0;
                for (l = 1; l <= n; l++) a[ll][l] -= a[icol][l] * dum;
                for (l = 1; l <= m; l++) b[ll][l] -= b[icol][l] * dum;
            }
    }
    for (l = n; l >= 1; l--) {
        if (indxr[l] != indxc[l])
            for (k = 1; k <= n; k++)
                SWAP(a[k][indxr[l]], a[k][indxc[l]]);
    }
    free_ivector(ipiv, 1, n);
    free_ivector(indxr, 1, n);
    free_ivector(indxc, 1, n);
}

/*******************************************************************************
 * 
 * 
 ******************************************************************************/
void covsrt(double **covar, int cofNum, int ia[], int mfit) {
    int i, j, k;
    double temp;

    for (i = mfit + 1; i <= cofNum; i++)
        for (j = 1; j <= i; j++) covar[i][j] = covar[j][i] = 0.0;
    k = mfit;
    for (j = cofNum; j >= 1; j--) {
        if (ia[j]) {
            for (i = 1; i <= cofNum; i++) SWAP(covar[i][k], covar[i][j])
                for (i = 1; i <= cofNum; i++) SWAP(covar[k][i], covar[j][i])
                    k--;
        }
    }
}

/*******************************************************************************
 * 
 * 
 ******************************************************************************/
double meand(double *data, int len) {
    int i;
    double mean = 0.0;
    for (i = 1; i <= len; i++) {
        mean += (data[i] / len);
    }
    return mean;
}

/*******************************************************************************
 * 
 * 
 ******************************************************************************/
void quickSort(int min, int max, double a[]) {
    double key = a[min];
    int i = min;
    int j = max;
    //float temp;
    if (min >= max)
        return;
    while (i < j) {

        while ((i < j) && (key <= a[j])) {
            j--;
        }
        if (key > a[j]) {
            a[i] = a[j];
            a[j] = key;
            i++;
        }

        while ((i < j) && (key >= a[i])) {
            i++;
        }
        if (key < a[i]) {
            a[j] = a[i];
            a[i] = key;
            j--;
        }

    }
    quickSort(min, i - 1, a);
    quickSort(i + 1, max, a);
}

/*******************************************************************************
 * 
 * 
 ******************************************************************************/
double mediand(double array[], int len) {

    double median = 0.0;
    if (len % 2 == 0) {
        median = (array[len / 2 - 1] + array[len / 2]) / 2;
    } else {
        median = array[len / 2];
    }
    return median;
}

/*******************************************************************************
 * 功能：二元高阶最小二乘功能函数，计算各系数所对应项的值，系数个数cofNum=order*(order+1)
 *      其中order为最高阶次。
 * 
 **输入：
 *      x1 输入数据x
 *      x2 输入数据y
 *      cofNum 系数个数
 * 
 **输出
 *      afunc 输出系数
 *          
 **返回值:
 *      0表示正确，其它值为错误码，
 *      蔡使用2001～2999，徐使用3001～3999，苑使用4001～4999，李使用5001～5999
 ******************************************************************************/
int cofun(double x1, double x2, double *afunc, int cofNum) {

    if (NULL == afunc) {
        printf("File %s line %d, afunc is NULL\n", __FILE__, __LINE__);
        return GWAC_FUNCTION_INPUT_NULL;
    }

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

double getLgdCof(double x, int n) {

    if (n < 0) {
        return 0;
    }

    double p0 = 1;
    double p1 = x;
    if (n == 0) {
        return p0;
    } else if (n == 1) {
        return p1;
    } else {
        n = n - 1;
        double pn = getLgdCof(x, n);
        double pn_1 = getLgdCof(x, n - 1);
        return ((2 * n + 1) * x * pn - n * pn_1) / (n + 1);
    }
}

int cofun_Legendre(double x1, double x2, double *afunc, int cofNum) {

    if (NULL == afunc) {
        printf("File %s line %d, afunc is NULL\n", __FILE__, __LINE__);
        return GWAC_FUNCTION_INPUT_NULL;
    }

    int order = sqrt(2 * cofNum);
    int j, k;
    int i = 1;
    for (j = 0; j < order; j++) {
        for (k = 0; k < order - j; k++) {
            double x = getLgdCof(x1, k);
            double y = getLgdCof(x2, j);
            afunc[i++] = x * y;
        }
    }

    return GWAC_SUCCESS;
}