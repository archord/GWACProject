/* 
 * File:   mathFunction.h
 * Author: xy
 *
 * Created on September 23, 2013, 8:12 PM
 */

#ifndef MATHFUNCTION_H
#define	MATHFUNCTION_H


void polynomialFit(double x1[], double x2[], double y[], double sig[], int ndat, double a[], int ia[],
        int ma, double **covar, double *chisq, int (*funcs)(double, double, double [], int));
double meand(double *data, int len);
void quickSort(int min, int max, double a[]);
double mediand(double array[], int len);
void tanPlaneToSphere(double ra0, double dec0, double xi, double eta, double &ra, double &dec);
void tanPlaneToSphere2(double ra0, double dec0, double xi, double eta, double &ra, double &dec);
int tanSphereToPlane(double ra0, double dec0, double ra, double dec, double &xi, double &eta);
void tanSphereToPlane2(double ra0, double dec0, double ra, double dec, double &xi, double &eta);
int cofun(double x1, double x2, double *afunc, int cofNum);
int cofun_Legendre(double x1, double x2, double *afunc, int cofNum);

#endif	/* MATHFUNCTION_H */

