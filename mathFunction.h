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


#endif	/* MATHFUNCTION_H */

