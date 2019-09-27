#ifndef OBS 
#define OBS


/* Data-related functions */
#include <gsl/gsl_matrix.h>
#include <stdio.h>

double CalcIfCensored(int censind, double x) ;


double CalcCovariateEffect(const double *param, const gsl_matrix * mtx_X, const int *inid);

double CalcCovariateEffect_mixture(const double *param, const gsl_matrix *mtx_X, const int *inid);

void SetNum(double *datat, int *datad, gsl_vector *vec_num, int *datan);

void SetNum2(double *datat, gsl_vector *vec_num2, int *datan);


#endif

