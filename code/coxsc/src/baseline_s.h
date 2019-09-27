#ifndef BASELINE_S
#define BASELINE_S

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


#include "math.h"


struct baseline_params_s
{
	const double *param;
	const gsl_matrix *mtx_t;
	gsl_vector *vec_baseline;
	const gsl_vector *vec_num;
	const double *datat;
	const int *datag;
	const int *datar;
	const int *datad;
	const double *datas;
	const int *datasidx;
	const gsl_matrix *mtx_X;
	const int *datan;
	const double *pr1;
	const double *pr3;
};


double CalculateBaseline_switch(const double value, void *params);

double CalculateBaseline_switch_mixture(const double value, void *params);

int SolveBaseline_switch(struct baseline_params_s *bparam);

void ConvertBaseline(gsl_vector *vec_baseline, gsl_vector *vec_num);

int SolveBaseline_switch_mixture(struct baseline_params_s *bparam);

double MixtureSumC(struct baseline_params_s *bparam, int xswidx, int swidx, int idx, int timeidx, int isterminal);

double CumBaseHaz(int time1, int time2, gsl_vector *vec_baseline);

double CumHaz(int time1, int time2, gsl_vector *vec_baseline, int switchidx, double psi1, double psi2);

int maxints(int v1, int v2);

int minints(int v1, int v2);

#endif
