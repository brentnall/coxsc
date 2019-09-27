const double MYZERO = 1.0e-20;
const double OUTOFBOUND = -100;
const double LIKPENALTY = -999;

/*data structure used by gsl optimiser function for pi */
struct ad_pidata
{
	double *param;
	gsl_vector *vec_event_time;
	gsl_vector_int *vec_event_type;
	int arm_n;
	double *pi;
	int lbnumswitch;
	int ubnumswitch;
};

/* Data structure used when fitting partial likelihood */
struct ad_optdata
{

	double *datat;
	int *datad;
	int *datag;
	int *datar;
	int *datan;
	int *datam;
	double *datas;
	double *datax;
	double unku;
	double *pr1;
	double *pr3;
	int *arm1_n;
	int *arm2_n;
	gsl_vector *vec_event_time;
	gsl_vector_int *vec_event_type;
	gsl_vector *vec_event_arm2_time;
	gsl_vector_int *vec_event_arm2_type;
	int lbnumswitch;
	int ubnumswitch;
};

// used for fitting partial likelihood when strata
struct ad_opdata_strata
{
	struct ad_optdata *optdata;
	int nstrata;
};

/* Prototypes */
/* Function to be called by R : partial likelihood (fitting pi)*/
void rint_pl(double *param, double *outlik, double *datat, int *datad, double *datax, int *datag, int *datar, int *datan, int *datam, double *datas, double *pr1, double *pitime, double *mypl, int *datasidx, double *pr3);

/* Testing full likelihood */
void testerfull(double *param, double *outlik, double *datat, int *datad, double *datax, int *datag, int *datar, int *datan, int *datam, double *datas, double *pr1, double *pitime, double *mypl, int *datasidx, double *pr3);

/* Function to be called by R: fitting model using Nelder-Mead */
//void testerpl(double *param, double *outlik, double *datat, int *datad, double *datax, int *datag, int *datar, int *datan, int *datam, double *datas, double *pr1, double *pitime, double *mypl );
void testerpl(double *param, double *outlik, double *datat, int *datad, double *datax, int *datag, int *datar, int *datan, int *datam, double *datas, double *pr1, double *pitime, double *mypl);

/* Estimates the partial likelihood given parameters, no fitting (other than pi) */
int est_partial_lik(double *param, double *datat, int *datad, int *datag, int *datar, int *datan, int *datam, double *datas, double *datax, double *unku, double *pr1, double *pr3, double *mypl);


/* This estimates the proportion of switchers at baseline and, optionally, the partial likelihood */
double partial_lik_all(double *param, double *datat, int *datad, int *datag, int *datar, int *datan, int *datam, double *datas, double *datax, double *unku, double *pr1, double *pr3, int *arm1_n, int *arm2_n, gsl_vector *vec_event_time, gsl_vector_int *vec_event_type, gsl_vector *vec_event_arm2_time, gsl_vector_int *vec_event_arm2_type, int lbnumswitch, int ubnumswitch, int plik);

/* Proportion of switchers amoungst those not had opportunity to switch yet, given number at start (unku) */
void propswitch(double *param, gsl_vector *vec_event_time, gsl_vector_int *vec_event_type, int *arm1_n, double *unku, double *pi1);

/* Log likelihood using struct object to pass the data */
double switch_loglik(double unku, void *indata);

/* 1D brent optimisation to find optimal number switchers at the start */
int find_startingu(double *unku, struct ad_pidata *indata);

/* This gets the number in each arm */
int setupdata_narms(double *datat, int *datad, int *datag, int *datar, int *datan, int *arm1_n, int *arm2_n);

/*Calculates the time and event data used to fit proportion of switchers */
int setupdata_tmev(double *datat, int *datad, int *datag, int *datar, int *datan, double *datas, gsl_vector * vec_event_time, gsl_vector_int * vec_event_type, gsl_vector * vec_event_arm2_time, gsl_vector_int * vec_event_arm2_type, int *arm1_n, int *ubnumswitch, int *lbnumswitch);

/* Fits the starting number of switcher in arm 1 (unku) */
int fitswitch(double *unku, struct ad_pidata *pidata);

/* Estimates the proportion of switchers in arms 1 and 2 given starting number (unku), switcher effect and other data such as survival time points */
int fitswitchprop(double *unku, struct ad_pidata *pidata, struct ad_pidata *pidata2, double *datat, int *datar, int *datan, double *pr1, double *pr3, double *allpi1, double *allpi3);

/* Initialise the hazards */
void inithazard(int *datan, int *datam, double *datax, double *param, gsl_matrix *mtx_X, gsl_matrix *mtx_t);

/* Calcs partial likelihood */
double partial_lik(double *param, double *datat, int *datad, int *datag, int *datar, int *datan, double *datas, double **pr, gsl_matrix *mtx_t);

/* GSL wrapper to partial_lik */
double partial_lik_all_gsl(const gsl_vector *vec_param, void *indata);

/* Fit the partial likelihood estimates */
int fit_partial_lik(double *param, double *datat, int *datad, int *datag, int *datar, int *datan, int *datam, double *datas, double *datax, double *unku, double *pr1, double *pr3, double *mypl);

/* Fit full likelihood */
double Lik3s_mixture(double *param, double *outlik, double *datat, int *datad, double *datax, int *datag, int *datar, int *datan, int *datam, double *datas, int *datasidx, double *pr1, double *pr3, gsl_vector *vec_outbase);

/* Estimate PL when strata */
int est_partial_lik_strata(double *param, double *datat, int *datad, int *datag, int *datar, int *datan, int *datam, double *datas, double *datax, double *unku, double *pr1, double *pr3, double *mypl, int *strata, int *nstrata);

void rint_pl_strata(double *param, double *outlik, double *datat, int *datad, double *datax, int *datag, int *datar, int *datan, int *datam, double *datas, double *pr1, double *pitime, double *mypl, int *strata, int *nstrata);

double partial_lik_all_strata_gsl(const gsl_vector *vec_param, void *indata);

int fit_partial_lik_strata(double *param, double *datat, int *datad, int *datag, int *datar, int *datan, int *datam, double *datas, double *datax, double *unku, double *pr1, double *pr3, double *mypl, int *strata, int *nstrata);

void testerpl_strata(double *param, double *outlik, double *datat, int *datad, double *datax, int *datag, int *datar, int *datan, int *datam, double *datas, double *pr1, double *pitime, double *mypl, int *strata, int *nstrata);

double Lik3s_mixture_strata(double *param, double *outlik, double *datat, int *datad, double *datax, int *datag, int *datar, int *datan, int *datam, double *datas, int *datasidx, double *pr1, double *pr3, gsl_vector *vec_outbase,  int *strata, int *nstrata);

void testerfull_strata(double *param, double *outlik, double *datat, int *datad, double *datax, int *datag, int *datar, int *datan, int *datam, double *datas, double *pr1, double *pitime, double *mypl, int *datasidx, double *pr3, int *strata, int *nstrata);

