#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_minmax.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>


#include "obs.h"
#include "math.h"
#include "Rmain.h"
#include "baseline_s.h"

/* Functions */

/*test interface to r.partial likelihood.*/
void testerpl(double *param, double *outlik, double *datat, int *datad, double *datax, int *datag, int *datar, int *datan, int *datam, double *datas, double *pr1, double *pitime, double *mypl)
{
	
	int myout;

	double unku = 250;
	
	double ghostpar[3];

	ghostpar[0] = param[0];
	ghostpar[1] = 0;
	ghostpar[2] = param[1];
	ghostpar[3] = 0;

//Fit PL
	myout = fit_partial_lik(ghostpar, datat, datad, datag, datar, datan, datam, datas, datax, &unku, pr1, pitime, mypl);

	param[0] = ghostpar[0];	
	param[1] = ghostpar[2];	
}

void testerpl_strata(double *param, double *outlik, double *datat, int *datad, double *datax, int *datag, int *datar, int *datan, int *datam, double *datas, double *pr1, double *pitime, double *mypl, int *strata, int *nstrata)
{
	
	int myout;

	double unku = 250;
	
	double ghostpar[3];

	ghostpar[0] = param[0];
	ghostpar[1] = 0;
	ghostpar[2] = param[1];
	ghostpar[3] = 0;

//Fit PL
	myout = fit_partial_lik_strata(ghostpar, datat, datad, datag, datar, datan, datam, datas, datax, &unku, pr1, pitime, mypl, strata, nstrata);

	param[0] = ghostpar[0];	
	param[1] = ghostpar[2];	
}


/*test interface to r. mixture model from time zero, full likelihood.*/
void testerfull(double *param, double *outlik, double *datat, int *datad, double *datax, int *datag, int *datar, int *datan, int *datam, double *datas, double *pr1, double *pitime, double *mypl, int *datasidx, double *pr3)
{

	
	int myout;

	double unku = 250;
	
	double ghostpar[3];

	ghostpar[0] = param[0];
	ghostpar[1] = 0;
	ghostpar[2] = param[1];
	ghostpar[3] = 0;


	//FULL LIK
	gsl_vector * vec_baseline = gsl_vector_alloc (*datan);

	Lik3s_mixture(ghostpar, outlik, datat, datad, datax, datag, datar, datan, datam, datas, datasidx, pr1, pr3, vec_baseline);


	param[0] = ghostpar[0];	
	param[1] = ghostpar[2];	


	gsl_vector_free (vec_baseline);


}

/*test interface to r. mixture model from time zero, full likelihood.*/
void testerfull_strata(double *param, double *outlik, double *datat, int *datad, double *datax, int *datag, int *datar, int *datan, int *datam, double *datas, double *pr1, double *pitime, double *mypl, int *datasidx, double *pr3, int *strata, int *nstrata)
{

	
	int myout;

	double unku = 250;
	
	double ghostpar[3];

	ghostpar[0] = param[0];
	ghostpar[1] = 0;
	ghostpar[2] = param[1];
	ghostpar[3] = 0;

	gsl_vector * vec_baseline = gsl_vector_alloc (*datan);

	Lik3s_mixture_strata(ghostpar, outlik, datat, datad, datax, datag, datar, datan, datam, datas, datasidx, pr1, pr3, vec_baseline, strata, nstrata);


	param[0] = ghostpar[0];	
	param[1] = ghostpar[2];	


	gsl_vector_free (vec_baseline);


}

/*R interface to partial likelihood */
void rint_pl(double *param, double *outlik, double *datat, int *datad, double *datax, int *datag, int *datar, int *datan, int *datam, double *datas, double *pr1, double *pitime, double *mypl, int *datasidx, double *pr3)
{
	
	int myout;

	double unku = 250;
	
	double ghostpar[3];

	ghostpar[0] = param[0];
	ghostpar[1] = 0;
	ghostpar[2] = param[1];
	ghostpar[3] = 0;

//Obtain PL
	myout = est_partial_lik(ghostpar, datat, datad, datag, datar, datan, datam, datas, datax, &unku, pr1, pitime, mypl);

	
}


/*R interface to partial likelihood */
void rint_pl_strata(double *param, double *outlik, double *datat, int *datad, double *datax, int *datag, int *datar, int *datan, int *datam, double *datas, double *pr1, double *pitime, double *mypl, int *strata, int *nstrata)
{
	
	int myout;

	double unku = 250;
	
	double ghostpar[3];

	ghostpar[0] = param[0];
	ghostpar[1] = 0;
	ghostpar[2] = param[1];
	ghostpar[3] = 0;
//Obtain PL
	myout = est_partial_lik_strata(ghostpar, datat, datad, datag, datar, datan, datam, datas, datax, &unku, pr1, pitime, mypl, strata, nstrata);

	
}


/* Partial likelihood calculations, allow for stratified model
strata - array end / start of strata index length nstrata. Everything assumed orded by strata, then time.
*/
int est_partial_lik_strata(double *param, double *datat, int *datad, int *datag, int *datar, int *datan, int *datam, double *datas, double *datax, double *unku, double *pr1, double *pr3, double *mypl, int *strata, int *nstrata)
{

	int * arr_arm1_n = malloc ( ((*nstrata)) * sizeof(int));
	int * arr_arm2_n = malloc ( ((*nstrata)) * sizeof(int));
	
	struct ad_optdata * arr_optdata = malloc( ((*nstrata)) * sizeof(struct ad_optdata));
	struct ad_optdata optdata;
	int * tempdatan = malloc  ( ((*nstrata)) * sizeof(int) );
	int tempdatan1  = 0;
	int tempdatan2 = 0;
	
	gsl_vector * vec_event_time; 
	gsl_vector_int * vec_event_type; 
	gsl_vector * vec_event_arm2_time;
	gsl_vector_int * vec_event_arm2_type;

	int i; int j; int k;
	int mytest = 1;
	int lbnumswitch = 0, ubnumswitch = 0;
	double thispl;
	gsl_vector * vec_param = gsl_vector_alloc(2); //will need to change for more x's in future

	/* ============================================================
	  Set up data structures, for each strata 
	============================================================ */
 	gsl_vector_set(vec_param, 0 ,param[0]);
	gsl_vector_set(vec_param, 1 ,param[2]);

	for(i=0; i < *nstrata; i++)
	{

		tempdatan[i] = strata[i+1] - strata[i];	

		/* FOR EACH STRATA */
		//Numer in each arm, needed to dynamically allocate next structures
		mytest = setupdata_narms(&datat[strata[i]] , &datad[strata[i]], &datag[strata[i]] , &datar[strata[i]], &tempdatan[i], &arr_arm1_n[i], &arr_arm2_n[i]);


		// create GSL vectors in PL strata	
		vec_event_time = gsl_vector_alloc ((size_t) arr_arm1_n[i]);
		vec_event_type = gsl_vector_int_alloc ( (size_t) arr_arm1_n[i]);
		vec_event_arm2_time = gsl_vector_alloc ( (size_t) arr_arm2_n[i]);
		vec_event_arm2_type = gsl_vector_int_alloc ( (size_t) arr_arm2_n[i]);

		// set the vectors
		mytest = setupdata_tmev(&datat[strata[i]], &datad[strata[i]], &datag[strata[i]], &datar[strata[i]], &tempdatan[i], &datas[strata[i]], vec_event_time, vec_event_type, vec_event_arm2_time, vec_event_arm2_type, &arr_arm1_n[i], &ubnumswitch, &lbnumswitch);


		//INITIALISE struct	
               optdata = (struct ad_optdata) {&datat[strata[i]], &datad[strata[i]], &datag[strata[i]], &datar[strata[i]], &tempdatan[i], datam, &datas[strata[i]], &datax[strata[i]], *unku, &pr1[tempdatan1], &pr3[tempdatan2],  &arr_arm1_n[i],  &arr_arm2_n[i],  vec_event_time, vec_event_type, vec_event_arm2_time, vec_event_arm2_type, lbnumswitch, ubnumswitch};

		//add struct to array
		arr_optdata[i] = optdata;

		//update starting point for pi records
		tempdatan1 +=arr_arm1_n[i];
 		tempdatan2 +=arr_arm2_n[i];

	}

	/* ============================================================
	   Estimate PL 
	============================================================ */
	//set up data object
	struct ad_opdata_strata optdata_strata = {arr_optdata, *nstrata};
	thispl = 0;
	thispl = partial_lik_all_strata_gsl(vec_param, &optdata_strata);
	printf("This PL = %g \n", thispl);

	*mypl = thispl;
	//clean up
      

	for(i=0; i<(*nstrata); i++)
	{
		gsl_vector_free(arr_optdata[i].vec_event_time);
	        gsl_vector_free(arr_optdata[i].vec_event_arm2_time);
		gsl_vector_int_free(arr_optdata[i].vec_event_type);
		gsl_vector_int_free(arr_optdata[i].vec_event_arm2_type);

	}

	free(arr_optdata);
	free(arr_arm1_n);
	free(arr_arm2_n);
	free(tempdatan);

	gsl_vector_free(vec_param);

	return 1;

}

/* Partial likelihood including estimation of prop switchers through time */
int est_partial_lik(double *param, double *datat, int *datad, int *datag, int *datar, int *datan, int *datam, double *datas, double *datax, double *unku, double *pr1, double *pr3, double *mypl)
{

	int arm1_n = 0;
	int arm2_n = 0;
	int i;
	int mytest = 1;
	int lbnumswitch = 0, ubnumswitch = 0;
	double thispl;
	gsl_vector * vec_param = gsl_vector_alloc(2); //will need to change for more x's in future

	/* ============================================================
	  Set up data structures 
	============================================================ */
	
	gsl_vector_set(vec_param, 0 ,param[0]);
	gsl_vector_set(vec_param, 1 ,param[2]);

	//Numer in each arm, needed to dynamically allocate next structures
	mytest = setupdata_narms(datat, datad, datag, datar, datan, &arm1_n, &arm2_n);

	//******
	//ARM 1
        gsl_vector * vec_event_time = gsl_vector_alloc (arm1_n);

        gsl_vector_int * vec_event_type = gsl_vector_int_alloc (arm1_n);

	//ARM 2
        gsl_vector * vec_event_arm2_time = gsl_vector_alloc (arm2_n);

        gsl_vector_int * vec_event_arm2_type = gsl_vector_int_alloc (arm2_n);

	//Setup data
	mytest = setupdata_tmev(datat, datad, datag, datar, datan, datas, vec_event_time, vec_event_type, vec_event_arm2_time, vec_event_arm2_type, &arm1_n, &ubnumswitch, &lbnumswitch);

	/* ============================================================
	   Estimate PL 
	============================================================ */
        struct ad_optdata optdata = {datat, datad, datag, datar, datan, datam, datas, datax, *unku, pr1, pr3, &arm1_n, &arm2_n, vec_event_time, vec_event_type, vec_event_arm2_time, vec_event_arm2_type, lbnumswitch, ubnumswitch};

	thispl = partial_lik_all_gsl(vec_param, &optdata);
	printf("PL2 %g \n", thispl);

	*mypl = thispl;
	//clean up
	gsl_vector_int_free (vec_event_type);

	gsl_vector_free (vec_event_time);

        gsl_vector_free(vec_event_arm2_time);

	gsl_vector_free(vec_param);

        gsl_vector_int_free(vec_event_arm2_type);

	return 1;

}

/* Partial likelihood given pi in each arms. arrpi is 2 x number time points array of pis (col 0 = arm1, col 1 = arm2) */
double partial_lik(double *param, double *datat, int *datad, int *datag, int *datar, int *datan, double *datas, double **arrpi, gsl_matrix *mtx_t)
{
/*
Based on mtx_t 0 = c xs, 1 = c ys, 2 = t xs, 3= t ys; 3 is hazard for switchers in control post switch 
doesn't account for ties
*/
	int idX = 0;
	int idY = 0;
	double mypl = 0;
	double etai = 0;
	double etaj=0;
	double eta = 0;
	double mypi = 0;

	for(idX=0; idX < *datan; idX++)
	{
		if(datad[idX]==1)
		{
			if(datag[idX]==0) //switching status unknown, mixture
			{
				etai = (1-arrpi[ datar[idX]-1 ][idX]) * gsl_matrix_get(mtx_t, idX, 2*(datar[idX]-1)) ; //not a switcher

				etai += (arrpi[ datar[idX]-1 ][idX]) * gsl_matrix_get(mtx_t, idX, 2*(datar[idX]-1)+1) ; //a switcher
				
			} else //control arm and switching status known
			{
				etai = gsl_matrix_get(mtx_t, idX, (int) (3*((double) (datag[idX]-1))) ) ; 

			}

				eta = etai;

			for(idY = idX+1; idY< *datan; idY++)	
			{

				if((datas[idY] > datat[idX]) | (datag[idY] == 0)) //switching status unknown, mixture
				{
					etaj = (1- arrpi[ datar[idY]-1 ][idX]) * gsl_matrix_get(mtx_t, idY, 2*(datar[idY]-1)) ; //not a switcher
		
					etaj += (arrpi[ datar[idY]-1 ][idX]) * gsl_matrix_get(mtx_t, idY, 2*(datar[idY]-1)+1) ; //a switchera
			
				} else //control arm and switching status known
				{
					etaj = gsl_matrix_get(mtx_t, idY, (int) ( 3 * (((double)  datag[idY])-1)));
				}
				
				eta += etaj;
			}
			
			mypl += log(etai/eta);
		}
		
	}

	return mypl;	

}

/* Number of individuals in each arm */
int setupdata_narms(double *datat, int *datad, int *datag, int *datar, int *datan, int *arm1_n, int *arm2_n)
{
	int i;
	int arm2n = 0;	

	// number in each arm
	for (i = 0; i < *datan; i++)
	{
	//	if(datad[i]==1)
		{
			if((datar[i]==2))
			{
				arm2n++;
			}
		}
	}

	*arm2_n = arm2n;

	*arm1_n = *datan - *arm2_n;

	
	return 1;

}

/* This returns some of the data structures used to fit proportion of switchers at start in arm1 */
int setupdata_tmev(double *datat, int *datad, int *datag, int *datar, int *datan, double *datas, gsl_vector * vec_event_time, gsl_vector_int * vec_event_type, gsl_vector * vec_event_arm2_time, gsl_vector_int * vec_event_arm2_type, int *arm1_n, int *ubnumswitch, int *lbnumswitch)
{

	int i;
	int mycounter = 0;
	int mycounter3 = 0;

	//vector of event times to calculate proportion of switchers
	//ARM 1
	gsl_permutation * perm = gsl_permutation_alloc(*arm1_n);

	//ARM 2
	mycounter = 0;

	mycounter3 = 0;

	*ubnumswitch = *arm1_n;
	*lbnumswitch = 0;

	for (i = 0; i < *datan; i++)
	{
		if(datar[i]==1) //arm 1 only
		{
			if( (datad[i]==0) & (datag[i]==0) ) //censoring before switch opp
			{
				gsl_vector_int_set(vec_event_type, mycounter,3 ); //event ( 0-2) or not (3) before switch
			}
			else
			{
				gsl_vector_int_set(vec_event_type, mycounter, datag[i] ); //event ( 0-2) or not (3) before switch
			}

			if( datat[i] >= datas[i] )
			{
				gsl_vector_set(vec_event_time, mycounter, datas[i]);
				
				if(datag[i]==2) 
				{
					(*lbnumswitch) ++; ////////////
				}
				else if(datag[i]==1)
				{
					(*ubnumswitch)--;
				}
			} 
			else
			{
				gsl_vector_set(vec_event_time, mycounter, datat[i]);
			}
		
			mycounter++;
		}
		else 	//ARM 2
		{
			gsl_vector_set(vec_event_arm2_time, mycounter3, datat[i]);

			if( (datad[i]==0) ) //censored
			{
				gsl_vector_int_set(vec_event_arm2_type, mycounter3, 3);
 			}
			else
			{
				gsl_vector_int_set(vec_event_arm2_type, mycounter3, 0); 
			}
			mycounter3++;
		}
	}


	// sort time points
	gsl_sort_vector_index (perm, vec_event_time);

	gsl_permute_vector(perm, vec_event_time);

	gsl_permute_vector_int(perm, vec_event_type);

  //clean up

       gsl_permutation_free (perm);

	return 1;

}

/* This fits the number of switchers at the start, using arm1 data */
int fitswitch(double *unku, struct ad_pidata *pidata)
{
        double lbll = 0, ubll = 0, mll = 0;
        double guesspi = 0.5;
	double thislik=0;
	double optlik = 0;
	int idx;

	//optimise

	guesspi = ((double) (pidata->lbnumswitch) / ((double) pidata->lbnumswitch + pidata->ubnumswitch));	

	*unku = guesspi * ((double) pidata->arm_n);

	//****************
	// 1D solver assumes that log likelihood at unku is less than at ends. Following loop makes sure of this before calling it.
	lbll = switch_loglik(pidata->lbnumswitch, pidata);

	ubll =  switch_loglik(pidata->ubnumswitch, pidata);

	mll = switch_loglik(*unku, pidata);

	if(mll>ubll | mll>lbll)
	{

		//find lowest lik between lb and ub and use that
		if(lbll<ubll)
		{
			optlik = lbll;
			*unku = pidata->lbnumswitch;
		}
		else
		{
                        optlik = ubll;
                        *unku = pidata->ubnumswitch;
		}
		for(idx=pidata->lbnumswitch+1; idx<=pidata->ubnumswitch-1; idx++)
		{
			thislik = switch_loglik(idx, pidata);
			if(thislik < optlik)
			{
				optlik = thislik;
				*unku = idx;
				break;
			}
		}
		if(idx>=pidata->ubnumswitch-1)
		{
			printf("NOTE: COULD NOT RUN SWITCHER PROPORTION ALGORITHM, ESTIMATE OF STARTING POINT MIGHT NOT BE VERY GOOD\n");
			return -1;
		}
	}

	//****************
	/* Call the algorithm */
	if(find_startingu(unku, pidata) == -1)
	{
		return -1;
	}
	else
	{
		return 1;
	}
	


}
 
/* This gives the proportion of switchers at each survival time point in data passed */
int fitswitchprop(double *unku, struct ad_pidata *pidata, struct ad_pidata *pidata2, double *datat, int *datar, int *datan, double *pr1, double *pr3, double *allpi1, double *allpi3)
{
	//save results for event times (not survival)
	propswitch(pidata->param, pidata->vec_event_time, pidata->vec_event_type, &(pidata->arm_n), unku, pidata->pi);

	propswitch(pidata2->param, pidata2->vec_event_time, pidata2->vec_event_type, &(pidata2->arm_n), unku, pidata2->pi);

	// For original times (survival)
	int i = 0;
	int mycounter = 0; //switch event time arm 1
	int mycounter2 = 0; //survival event time arm 1
	int mycounter4 = 0;

	for (i = 0; i < *datan; i++)
	{
		if(datar[i]==1) //arm 1 
		{
			if(mycounter < pidata->arm_n) 
			{
				while( (gsl_vector_get(pidata->vec_event_time, mycounter) < datat[i]) & (mycounter < ((pidata->arm_n)-1)) )
				{
					mycounter++;
				} 
			}
				

			pr1[mycounter2] = pidata->pi[mycounter];


			mycounter2++;

		}
		else //arm 2
		{

			pr3[mycounter4] = pidata2->pi[mycounter4];

			mycounter4++;

		}

	}

	/* Now get pi values as originally ordered */

	double lastpi1;
	double lastpi3;

	lastpi1 = pr1[0];
	lastpi3 = pr3[0];

	mycounter = 0;
	mycounter2 = 0;

	for(i = 0; i< *datan; i++)
	{
		if(datar[i]==1)
		{
			lastpi1 = pr1[mycounter];

			mycounter ++;
		}
		else
		{
			lastpi3 = pr3[mycounter2];

			mycounter2++;
		}		

		allpi1[i] = lastpi1;

		allpi3[i] = lastpi3;

		
	}

  

	return 1;

}

/* Initialise the hazard and covariate matricies */
void inithazard(int *datan, int *datam, double *datax, double *param, gsl_matrix *mtx_X, gsl_matrix *mtx_t)
{
	int i = 0;
	int j = 0;

        /* 2.1 put covariates into GSL matrix to help avoid programming mistakes */
        for (i = 0; i < *datan; i++)
                 for (j = 0; j < *datam; j++)
                        gsl_matrix_set (mtx_X, i, j, datax[i + (j * (*datan))]);

        /* 2.2 Possible hazards for each i=1,..,n*/

        for (i = 0; i < *datan; i++)
        {
                double coveffect;

                coveffect = CalcCovariateEffect(param, mtx_X, &i);

                gsl_matrix_set(mtx_t, i, 0, exp(coveffect)); //control, don't switch pre switch

                gsl_matrix_set(mtx_t, i, 1, exp(param[0] + coveffect)); //control, switcher, preswitch

                gsl_matrix_set(mtx_t, i, 2, exp(param[2] + coveffect)); // treatment, not switcher

                gsl_matrix_set(mtx_t, i, 3, exp(param[0] + param[2] + coveffect)); //control, switcher, postswitch, treatment switcher

        }


}

/* Partial likelihood including estimation of prop switchers through time */
int fit_partial_lik(double *param, double *datat, int *datad, int *datag, int *datar, int *datan, int *datam, double *datas, double *datax, double *unku, double *pr1, double *pr3, double *mypl)
{

	int arm1_n = 0;
	int arm2_n = 0;
	int i;
	int mytest = 1;
	int lbnumswitch = 0, ubnumswitch = 0;
	double thispl;
	gsl_vector * vec_param = gsl_vector_alloc(2); //will need to change for more x's in future

	/* ============================================================
	 1. Set up data structures 
	============================================================ */
	
//	for(i=0; i<4; i++)
//	{
	gsl_vector_set(vec_param, 0 ,param[0]);
	gsl_vector_set(vec_param, 1 ,param[2]);

//	}

	//Numer in each arm, needed to dynamically allocate next structures
	mytest = setupdata_narms(datat, datad, datag, datar, datan, &arm1_n, &arm2_n);

	//******
	//ARM 1
        gsl_vector * vec_event_time = gsl_vector_alloc (arm1_n);

        gsl_vector_int * vec_event_type = gsl_vector_int_alloc (arm1_n);

	//ARM 2
        gsl_vector * vec_event_arm2_time = gsl_vector_alloc (arm2_n);

        gsl_vector_int * vec_event_arm2_type = gsl_vector_int_alloc (arm2_n);

	//Setup data
	mytest = setupdata_tmev(datat, datad, datag, datar, datan, datas, vec_event_time, vec_event_type, vec_event_arm2_time, vec_event_arm2_type, &arm1_n, &ubnumswitch, &lbnumswitch);

	/* ============================================================
	 4. Fit 
	============================================================ */
        struct ad_optdata optdata = {datat, datad, datag, datar, datan, datam, datas, datax, *unku, pr1, pr3, &arm1_n, &arm2_n, vec_event_time, vec_event_type, vec_event_arm2_time, vec_event_arm2_type, lbnumswitch, ubnumswitch};


	thispl = partial_lik_all_gsl(vec_param, &optdata);

	*mypl = thispl;

	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
       	gsl_vector *ss;
       	gsl_multimin_function minex_func;
     
       	size_t iter = 0;
       	int status;
       	double size;
     
       //Set initial step sizes to 1 
       	ss = gsl_vector_alloc (2);
       	gsl_vector_set_all (ss, 0.1);
     

       //Initialize method and iterate 
       	minex_func.n = 2;
       	minex_func.f = &partial_lik_all_gsl;
       	minex_func.params = &optdata;

       	s = gsl_multimin_fminimizer_alloc (T, 2);
	gsl_multimin_fminimizer_set (s, &minex_func, vec_param, ss);    	
       	do
	{
		iter++;
        	status = gsl_multimin_fminimizer_iterate(s);
           	
	        if (status) 
	             break;
     
        	size = gsl_multimin_fminimizer_size (s);
	        status = gsl_multimin_test_size (size, 1e-2);
     
        	if (status == GSL_SUCCESS)
	        {
               		printf ("converged to minimum at\n");
  	        }
     
        }
 	while (status == GSL_CONTINUE && iter < 100);

	//printf("Fit: %g, %g; iter %u \n", exp(gsl_vector_get (s->x, 0)), exp(gsl_vector_get (s->x, 1)), iter);

	param[0] = gsl_vector_get (s->x, 0);

	param[2] = gsl_vector_get (s->x, 1);

       	

	/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	//gsl_vector_free(x);
       	gsl_vector_free(ss);
       	gsl_multimin_fminimizer_free (s);

	//clean up
	gsl_vector_int_free (vec_event_type);

	gsl_vector_free (vec_event_time);

        gsl_vector_free(vec_event_arm2_time);

	gsl_vector_free(vec_param);

        gsl_vector_int_free(vec_event_arm2_type);

	return 1;

}

/* Partial likelihood including estimation of prop switchers through time */
int fit_partial_lik_strata(double *param, double *datat, int *datad, int *datag, int *datar, int *datan, int *datam, double *datas, double *datax, double *unku, double *pr1, double *pr3, double *mypl, int *strata, int *nstrata)
{

	int * arr_arm1_n = malloc ( ((*nstrata)) * sizeof(int));
	int * arr_arm2_n = malloc ( ((*nstrata)) * sizeof(int));
	
	struct ad_optdata * arr_optdata = malloc( ((*nstrata)) * sizeof(struct ad_optdata));
	struct ad_optdata optdata;
	int * tempdatan = malloc  ( ((*nstrata)) * sizeof(int) );
	int tempdatan1  = 0;
	int tempdatan2 = 0;
	
	gsl_vector * vec_event_time; 
	gsl_vector_int * vec_event_type; 
	gsl_vector * vec_event_arm2_time;
	gsl_vector_int * vec_event_arm2_type;

	int i; int j; int k;
	int mytest = 1;
	int lbnumswitch = 0, ubnumswitch = 0;
	double thispl;
	gsl_vector * vec_param = gsl_vector_alloc(2); //will need to change for more x's in future

	/* ============================================================
	  Set up data structures, for each strata 
	============================================================ */
 	gsl_vector_set(vec_param, 0 ,param[0]);
	gsl_vector_set(vec_param, 1 ,param[2]);

	for(i=0; i < *nstrata; i++)
	{

		tempdatan[i] = strata[i+1] - strata[i];	

		/* FOR EACH STRATA */
		//Numer in each arm, needed to dynamically allocate next structures
		mytest = setupdata_narms(&datat[strata[i]] , &datad[strata[i]], &datag[strata[i]] , &datar[strata[i]], &tempdatan[i], &arr_arm1_n[i], &arr_arm2_n[i]);

		// create GSL vectors in PL strata	
		vec_event_time = gsl_vector_alloc ((size_t) arr_arm1_n[i]);
		vec_event_type = gsl_vector_int_alloc ( (size_t) arr_arm1_n[i]);
		vec_event_arm2_time = gsl_vector_alloc ( (size_t) arr_arm2_n[i]);
		vec_event_arm2_type = gsl_vector_int_alloc ( (size_t) arr_arm2_n[i]);

		// set the vectors
		mytest = setupdata_tmev(&datat[strata[i]], &datad[strata[i]], &datag[strata[i]], &datar[strata[i]], &tempdatan[i], &datas[strata[i]], vec_event_time, vec_event_type, vec_event_arm2_time, vec_event_arm2_type, &arr_arm1_n[i], &ubnumswitch, &lbnumswitch);


		//INITIALISE struct	
               optdata = (struct ad_optdata) {&datat[strata[i]], &datad[strata[i]], &datag[strata[i]], &datar[strata[i]], &tempdatan[i], datam, &datas[strata[i]], &datax[strata[i]], *unku, &pr1[tempdatan1], &pr3[tempdatan2],  &arr_arm1_n[i],  &arr_arm2_n[i],  vec_event_time, vec_event_type, vec_event_arm2_time, vec_event_arm2_type, lbnumswitch, ubnumswitch};

		//add struct to array
		arr_optdata[i] = optdata;

		//update starting point for pi records
		tempdatan1 +=arr_arm1_n[i];
 		tempdatan2 +=arr_arm2_n[i];

	}

	/* ============================================================
	   Estimate PL 
	============================================================ */
	//set up data object
	struct ad_opdata_strata optdata_strata = {arr_optdata, *nstrata};
	thispl = 0;
	thispl = partial_lik_all_strata_gsl(vec_param, &optdata_strata);
	

	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
       	gsl_vector *ss;
       	gsl_multimin_function minex_func;
     
       	size_t iter = 0;
       	int status;
       	double size;
         
       //Set initial step sizes to 1 
       	ss = gsl_vector_alloc (2);
       	gsl_vector_set_all (ss, 0.1);
     

       //Initialize method and iterate 
       	minex_func.n = 2;
       	minex_func.f = &partial_lik_all_strata_gsl;
       	minex_func.params = &optdata_strata;

       	s = gsl_multimin_fminimizer_alloc (T, 2);
	
	gsl_multimin_fminimizer_set (s, &minex_func, vec_param, ss);    	
       	
	do
	{
		iter++;
        	status = gsl_multimin_fminimizer_iterate(s);
           	
	        if (status) 
	             break;
     
        	size = gsl_multimin_fminimizer_size (s);
	        status = gsl_multimin_test_size (size, 1e-2);
     
        	if (status == GSL_SUCCESS)
	        {
               		printf ("converged to minimum at\n");
  	        }
     
//           	printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", 
 //                  iter,
  //                 gsl_vector_get (s->x, 0), 
   //                gsl_vector_get (s->x, 1), 
    //               s->fval, size);
        }
 	while (status == GSL_CONTINUE && iter < 100);

	//printf("Fit: %g, %g; iter %u \n", exp(gsl_vector_get (s->x, 0)), exp(gsl_vector_get (s->x, 1)), iter);

	param[0] = gsl_vector_get (s->x, 0);

	param[2] = gsl_vector_get (s->x, 1);

       	

	/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	//gsl_vector_free(x);
       	gsl_vector_free(ss);
       	gsl_multimin_fminimizer_free (s);

	//clean up
        gsl_vector_free(vec_param);

	for(i=0; i<(*nstrata); i++)
	{
		gsl_vector_free(arr_optdata[i].vec_event_time);
	        gsl_vector_free(arr_optdata[i].vec_event_arm2_time);
		gsl_vector_int_free(arr_optdata[i].vec_event_type);
		gsl_vector_int_free(arr_optdata[i].vec_event_arm2_type);

	}

	free(arr_optdata);
	free(arr_arm1_n);
	free(arr_arm2_n);
	free(tempdatan);

	//end

	return 1;

}


/* Wrapper to the partial likelihood function, using form required by GSL optimisation */
double partial_lik_all_gsl(const gsl_vector *vec_param, void *indata)
{
	double thispl=0;	

	struct ad_optdata *od = (struct ad_optdata *) indata;

	double param[4];

	int i = 0;
	param[0] = gsl_vector_get(vec_param, 0);
	param[1] = 0;
	param[2] = gsl_vector_get(vec_param, 1);
	param[3] = 0;

	thispl= partial_lik_all(param, od->datat, od->datad, od->datag, od->datar, od->datan, od->datam, od->datas, od->datax, &(od->unku), od->pr1, od->pr3, od->arm1_n, od->arm2_n, od->vec_event_time, od->vec_event_type, od->vec_event_arm2_time, od->vec_event_arm2_type, od->lbnumswitch, od->ubnumswitch, (int) 1);

	return -thispl;

}

/* Wrapper to the partial likelihood function allowing for strata, using form required by GSL optimisation */
double partial_lik_all_strata_gsl(const gsl_vector *vec_param, void *indata)
{
	double thispl=0;	

	struct ad_opdata_strata *optdata_strata = (struct ad_opdata_strata *) indata;
	struct ad_optdata *arr_optdata = (struct ad_optdata *) optdata_strata->optdata; // first part is pointer to array of optdata
	struct ad_optdata * od;
	int nstrata = (int) optdata_strata->nstrata; // second part is number of strata

	double param[4];

	param[0] = gsl_vector_get(vec_param, 0);
	param[1] = 0;
	param[2] = gsl_vector_get(vec_param, 1);
	param[3] = 0;

	int i = 0;
	
	for(i=0; i<nstrata; i++)
	{
		od = &arr_optdata[i];

		thispl += partial_lik_all(param, od->datat, od->datad, od->datag, od->datar, od->datan, od->datam, od->datas, od->datax, &(od->unku), od->pr1, od->pr3, od->arm1_n, od->arm2_n, od->vec_event_time, od->vec_event_type, od->vec_event_arm2_time, od->vec_event_arm2_type, od->lbnumswitch, od->ubnumswitch, (int) 1);

	//	printf("PL %g \n", thispl);
	}

	return -thispl;

}


/* Partial likelihood including estimation of prop switchers through time */
double partial_lik_all(double *param, double *datat, int *datad, int *datag, int *datar, int *datan, int *datam, double *datas, double *datax, double *unku, double *pr1, double *pr3, int *arm1_n, int *arm2_n, gsl_vector *vec_event_time, gsl_vector_int *vec_event_type, gsl_vector *vec_event_arm2_time, gsl_vector_int *vec_event_arm2_type, int lbnumswitch, int ubnumswitch, int plik)
{




	int i;
	int mytest = 1;
	double * allpi1 = malloc ( (*datan) * sizeof(double));
	double * allpi3 = malloc ( (*datan) * sizeof(double));
	double *myallpi[2];
	double mypl;

	
	/* ============================================================
	 1. Set up data structures for proportion of switchers
	============================================================ */

	double * pi1 = malloc ((*arm1_n) * sizeof(double)); //proportion switcher at event times (not survival)

	double * pi3 = malloc ((*arm2_n) * sizeof(double));

	//Setup data
	mytest = setupdata_tmev(datat, datad, datag, datar, datan, datas, vec_event_time, vec_event_type, vec_event_arm2_time, vec_event_arm2_type, arm1_n, &ubnumswitch, &lbnumswitch);


        struct ad_pidata pidata= {param, vec_event_time, vec_event_type, *arm1_n, pi1, lbnumswitch, ubnumswitch};

        struct ad_pidata pidata2= {param, vec_event_arm2_time, vec_event_arm2_type, *arm2_n, pi3, lbnumswitch, ubnumswitch};

	/* ============================================================
	 2. Set up data structures for fitting
	============================================================ */

	gsl_matrix * mtx_t = gsl_matrix_alloc (*datan, 4); //

        gsl_matrix * mtx_X = gsl_matrix_alloc (*datan, *datam);

	inithazard(datan, datam, datax, param, mtx_X, mtx_t);

	/* ============================================================
	 3. Estimate proportion of switchers
	============================================================ */
	
	//3.1 find starting point
	if(fitswitch(unku, &pidata)==-1) //current parameters not compatible
	{
                mypl = LIKPENALTY;
	}
	else
	{
		mypl = 0;
	}


	//3.2 get proportion switcher estimates
	mytest = fitswitchprop(unku, &pidata, &pidata2, datat, datar, datan, pr1, pr3, allpi1, allpi3);

	//3.3 point to proportion through time in 2-d pointer array
	myallpi[0] = allpi1;
	myallpi[1] = allpi3;

	/* ============================================================
	 4. Partial likelihood
	============================================================ */

	if(plik ==1)
		{
			mypl += partial_lik(param, datat, datad, datag, datar, datan, datas, myallpi, mtx_t);
		}
	else
		{
			mypl = 0;
		}
	

	exit_x:

	//clean up
	gsl_matrix_free(mtx_t);

	gsl_matrix_free(mtx_X);

	free(pi1);

	free(pi3);

	free(allpi1);

	free(allpi3);

	return mypl;

}

/* 1D brent optimisation to find optimal number switchers at the start */
int find_startingu(double *unku, struct ad_pidata *indata)
{
	gsl_set_error_handler_off(); //large likelihood penalty for such an outcome

	struct ad_pidata *pidata = (struct ad_pidata *) indata;

       	int status;
       	int iter = 0, max_iter = 100;
       	const gsl_min_fminimizer_type *T;
       	gsl_min_fminimizer *s;

       	double a;
	a = (double) pidata->lbnumswitch;
	double b;
	b = (double) pidata->ubnumswitch;

	//printf("Lower %g Upper %g", a, b);

	double m;
	m = *unku;
//	m = (int) (a+b)/2;

       gsl_function F;
       F.function = &switch_loglik;
       F.params = indata;
     
       T = gsl_min_fminimizer_brent;
       s = gsl_min_fminimizer_alloc (T);
//       printf ("using %s method\n", gsl_min_fminimizer_name (s));

	status = gsl_min_fminimizer_set (s, &F, m, a, b);

	if(status==4)
	{
		printf("NOTE: COULD NOT RUN SWITCHER PROPORTION ALGORITHM, ESTIMATE OF STARTING POINT MIGHT NOT BE VERY GOOD\n");
		printf("a %g, b%g, m %g, la %g, lb %g, lm %g\n", a, b, m, switch_loglik(a, indata), switch_loglik(b, indata), switch_loglik(m, indata));
		m = *unku; 
		status = -1;
		goto exit_x;
	}
//       printf ("STATUS %s\n", status);

    
     
       //printf ("%5s [%9s, %9s] %9s %10s %9s\n",
        //       "iter", "lower", "upper", "min",
          //     "err", "err(est)");
     
       //printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
            //   iter, a, b,
             //  m, m - m_expected, b - a);
     
       do
         {
           iter++;
           status = gsl_min_fminimizer_iterate (s);
     
           m = gsl_min_fminimizer_x_minimum (s);
           a = gsl_min_fminimizer_x_lower (s);
           b = gsl_min_fminimizer_x_upper (s);
     
           status = gsl_min_test_interval (a, b, 0.001, 0.0);
     
//           if (status == GSL_SUCCESS)
       //      printf ("Converged:\n");
     
//           printf ("%5d [%.7f, %.7f] "
  //                 "%.7f %+.7f %.7f\n",
    //               iter, a, b,
      //             m, m - m_expected, b - a);
         }
       while (status == GSL_CONTINUE && iter < max_iter);
//	printf("Sucess \n");

	exit_x:   
 
	*unku = m;	

       gsl_min_fminimizer_free (s);   
 
       return status;


}

/* likelihood of switch given unku and param (switch effect at position 0) */
double switch_loglik(double unku, void *indata)
{

	int i;
	double mylik =0;
	
	double thispi = 0;

	double lasttime = 0;

	struct ad_pidata *pidata = (struct ad_pidata *) indata;

	double *param;
	gsl_vector *vec_event_time;
	gsl_vector_int *vec_event_type;
	int arm1_n;
	double *pi1;
	
	param = pidata->param;

	vec_event_time = pidata->vec_event_time;

	vec_event_type = pidata->vec_event_type;

	arm1_n = pidata->arm_n;

	pi1 = pidata->pi;

	//calculate proportion switchers given starting number unku
	propswitch(param, vec_event_time, vec_event_type, &arm1_n, &unku, pi1);

	int mynswitch = 0, mynnswitch=0;
	
	for(i=0; i<arm1_n-1; i++)
	{
		if(gsl_vector_get(vec_event_time, i) > lasttime) //update pi estimate and time estimate (otherwise multiple events at same time)
		{
			thispi = pi1[i];
			lasttime = gsl_vector_get(vec_event_time, i);
		}

		if( gsl_vector_int_get(vec_event_type, i) == 1) //stay in arm
		{
			mylik += log(GSL_MIN(1, GSL_MAX( 1-thispi, MYZERO ) ) ); //estimate of proportion of non-switchers just before switching
			if( (thispi>1) | (thispi<0)) mylik += OUTOFBOUND; //add penalty for out of bounds

			mynnswitch++;

		}
		else if( gsl_vector_int_get(vec_event_type, i) ==2) //switch
		{
			mylik += log(GSL_MIN(1,GSL_MAX(thispi,MYZERO))); //estimate of proportion switchers just before switching
	
			if( (thispi>1) | (thispi<0)) mylik += OUTOFBOUND; //add penalty for out of bounds
		
			mynswitch++;


		}

	}	


	//clean up
//	free(pi1);
//	gsl_vector_free(vec_event_time);
//	gsl_vector_int_free(vec_event_type);

	return(-mylik);

}

/* Proportion of switchers amoungst those not had opportunity to switch yet, given number at start (unku) */
void propswitch(double *param, gsl_vector *vec_event_time, gsl_vector_int *vec_event_type, int *arm1_n, double *unku, double *pi1)
{
	int i;
	double xi;
	double thisu;
	int thisn;
	double switcheff;
	double thispi;

	switcheff = exp(param[0]);

	thisu = *unku;

	thisn = *arm1_n;

	thispi = thisu / thisn;

	pi1[0] = thispi;
 
	for (i = 0; i < *arm1_n-1; i++)
	{
		if(gsl_vector_int_get(vec_event_type, i) == 0)	// event time before switching
		{
			xi = thispi * switcheff / ( (1-thispi) + thispi*switcheff);	
		}
		else if(gsl_vector_int_get(vec_event_type, i) == 1) //switch opportunity, stay in arm 1
		{
			xi = 0;
		}
		else if(gsl_vector_int_get(vec_event_type, i) == 2) //switch opportunity, cross over
		{
			xi = 1;	
		}
		else if(gsl_vector_int_get(vec_event_type, i) == 3) //censored before switch switch opportunity
		{
			xi = thispi;	

		}

		thisu = thisu - xi;

		thispi = gsl_max(0,(gsl_min(1,thisu / (*arm1_n - (i+1)))));


		pi1[i+1] = thispi;
	}


}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Full likelihood */

/* Switching likelihood. Mixture model from time zero. */
double Lik3s_mixture(double *param, double *outlik, double *datat, int *datad, double *datax, int *datag, int *datar, int *datan, int *datam, double *datas, int *datasidx, double *pr1, double *pr3, gsl_vector *vec_outbase)
{
	/*1. Initialise variables */
	int i, j;
	
	int arm1_n = 0;

	int arm2_n = 0;

	int lbnumswitch = 0, ubnumswitch = 0;

	int mytest = 0;	

	double mypl;

	double lik = 0;

	double sum = 0;

	double sum2 = 0;
	
	double unku = 0;

	double p1 = 0; double p2 = 0; double p3 = 0;

	gsl_matrix * mtx_X = gsl_matrix_alloc (*datan, *datam);

	gsl_vector * vec_baseline = gsl_vector_alloc (*datan);

	gsl_vector * vec_baseline_cuml = gsl_vector_alloc (*datan); //cumulative baseline hazard

	gsl_vector * vec_num = gsl_vector_alloc (*datan);

	gsl_vector * vec_num2 = gsl_vector_alloc (*datan);

	gsl_matrix * mtx_t = gsl_matrix_alloc (*datan, 4);
 
	gsl_matrix_set_zero(mtx_t);
	


	/*2. Populate vectors and matricies*/

	
	//Numer in each arm, needed to dynamically allocate next structures
	mytest = setupdata_narms(datat, datad, datag, datar, datan, &arm1_n, &arm2_n);

	//******
	//ARM 1
        gsl_vector * vec_event_time = gsl_vector_alloc (arm1_n);

        gsl_vector_int * vec_event_type = gsl_vector_int_alloc (arm1_n);

	//ARM 2
        gsl_vector * vec_event_arm2_time = gsl_vector_alloc (arm2_n);

        gsl_vector_int * vec_event_arm2_type = gsl_vector_int_alloc (arm2_n);

	//Setup data
	mytest = setupdata_tmev(datat, datad, datag, datar, datan, datas, vec_event_time, vec_event_type, vec_event_arm2_time, vec_event_arm2_type, &arm1_n, &ubnumswitch, &lbnumswitch);

	//use same definitions as partial likelihood
	inithazard(datan, datam, datax, param, mtx_X, mtx_t);


	//get the proportion of switcher through time
	//mypl = partial_lik_all(param, datat, datad, datag, datar, datan, datam, datas, datax, &unku, pr1, pr3, &arm1_n, &arm2_n, vec_event_time, vec_event_type, vec_event_arm2_time, vec_event_arm2_type, lbnumswitch, ubnumswitch, (int) 0);

    /*======================================================================================================*/
	double * allpi1 = malloc ( (*datan) * sizeof(double)); //not actually used in full lik, but part of object. not much extra mem reqmt so kept for moment
	double * allpi3 = malloc ( (*datan) * sizeof(double)); //not actuall used in full lik
//	double *myallpi[2];

	/* ============================================================
	 1. Set up data structures for proportion of switchers
	============================================================ */

	//double * pi1 = malloc ((arm1_n) * sizeof(double)); //proportion switcher at event times (not survival)
	//double * pi3 = malloc ((arm2_n) * sizeof(double));
	
//	double * pi1 = malloc ( (*datan) * sizeof(double));
//	double * pi3 = malloc ( (*datan) * sizeof(double));

	double pistart;
	//******

	//Setup data
        struct ad_pidata pidata= {param, vec_event_time, vec_event_type, arm1_n, allpi1, lbnumswitch, ubnumswitch};

        struct ad_pidata pidata2= {param, vec_event_arm2_time, vec_event_arm2_type, arm2_n, allpi3, lbnumswitch, ubnumswitch};

	/* ============================================================
	 3. Estimate proportion of switchers
	============================================================ */
	
	//3.1 find starting point
	if(fitswitch(&unku, &pidata)==-1) // no feasible starting u given current parameters, so lik out of bounds, no need for rest of loop
	{
		printf("HERE\n");
		lik += LIKPENALTY;
//		goto exit_x;
	}

	//3.2 get proportion switcher estimates
//	mytest = fitswitchprop(&unku, &pidata, &pidata2, datat, datar, datan, pr1, pr3, allpi1, allpi3); //not needed FL

	//3.3 point to proportion through time in 2-d pointer array. Not NEEDED FOR FULL LIK
	//myallpi[0] = allpi1;
	//myallpi[1] = allpi3;

	/* Full lik uses pi at start. no covariates here so same for all */
//	pistart = allpi1[0];
	pistart = unku / arm1_n;
	pr1[0] = pistart;
	
	printf("pi1, %g\n", pistart);

	/* Set number of events that take place at each time point 
	Would be more efficient to do this outside likelihood*/
	SetNum(datat, datad, vec_num, datan);

	SetNum2(datat, vec_num2, datan);

				
	/* 3. Calculate the baseline hazard function */
//	struct baseline_params_s bparam = {param, mtx_t, vec_baseline, vec_num, datat, datag, datar, datad, datas, datasidx, mtx_X, datan, allpi1, allpi3}; // pass all proportion points
	struct baseline_params_s bparam = {param, mtx_t, vec_baseline, vec_num, datat, datag, datar, datad, datas, datasidx, mtx_X, datan, &pistart, &pistart}; //only use proportion at baseline so just pass that


	int solv_code = 0;


/* Solve Basline*/
	solv_code = SolveBaseline_switch_mixture(&bparam);

	printf("Solver code = %d\n", solv_code);

	if(solv_code == 0) //could get baseline hazard
	{

		ConvertBaseline(vec_baseline, vec_num2);

		gsl_vector_memcpy(vec_outbase, vec_baseline);

		sum = 0;

 		printf("params %g, %g\n", param[0], param[1]);

		for(i=0; i < *datan; i++)
		{
			if(gsl_vector_get(vec_num2, i) > 0)
			{
				sum += gsl_vector_get(vec_baseline, i);
			}

			gsl_vector_set(vec_baseline_cuml, i, sum);
	

		}


		/* 5 Calculate the likelihood */
		for(i=0; i < *datan; i++)
	        {

                	if(datar[i]==1) //Control, arm with switching
	                {

				if(datag[i]==0)
				{
					sum = gsl_vector_get(vec_baseline_cuml, i); 

					//C NS
					p1 = CalcIfCensored( datad[i], gsl_vector_get(vec_baseline, i) * gsl_matrix_get(mtx_t, i, 0)) * exp( -gsl_matrix_get(mtx_t, i, 0) * sum);
					//C S
					p2 = CalcIfCensored( datad[i], gsl_vector_get(vec_baseline, i) * gsl_matrix_get(mtx_t, i, 1)) * exp( -gsl_matrix_get(mtx_t, i, 1) * sum);

					lik+= log( (1 - pistart) * p1 + pistart * p2);
				}
				else if(datag[i]==1) //post-switch, don't switch
				{
					sum = gsl_vector_get(vec_baseline_cuml, datasidx[i]); //pre-switch
 
					sum2 = gsl_vector_get(vec_baseline_cuml, i) - sum; //post-switch

					//C NS, pre-switch
					p1 = exp( -gsl_matrix_get(mtx_t, i, 0) * sum );

					//C S, pre-switch
					p2 = exp( -gsl_matrix_get(mtx_t, i, 1) * sum );
			
					//C NS, post-switch
					p3 = CalcIfCensored( datad[i], gsl_vector_get(vec_baseline, i) * gsl_matrix_get(mtx_t, i, 0)) * exp( -gsl_matrix_get(mtx_t, i, 0) * sum2);

					lik+= log( (1 - pistart) * p1 + pistart * p2) + log(p3);
				}
				else if(datag[i]==2) //post-switch, switch
				{

					sum = gsl_vector_get(vec_baseline_cuml, datasidx[i]); //pre-switch
 
					sum2 = gsl_vector_get(vec_baseline_cuml, i) - sum; //post-switch

					//C NS, pre-switch
					p1 = exp( -gsl_matrix_get(mtx_t, i, 0) * sum );

					//C S, pre-switch
					p2 = exp( -gsl_matrix_get(mtx_t, i, 1) * sum );
			
					//C S, post-switch
					p3 =  CalcIfCensored( datad[i], gsl_vector_get(vec_baseline, i) * gsl_matrix_get(mtx_t, i, 3)) * exp( -gsl_matrix_get(mtx_t, i, 3) * sum2);

					lik+= log( (1 - pistart) * p1 + pistart * p2) + log(p3);
				}
			}
	                else if(datar[i]==2) //unobserved switchers
			{
				sum = gsl_vector_get(vec_baseline_cuml, i); 
		
				//T NS
				p1 = CalcIfCensored( datad[i], gsl_vector_get(vec_baseline, i) * gsl_matrix_get(mtx_t, i, 2)) * exp( -gsl_matrix_get(mtx_t, i, 2) * sum);

				//T S
				p2 = CalcIfCensored( datad[i], gsl_vector_get(vec_baseline, i) * gsl_matrix_get(mtx_t, i, 3)) * exp( -gsl_matrix_get(mtx_t, i, 3) * sum);
				lik+= log( (1 - pistart) * p1 + pistart * p2);

			}
			else // TC -- no defined at present, datag[i]==3 should not exist!
			{
				printf("WARNING -- Only 2 arms allowed \n");
//				lik += log( CalcIfCensored(datad[i], gsl_vector_get(vec_baseline, i)) * gsl_matrix_get(mtx_p, i, 2));
				lik += log(0);
			}

	}

       		printf ("llikelihood = %g\n", lik);
	}
	else 
	{
		lik = -999999;
	}


	exit_x:
	
	//Clean up
	gsl_matrix_free (mtx_X);

	gsl_matrix_free (mtx_t);

	gsl_vector_free (vec_baseline);

	gsl_vector_free (vec_baseline_cuml);

	gsl_vector_free (vec_num);

	gsl_vector_free (vec_num2);

	//clean up
	gsl_vector_int_free (vec_event_type);

	gsl_vector_free (vec_event_time);

        gsl_vector_free(vec_event_arm2_time);

        gsl_vector_int_free(vec_event_arm2_type);

	free(allpi1);

	free(allpi3);

	//return -lik;
	*outlik = -lik;

	return -lik;


}
/* Switching likelihood, stratified. Mixture model from time zero. */
double Lik3s_mixture_strata(double *param, double *outlik, double *datat, int *datad, double *datax, int *datag, int *datar, int *datan, int *datam, double *datas, int *datasidx, double *pr1, double *pr3, gsl_vector *vec_outbase,  int *strata, int *nstrata)
{
	double thislik = 0;
	int i = 0;
        int * arr_arm1_n = malloc ( ((*nstrata)) * sizeof(int));
        int * arr_arm2_n = malloc ( ((*nstrata)) * sizeof(int));
        int * tempdatan = malloc  ( ((*nstrata)) * sizeof(int) );
        int tempdatan1  = 0;
        int tempdatan2 = 0;
	int mytest = 0;
        gsl_vector * vec_baseline;


        /* FOR EACH STRATA */
	for(i=0; i<*nstrata; i++)
	{

		tempdatan[i] = strata[i+1] - strata[i];

                //Numer in each arm, needed to for pr1 and pr3
   		vec_baseline = gsl_vector_alloc (tempdatan[i]);

		thislik += Lik3s_mixture(param, outlik, &datat[strata[i]], &datad[strata[i]], &datax[strata[i]], &datag[strata[i]], &datar[strata[i]], &tempdatan[i], datam, &datas[strata[i]], &datasidx[strata[i]], &pr1[tempdatan1], &pr3[tempdatan2], vec_baseline);

		gsl_vector_free(vec_baseline);

                //update starting point for pi records
           	mytest = setupdata_narms(&datat[strata[i]] , &datad[strata[i]], &datag[strata[i]] , &datar[strata[i]], &tempdatan[i], &arr_arm1_n[i], &arr_arm2_n[i]);
                tempdatan1 +=arr_arm1_n[i];
                tempdatan2 +=arr_arm2_n[i];
	}

	*outlik = thislik;


	//clean up
	free(arr_arm1_n);
	free(arr_arm2_n);
	free(tempdatan);

	return thislik;

}

