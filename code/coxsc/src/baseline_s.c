#include "baseline_s.h"
#include "obs.h"

/* Baseline hazard where all the data, except the first hazard (value) are passed as a struct. Type of the struct is void because of the implementation of the root finding algorithms in the GSL. */
double CalculateBaseline_switch(const double value, void *params)
{

// param[0] = insistor effect
// param[1] = refuser effect
// param[2] = treatment effect
// param[3] = covariate effect

	struct baseline_params_s *bparam = (struct baseline_params_s *) params;
	
	double temp = value;

	double sum = 0.0;

	double censored_sum = 0.0;

	for(int i=0; i<(*bparam->datan); i++)
	{	

		if(temp <= 0.0)
		{
			return (temp - 1000.0*( (*bparam->datan) -i));
		}

		if(bparam->datad[i] == 0) //if censored, skip onto next
		{

			gsl_vector_set(bparam->vec_baseline, i, 0.0);

		}
		else
		{

			if(gsl_vector_get(bparam->vec_num, i) > 0) // vec_num indicates how many observations are events (not censored) at the same time point t, in the following obs[i] obs[i+1], ... with the same event time.
			{
				temp = temp - censored_sum;
				//printf("tempb %g\n", temp); //debug

				censored_sum = 0.0;

				if(temp <= 0.0)
				{
					return (temp - 1000.0*( (*bparam->datan)-i));
				}

				gsl_vector_set(bparam->vec_baseline, i, gsl_vector_get(bparam->vec_num, i) / temp);

				sum += gsl_vector_get(bparam->vec_baseline, i);

			}
			else
			{
				gsl_vector_set(bparam->vec_baseline, i , 0);
			}

			if(bparam->datag[i] == 1) //CT, switcher
			{
				if( bparam->datas[i] > bparam->datat[i]) //control, pre-switch
				{
					censored_sum += gsl_matrix_get(bparam->mtx_t, i, 0);
				}
				else //control, after switch
				{
					censored_sum += gsl_matrix_get(bparam->mtx_t, i, 2);
				}
				
			}

			else if(bparam->datag[i] == 0) //CC, non-switcher
			{				
				if( bparam->datas[i] > bparam->datat[i]) //control, pre-switch
				{
					censored_sum += gsl_matrix_get(bparam->mtx_t, i, 0);
				}
				else //control, after switch time
				{
					censored_sum += gsl_matrix_get(bparam->mtx_t, i, 1);
				}
			}
			else if(bparam->datag[i] == 2) // TT, treatment don't know if switcher or not
			{

				if( bparam->datas[i] > bparam->datat[i]) //treatment, pre-switch
				{

					censored_sum += gsl_matrix_get(bparam->mtx_t, i, 3);
				}
				else //treatment, after switch, mixture
				{

					double adC, adC1, adC2;

					//Don't switch part of C
					adC1 = (1 - bparam->pr1[i]) * CalcIfCensored(bparam->datad[i], gsl_matrix_get(bparam->mtx_t, i, 4)) * (exp(-CumHaz(0, i, bparam->vec_baseline, bparam->datasidx[i], gsl_matrix_get(bparam->mtx_t, i, 3), gsl_matrix_get(bparam->mtx_t, i, 4))));

					//Switch part of C
					adC2 = (bparam->pr1[i]) * CalcIfCensored(bparam->datad[i], gsl_matrix_get(bparam->mtx_t, i, 5)) * (exp( -CumHaz(0, i, bparam->vec_baseline, bparam->datasidx[i], gsl_matrix_get(bparam->mtx_t, i, 3), gsl_matrix_get(bparam->mtx_t, i, 5))));

					adC = 1 / (adC1 + adC2);

					censored_sum += ((adC1 * gsl_matrix_get(bparam->mtx_t, i, 4)) + (adC2 * gsl_matrix_get(bparam->mtx_t, i, 5))) * adC;

	//			printf("%g; %g; %g; %g; %g; %g \n", CalcCovariateEffect(bparam->param, bparam->mtx_X, &i), bparam->param[0], bparam->pr1[i], p1, bparam->param[2], p4);
				}
			}
			else //TC -- not implemented yet
			{
				printf("WARNING (baseline_s.c): no switching from treatment to control yet!\n");			
		
				censored_sum += 0 ;
			}
		
		}

	}

	temp = temp - censored_sum;

	if(temp < -1.0e7)
	{
		temp = -1.0e7;
	}


	return temp;
}

/* Base line hazard for mixture model at time zero */
/* NOTE UPDATED TO USE PL PARAMS 1/6/12 */
double CalculateBaseline_switch_mixture(const double value, void *params)
{

// param[0] = insistor effect
// param[1] = refuser effect
// param[2] = treatment effect
// param[3] = covariate effect

	struct baseline_params_s *bparam = (struct baseline_params_s *) params;
	
	double temp = value;

	double sum = 0.0;

	double censored_sum = 0.0;

	int lasti=0;

	for(int i=0; i<(*bparam->datan); i++)
	{	

		
		if(temp <= 0.0)
		{
			return (temp - 1000.0*( (*bparam->datan) -i));
		}

		if(bparam->datad[i] == 0) //if censored then will not cancel out in Delta_{i+k} - Delta_i for two non-censored events, so need to add to censored_sum (D_i in my notation), then skip onto next
		{
			gsl_vector_set(bparam->vec_baseline, i, 0.0);

		}
		else
		{

			if(gsl_vector_get(bparam->vec_num, i) > 0) // vec_num indicates how many observations are events (not censored) at the same time point t, in the following obs[i] obs[i+1], ... with the same event time.
			{
				/* new b-t ---------->*/

				for(int j = i; j< (*bparam->datan); j++)
				{
					if( bparam->datag[j]>0 ) // arm 0 observed switch 
					{

						if( (bparam->datasidx[j] > lasti) & (bparam->datasidx[j] <= i)) //then switched between last event and this one
						{
							if(bparam->datag[j] == 1) // don't switch
							{
								censored_sum += MixtureSumC(bparam, 0, 1, j,  bparam->datasidx[j], 0);
								censored_sum -= gsl_matrix_get(bparam->mtx_t, j, 0);
							}
							else //switch
							{
								censored_sum += MixtureSumC(bparam, 0, 1, j,  bparam->datasidx[j], 0);
								censored_sum -=  gsl_matrix_get(bparam->mtx_t, j, 3);
							}

						}
					}
				}


				lasti = i;

				/* end new b-t ---------->*/

				temp = temp - censored_sum;
				//printf("tempb %g\n", temp); //debug

				censored_sum = 0.0;

				if(temp <= 0.0)
				{
					if(isinf(-temp))
					{
						return (-2*1000.0*( (*bparam->datan)-i));
					}
					else
					{
						return (temp - 1000.0*( (*bparam->datan)-i));
					}
				}

				gsl_vector_set(bparam->vec_baseline, i, gsl_vector_get(bparam->vec_num, i) / temp);

				sum += gsl_vector_get(bparam->vec_baseline, i);

			}
			else
			{
				gsl_vector_set(bparam->vec_baseline, i , 0);
			}
			
		}
		if(bparam->datar[i]==1) // control
		{
			if(bparam->datag[i] == 0)
			{
				censored_sum += MixtureSumC(bparam, 0, 1, i,  i, 1);
			}
			else if (bparam->datag[i] == 1)	// don't switch given opp.
			{
				censored_sum += gsl_matrix_get(bparam->mtx_t, i, 0) ;
			}
			else if (bparam->datag[i] == 2)	// DO switch given opp.
			{
				censored_sum += gsl_matrix_get(bparam->mtx_t, i, 3);
			}
		}
		if (bparam->datar[i]==2) // trt
		{
			censored_sum += MixtureSumC(bparam, 2, 3, i,  i, 1);
		}
					

	}

	temp = temp - censored_sum;

	if(temp < -1.0e7)
	{
		temp = -1.0e7;
	}


	return temp;
}

/* Mixture calculation used in Calculating baseline 
preidx - pre-switching hazard

*/
double MixtureSumC(struct baseline_params_s *bparam, int xswidx, int swidx, int idx, int timeidx, int isterminal)
{
	double adC, adC1, adC2, myout;

	if(isterminal==1)
	{
		//Don't switch part of C
		adC1 = (1 - bparam->pr1[0]) * CalcIfCensored(bparam->datad[idx], gsl_matrix_get(bparam->mtx_t, idx, xswidx)) * exp(-CumHaz(0, timeidx, bparam->vec_baseline, bparam->datasidx[idx], gsl_matrix_get(bparam->mtx_t, idx, xswidx), gsl_matrix_get(bparam->mtx_t, idx, xswidx)));

		//Switch part of C
		adC2 =  bparam->pr1[0] * CalcIfCensored(bparam->datad[idx], gsl_matrix_get(bparam->mtx_t, idx, swidx)) * exp(-CumHaz(0, timeidx, bparam->vec_baseline, bparam->datasidx[idx], gsl_matrix_get(bparam->mtx_t, idx, swidx), gsl_matrix_get(bparam->mtx_t, idx, swidx)));


	}
	else
	{
		adC1 = (1 - bparam->pr1[0]) * exp(-CumHaz(0, timeidx, bparam->vec_baseline, bparam->datasidx[idx], gsl_matrix_get(bparam->mtx_t, idx, xswidx), gsl_matrix_get(bparam->mtx_t, idx, xswidx)));

		//Switch part of C
		adC2 = bparam->pr1[0] * exp(-CumHaz(0, timeidx, bparam->vec_baseline, bparam->datasidx[idx], gsl_matrix_get(bparam->mtx_t, idx, swidx), gsl_matrix_get(bparam->mtx_t, idx, swidx)));

	}

	
	if(adC1+adC2 == 0.0)
	{
		adC = 999999999999999;
	}
	else
	{
		adC = 1 / (adC1 + adC2);
	}


	myout = ((adC1 * gsl_matrix_get(bparam->mtx_t, idx, xswidx)) + (adC2 * gsl_matrix_get(bparam->mtx_t, idx, swidx))) * adC;

	return myout;
}

/* Brent-Dekker root finding algorithm to obtain the first hazard step.
The upper and lower bounds are the inverse of the hazard. i.e. if the function is evalued at 2.0, the hazard is 1/2.0 = 0.5 for the first time point */
int SolveBaseline_switch(struct baseline_params_s *bparam)
{

	gsl_set_error_handler_off(); //large likelihood penalty for such an outcome

       int status;
       int iter = 0, max_iter = 100;
       const gsl_root_fsolver_type *T;
       gsl_root_fsolver *s;
       double r = 0, r_expected = sqrt (5.0);
       double x_lo = 0.001, x_hi = 1000000.0; //arb - investigate choice x_hi
       
      	gsl_function F;

       F.function = &CalculateBaseline_switch;
       F.params = bparam;
     
       T = gsl_root_fsolver_brent;
       s = gsl_root_fsolver_alloc (T);
       gsl_root_fsolver_set (s, &F, x_lo, x_hi);
     
//       printf ("using %s method\n", 
   //            gsl_root_fsolver_name (s));
     
  //     printf ("%5s [%9s, %9s] %9s %10s %9s\n",
    //           "iter", "lower", "upper", "root", 
      //         "err", "err(est)");
     
       do
         {
           iter++;
           status = gsl_root_fsolver_iterate (s);
           r = gsl_root_fsolver_root (s);
           x_lo = gsl_root_fsolver_x_lower (s);
           x_hi = gsl_root_fsolver_x_upper (s);
           status = gsl_root_test_interval (x_lo, x_hi,
                                            0.0001,0); 
     
//           if (status == GSL_SUCCESS)
  //           printf ("Converged:\n");
     
        //   printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
          //         iter, x_lo, x_hi,
            //       r, r - r_expected, 
              //     x_hi - x_lo);
         }
       while (status == GSL_CONTINUE && iter < max_iter);
     
       gsl_root_fsolver_free (s);
     
       return status;

}

/* Brent-Dekker root finding algorithm to obtain the first hazard step.
The upper and lower bounds are the inverse of the hazard. i.e. if the function is evalued at 2.0, the hazard is 1/2.0 = 0.5 for the first time point */
int SolveBaseline_switch_mixture(struct baseline_params_s *bparam)
{

	gsl_set_error_handler_off(); //large likelihood penalty for such an outcome

       int status;
       int iter = 0, max_iter = 100;
       const gsl_root_fsolver_type *T;
       gsl_root_fsolver *s;
       double r = 0, r_expected = sqrt (5.0);
       double x_lo = 0.001, x_hi = 10000000.0; //arb - investigate choice x_hi
       
      	gsl_function F;

       F.function = &CalculateBaseline_switch_mixture;
       F.params = bparam;
     
       T = gsl_root_fsolver_brent;
       s = gsl_root_fsolver_alloc (T);
       gsl_root_fsolver_set (s, &F, x_lo, x_hi);
     
//       printf ("using %s method\n", 
   //            gsl_root_fsolver_name (s));
     
  //     printf ("%5s [%9s, %9s] %9s %10s %9s\n",
    //           "iter", "lower", "upper", "root", 
      //         "err", "err(est)");
     
       do
         {
           iter++;
           status = gsl_root_fsolver_iterate (s);
           r = gsl_root_fsolver_root (s);
           x_lo = gsl_root_fsolver_x_lower (s);
           x_hi = gsl_root_fsolver_x_upper (s);
           status = gsl_root_test_interval (x_lo, x_hi,
                                            0.0001,0);     
  //         if (status == GSL_SUCCESS)
    //         printf ("Converged:\n");
     
    //       printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
     //              iter, x_lo, x_hi,
      //             r, r - r_expected, 
       //            x_hi - x_lo);
         }
       while (status == GSL_CONTINUE && iter < max_iter);
     
       gsl_root_fsolver_free (s);
     
       return status;

}


/*
Fills in the hazards for time point ties in the data
*/
void ConvertBaseline(gsl_vector *vec_baseline, gsl_vector *vec_num)
{
        for(int i=0; i < vec_baseline->size; i++)
        {
                if(gsl_vector_get(vec_num, i) > 0)
                {
                        for(int j=0; j<gsl_vector_get(vec_num, i); j++)
                        {
                                gsl_vector_set(vec_baseline, i+j, gsl_vector_get(vec_baseline,i));
                        }
                }
        }
}


/*
Function to return the cumulative baseline hazard between two time point time1 and time2 (indicies)
*/
double CumBaseHaz(int time1, int time2, gsl_vector *vec_baseline)
{
	double mysum;

	mysum = 0;

	for(int i=time1; i<time2; i++)
	{
		mysum += gsl_vector_get(vec_baseline, i);
		

	}


	return mysum;
}

/*
Function to return the cumulative hazard between two time point time1 and time2 (indicies)
*/
double CumHaz(int time1, int time2, gsl_vector *vec_baseline, int switchidx, double psi1, double psi2)
{

	double mychaz1, mychaz2;
	
	mychaz1 = 0;

	mychaz2 = 0;

	if( time1< switchidx )
	{
		mychaz1 = CumBaseHaz(time1, minints(switchidx, time2), vec_baseline) * psi1;	
	}

	if (time2> switchidx)
	{
		mychaz2 = CumBaseHaz(switchidx, time2, vec_baseline) * psi2;

	}
	return mychaz1 + mychaz2;
}

int maxints(int v1, int v2)
{
        if(v1>v2)
        {
                return v1;
        }
        else
        {
                return v2;
        }
}

int minints(int v1, int v2)
{
        if(v1<v2)
        {
                return v1;
        }
        else
        {
                return v2;
        }
}

