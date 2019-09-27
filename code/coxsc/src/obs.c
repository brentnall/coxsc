#include "obs.h"

double CalcIfCensored(int censind, double x)
{
        return (censind == 1 ? x : 1);
}


/* bX for individual i */
double CalcCovariateEffect(const double *param, const gsl_matrix *mtx_X, const int *inid)
{
        double result = 0.0;
        
        for(int j=0; j< mtx_X->size2; j++)
        {
                result += gsl_matrix_get(mtx_X, *inid, j)  * param[j+3];
        }

        return result;


}

/* bX for individual i */
double CalcCovariateEffect_mixture(const double *param, const gsl_matrix *mtx_X, const int *inid)
{
        double result = 0.0;
        
        for(int j=0; j< mtx_X->size2; j++)
        {
                result += gsl_matrix_get(mtx_X, *inid, j)  * param[j+2];
        }

        return result;


}
/* 
This sets class variable num, to indicate how many observations are event sat the same value
The vector num() is the same length as the observations. Information on the number of events at the same value is put at the first observation with the same time value, the following time points are zero.
*/
void SetNum(double *datat, int *datad, gsl_vector *vec_num, int *datan)
{
        int i;

        gsl_vector_set_zero(vec_num);


        for(i=0; i < *datan;)
        {
                int temp = 0;

                int same = 0;
        
                int orig = i;
                
                do
                {
                        temp += datad[i];

                        same ++;
                        
                        i++;
                }
                while(i< *datan && datat[i] == datat[i-1]);

                gsl_vector_set(vec_num, orig, temp);
        }
}
/* 
This sets class variable num2, to indicate how many observations are tied at the same value
*/
void SetNum2(double *datat, gsl_vector *vec_num2, int *datan)
{
        int i;

	gsl_vector_set_zero(vec_num2);

        for(i=0; i < *datan;)
        {

                int same = 0;
        
                int orig = i;
                
                do
                {
                        same ++;
                        
                        i++;
                }
                while(i< *datan && datat[i] == datat[i-1]);

                gsl_vector_set(vec_num2, orig, same);
        }
}





