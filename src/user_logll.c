
#include "nbody.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "alloc.h"

extern int ndim_data; 




extern double stellar_mass; //the stellar mass in unit of Solar mass

// check ttvFast stellar mass

extern double t_start_yr;

extern double t_total_yr;



// do Nbody integration and calc TTV
int nbody_and_ttv(double *ptr_one_chain, double stellar_mass, double t_start_yr, double t_total_yr, int *ptr_n_transit_p0, int *ptr_n_transit_p1, double *t_transit_p0, double *t_transit_p1);



double logll_beta(double *ptr_one_chain, int nline_data, double *data_NlineNdim, double beta_one)
{


    int n_transit_p0 = 0;
    int n_transit_p1 = 0;
    //
    int *ptr_n_transit_p0;
    int *ptr_n_transit_p1;
    ptr_n_transit_p0 = &n_transit_p0;
    ptr_n_transit_p1 = &n_transit_p1;


    double *t_transit_p0; 
    double *t_transit_p1; 
    t_transit_p0 = alloc_1d_double(N_TRANSIT_MAX);
    t_transit_p1 = alloc_1d_double(N_TRANSIT_MAX);

    nbody_and_ttv(ptr_one_chain, stellar_mass, t_start_yr, t_total_yr, ptr_n_transit_p0, ptr_n_transit_p1, t_transit_p0, t_transit_p1);

/*
    printf("%lf  \n", beta_one);
if (beta_one > 0.5)
{
    printf("%d %d \n", *ptr_n_transit_p0, *ptr_n_transit_p1);



    for (int k = 0; k < *ptr_n_transit_p0; k++)
    {
        printf("%d %d %lf \n", *ptr_n_transit_p0, k, t_transit_p0[k]);
    }

    for (int k = 0; k < *ptr_n_transit_p1; k++)
    {
        printf("%d %d %lf \n", *ptr_n_transit_p1, k, t_transit_p1[k]);
    }
}
*/


    // calc log likelihood
    // 
    double logll;
    logll = 0;

    
    // real TTV: data_NlineNdim[i*ndim_data+j]
    // model: t_transit_p0, t_transit_p1
    
    // if the first is starting from 0.0, then subtract t_0 is not needed
    // t_0=data_NlineNdim[2];
    
    for (int i = 0; i < nline_data; i++)
    {
        double sigma_i, t_transit_i;
	int planet_id, epoch_id;
        sigma_i = data_NlineNdim[i*ndim_data+3]/1440.0;
        t_transit_i = data_NlineNdim[i*ndim_data+2];
	planet_id = data_NlineNdim[i*ndim_data+0];
	epoch_id = data_NlineNdim[i*ndim_data+1];

	double residual;

        if (planet_id == 0)
	{
	    if (epoch_id <= *ptr_n_transit_p0)
	    {
	        residual = t_transit_i - t_transit_p0[epoch_id];
                logll = logll -0.5*log(TPI) - log(sigma_i) - 0.5*pow((residual/sigma_i), 2.0);
                 //printf("%lf,%lf\n",t_transit_i,epoch_TT_0[epoch_id]);
	    }
	    else
	    {
	        // TODO: find a better way for residual in this case
	        residual = (t_transit_i - t_transit_p0[*ptr_n_transit_p0-1])*2.0;
                //logll = logll - sigma_i - 0.5*pow((residual/sigma_i), 2.0);
                logll = logll -0.5*log(TPI) - log(sigma_i) - 0.5*pow((residual/sigma_i), 2.0);
	    }
	}
        //
        if (planet_id == 1)
	{
	    if (epoch_id <= *ptr_n_transit_p1)
	    {
	        residual = t_transit_i - t_transit_p1[epoch_id];
                //logll = logll - sigma_i - 0.5*pow((residual/sigma_i), 2.0);
                logll = logll -0.5*log(TPI) - log(sigma_i) - 0.5*pow((residual/sigma_i), 2.0);
	    }
	    else
	    {
	        // TODO: find a better way for residual in this case
	        residual = (t_transit_i - t_transit_p1[*ptr_n_transit_p1-1])*2.0;
                logll = logll -0.5*log(TPI) - log(sigma_i) - 0.5*pow((residual/sigma_i), 2.0);
                //logll = logll - sigma_i - 0.5*pow((residual/sigma_i), 2.0);
	    }
	}

	// for debug
	//if (beta_one == 1.0)
	
	//{
           // printf("%lf\n", logll);
	//}
	//
    }
    

    /////////////////////////////////////////////////////
    free_1d_double(t_transit_p0);
    t_transit_p0 = NULL;
    free_1d_double(t_transit_p1);
    t_transit_p1 = NULL;


    ////////////////////////////
    // final tempering
    //logll = pow(logll, beta_one);
    // or
    logll = logll*beta_one;

    //printf("%lf\n", logll);
    return logll;
    //
}







