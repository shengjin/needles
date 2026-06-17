#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nbody.h"

#include <time.h>


extern double err_tolerance_intg; // = 1e-14;

const double h_init_scale = 0.01;  // init step size in days


// linear interpolate
double LinIntP(double x, double x1, double y1, double x2, double y2);

// RKF78 in case derivative function does not explicitly depend on time t.
int rkf78_no_t(double* ptr_t, double* ptr_dt, double err_tol, double* x, int n);
// RKF78 in case derivative function depend on time t, e.g., Earth's satelites.
//int rkf78(double* ptr_t, double* ptr_dt, double err_tol, double* x, int n);

//  angular conversion 
double zatan(double sina, double cosa);

// orbit elements to Cartesian coordinates
int ele2rat(double a, double e, double ci, double w, double omega, double cm0, int ii, double xout[6]);

// Cartesian to elements
int rat2ele(double xin[6], double* a, double* e, double* ci, double* w, double* omega, double* cm, int ii, int ito);

// do Nbody integration and calc TTV
int nbody_and_ttv(double *ptr_one_chain, double stellar_mass, double t_start_yr, double t_total_yr, int *ptr_n_transit_p0, int *ptr_n_transit_p1, double *t_transit_p0, double *t_transit_p1);

double xmass[MAX_BODIES] = { 0 };

// array to store Masses in unit of Earth Mass


int nbody_and_ttv(double *ptr_one_chain, double stellar_mass, double t_start_yr, double t_total_yr, int *ptr_n_transit_p0, int *ptr_n_transit_p1, double *t_transit_p0, double *t_transit_p1)
{

    
    // max_bodies = 2 planets + star
    double p0[MAX_BODIES], e0[MAX_BODIES], ci0[MAX_BODIES], w0[MAX_BODIES];
    double omg0[MAX_BODIES], cm00[MAX_BODIES];

    // i=0 refer to the center star
    xmass[0] = stellar_mass;
    //
    xmass[1] = (*ptr_one_chain)*MEARTH_2_MSUN; // M1
    p0[1] = *(ptr_one_chain+1); // p1
    e0[1] = *(ptr_one_chain+2); // e1
    ci0[1] = *(ptr_one_chain+3); // i1
    omg0[1] = 0;
    w0[1] = *(ptr_one_chain+4); // omg1
    cm00[1] = *(ptr_one_chain+5);  // M01
    // 
    xmass[2] = *(ptr_one_chain+6)*MEARTH_2_MSUN; // M2
    p0[2] = *(ptr_one_chain+7); // p2
    e0[2] = *(ptr_one_chain+8); // e2
    ci0[2] = *(ptr_one_chain+9); // i2
    omg0[2]=  *(ptr_one_chain+10); // Omg2-Omg1
    w0[2] = *(ptr_one_chain+11); // omg2
    cm00[2] = *(ptr_one_chain+12); // M02
    
    //jacobi coordinate convertion mimicing TTVFast
    //the input are orbital elements coordinates in jacobi coords we just need to modify the central mass for the planet 2
    //planet 2 orbits CoM(star + pl.1)   (central mass = M_star + M_1)
    //planet i orbits CoM(star + pl.1 + ... + pl.i-1)
    //we need to recompute semi-major axis using the Jacobi central mass then transform Jacobi Cartesian -> Barycentric Cartesian
    //jacobi central mass: mu_jac[i] = sum of masses of bodies 0..i-1
    //
    //
    double mu_jac[MAX_BODIES];
    mu_jac[0] = 0;                          // not used anywhere, initialize to 0.
    //
    // the index is speciallized as it is, to make functions (ele2rat, and the following part...) use same index i for each body.
    mu_jac[1] = xmass[0];                          // planet 1: mass of the star
    for (int i = 2; i < NBOD; i++)
    {
        mu_jac[i] = mu_jac[i-1] + xmass[i-1];     // add previous planet
    }




    // 1: calc semi_major axes; 
    // 2: convert degree to radian 
    double a0[MAX_BODIES], a[MAX_BODIES];
    double e[MAX_BODIES], ci[MAX_BODIES], w[MAX_BODIES], omg[MAX_BODIES], cm0[MAX_BODIES], cm[MAX_BODIES];


    // NOTE should use center of Mass rather than star only
    for (int i = 1; i < NBOD; i++)
    {
        // Kepler's 3rd law: a^3 = P^2 * (mu_jac + m_planet)
        a0[i] = pow(pow(p0[i] / YR_DAYS, 2) * (mu_jac[i] + xmass[i]), 1.0 / 3.0);
        a[i] = a0[i];
        e[i] = e0[i];
        ci[i] = ci0[i] * TPI / DEGREE_2_PI;  // 转换为弧度
        w[i] = w0[i] * TPI / DEGREE_2_PI;
        omg[i] = omg0[i] * TPI / DEGREE_2_PI;
        cm0[i] = cm00[i] * TPI / DEGREE_2_PI;
    }


    // initialzie the state vector for positions and velocities
    double x[MAX_BODIES][6] = { {0} };  // x,y,z and velocity x, v_y, v_z
    //
    double xp[6];                       // temporary store the positions of one body
    double xrk[MAX_BODIES * 6];         // 1-dimension array used for rkf78 integrator


    //converting Jacobi elements to Jacobi Cartesian 
    //ele2rat uses (xmass[0] + xmass[ii]) internally as the gravitational parameter to keep it consistent with the code and use (mu_jac[i] + xmass[i]) we temporarily set xmass[0] = mu_jac[i] before calling ele2rat.
    double xjac[MAX_BODIES][6] = {0};

    double save_m0 = xmass[0];
    for (int i = 1; i < NBOD; i++) 
    {
    //
    // because ele2rat use xmass[0] and xmass[ii] (extern variable)
        xmass[0] = mu_jac[i];   // make ele2rat use Jacobi central mass.
                                // mu_jac[i] = m_star + ... + m_p_(i-1)
        //
        ele2rat(a[i], e[i], ci[i], w[i], omg[i], cm0[i], i, xp);
	//
        for (int j = 0; j < NDIM; j++) 
	{
	    xjac[i][j] = xp[j];
	}
    }
    xmass[0] = save_m0;         // restore real star mass


    //from jacobi Cartesian to Barycentric Cartesian
    //jacobi position r_jac[i] = position of body i relative to the center-of-mass of bodies 0, 1, ..., i-1.


    //First convert to astrocentric (relative to star, heliocentric notation used in TTVFast), then to barycentric.
    //  Solar System Dynamics P441
    //
    //Astrocentric from jacobi
    //r_astro[1] = r_jac[1]
    //r_astro[i] = r_jac[i] + CoM_astro(bodies 0..i-1) ; CoM_astro = center of mass of bodies 0..i-1 measured from star
    //eta[i] = m_0 + m_1 + ... + m_i tracking the total com

    double eta[MAX_BODIES];
    double com_astro[6];   // CoM of bodies 0..i-1
    for (int j = 0; j < NDIM; j++)
    {
        com_astro[j] = 0.0;  // star at origin, astrocentric, heliocentric
    }

    eta[0] = xmass[0];

    for (int i = 1; i < NBOD; i++)
    {
        eta[i] = eta[i-1] + xmass[i];
        //
        // Astrocentric position of planet i = Jacobi position + CoM of inner bodies
        for (int j = 0; j < NDIM; j++)
	{
            x[i][j] = xjac[i][j] + com_astro[j];
        }
        //
        // Update running CoM
        for (int j = 0; j < NDIM; j++)
	{
            com_astro[j] = (eta[i-1] * com_astro[j] + xmass[i] * x[i][j]) / eta[i];
        }
    }


    // Second: x[1..NBOD-1] are astrocentric now
    //  then to barycentric:
    /*
    double tmass = 0.0;
    for (int i = 0; i < NBOD; i++) tmass += xmass[i];
    for (int j = 0; j < NDIM; j++) {
        x[0][j] = 0.0;
        for (int i = 1; i < NBOD; i++) {
            x[0][j] -= xmass[i] * x[i][j];
        }
        x[0][j] /= tmass;
    }
    */
    // comment above section because acutually: com_astro = -x[0][j] 
    for (int j = 0; j < NDIM; j++) 
    {
        x[0][j] = -com_astro[j];
    }
    //
    //
    //
    // re-calculate the positions relative to the center of mass
    for (int i = 1; i < NBOD; i++) 
    {
        for (int j = 0; j < NDIM; j++) {
            x[i][j] += x[0][j];
        }
    }



    // convert to 1-dimension array for the rkf78 integrator
    for (int i = 0; i < NBOD; i++)
    {
        for (int j = 0; j < NDIM; j++)
        {
            xrk[i * NDIM + j] = x[i][j];
        }
    }

    // control parameters in integration
    double h = TPI / YR_DAYS * h_init_scale;    // initial step siz
    double err = err_tolerance_intg; // err tolerance
    double t = t_start_yr * TPI;                     // initial integration time
    double tstop = t_total_yr * TPI;       // total integration time


    /* ---------- TTV state ---------- */

    // assuming we are observing from the +z axis (the sky is in the x-y plane)
    const int OBS_Z_SIGN = +1;


    // the main loop
    while (t < tstop) 
    {  //start of main loop
        // updat the coordinates
        for (int i = 0; i < NBOD; i++)
        {
            for (int j = 0; j < NDIM; j++)
            {
                x[i][j] = xrk[i * NDIM + j];
            }
        }

        // convert the Cartesian coordinates to orbit elements
	/* speed up, skip this and stability check at the end of each iteration
        for (int i = 1; i < NBOD; i++)
        {
            for (int j = 0; j < NDIM; j++)
            {
                xp[j] = x[i][j] - x[0][j];
            }
            rat2ele(xp, &a[i], &e[i], &ci[i], &w[i], &omg[i], &cm[i], i, 0);
        }
	*/


        double t_old = t;              

	///////////////////////////////////
        // Do integration once
        rkf78_no_t(&t, &h, err, xrk, NBOD * NDIM);

        /* start TTV detection (every adaptive step) ---------- */
        double t_new = t;
        if (t_new > t_old) // somettime  rkf78 donot integrate due to large err, i.e., it will re-integrate using smaller step 
	{ /* start TTV detection (every adaptive step) ---------- */
            for (int i = 1; i < NBOD; i++) 
	    { // start for each body

                // old vector_r_perp (r_pl_star) before this integration 
                double rx0 = x[i][0] - x[0][0];
                double ry0 = x[i][1] - x[0][1];
                double rz0 = x[i][2] - x[0][2];
                // old vector_v_perp at x_y plane (r_pl_star) before this integration 
                double vx0 = x[i][3] - x[0][3];
                double vy0 = x[i][4] - x[0][4];

                // new vector_r_perp (r_pl_star) after this integration 
                double rx1 = xrk[i*NDIM + 0] - xrk[0*NDIM + 0];
                double ry1 = xrk[i*NDIM + 1] - xrk[0*NDIM + 1];
                double rz1 = xrk[i*NDIM + 2] - xrk[0*NDIM + 2];
                // new vector_v_perp at x_y plane (r_pl_star) before this integration 
                double vx1 = xrk[i*NDIM + 3] - xrk[0*NDIM + 3];
                double vy1 = xrk[i*NDIM + 4] - xrk[0*NDIM + 4];

                // source function = vector_r_perp · vector_v_perp 
		//  Math explanation:  We would like to find the: min(vector_r_perp)
		//       We construct a target function which is:  1/2 * vector_r_perp**2.0
		//       Then, the derivative of the target function is the source function
		//       At the transit middle time, the source function will be passing the zero point, i.e., changing from negative (s0, before) to positive (s1, after)
                double s0 = rx0*vx0 + ry0*vy0;
                double s1 = rx1*vx1 + ry1*vy1;

                // Transit middle time: source function is passing 0 point, i.e., changing from negative value to positive
                if (s0 < 0.0 && s1 >= 0.0)
		{
	
                    double t_tr = LinIntP(0.0, s0, t_old, s1, t_new);
                    double t_tr_days = (t_tr / TPI) * YR_DAYS;

                    // rxm
                    //double rxm = LinIntP(0.0, s0, rx0, s1, rx1);
                    //double rym = LinIntP(0.0, s0, ry0, s1, ry1);
                    double rzm = LinIntP(0.0, s0, rz0, s1, rz1);

                    //double dperp = sqrt(rxm*rxm + rym*rym);

                    // skip the case that host star is between planet and the observers
                    if (OBS_Z_SIGN * rzm <= 0.0) continue;


                    if (i == 1)
		    {
                        int idx = (*ptr_n_transit_p0);
		        t_transit_p0[idx] = t_tr_days;
                        //printf("%d %.10f \n", idx, t_tr_days);
                        (*ptr_n_transit_p0)++;
		    }

                    if (i == 2)
		    {
                        int idx = (*ptr_n_transit_p1);
                        //printf("%d %.10f \n", idx, t_tr_days);
		        t_transit_p1[idx] = t_tr_days;
                        (*ptr_n_transit_p1)++;
		    }

                }

           }  // end for each body

        } /* end TTV detection (every adaptive step) ---------- */


        /* 稳定性检查 */
	/*
        for (int i = 1; i < NBOD; i++) {
            if (fabs(a[i] - a0[i]) > 0.1 * a0[i]) { break; }
        }
	*/

    }   //end of main loop


    return 0;
}




double LinIntP(double x, double x1, double y1, double x2, double y2)
{
    // makes linear interpolation for f(x) at the point x, if (x1,f(x1)=y1) and (x2,f(x2)=y2)
    // using Lagrange's formula
    if (x == x1)
    {
        return y1;
    }
    else if (x == x2)
    {
        return y2;
    }
    else if (x1 == x2)
    {
        return y2;
    }
    else
    {
        //printf("%lf %lf %lf %lf %lf %lf \n", x1,y1,x2,y2,x, ((x-x2)/(x1-x2))*y1+((x-x1)/(x2-x1))*y2);
        return ((x-x2)/(x1-x2))*y1+((x-x1)/(x2-x1))*y2;
    }
}



