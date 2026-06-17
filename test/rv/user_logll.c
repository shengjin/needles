#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846
#define PI2 6.28318530717958647692

extern int ndim_data;

double newton_solver(double M, double e);

/*
 * 51 Peg b model used by the EMPEROR/reddemcee reference script:
 *   theta[0] = P      [days]
 *   theta[1] = K      [m/s]
 *   theta[2] = lambda [rad], lambda = phase + omega
 *   theta[3] = h      sqrt(e) sin(omega)
 *   theta[4] = k      sqrt(e) cos(omega)
 *   theta[5] = gamma_dot [m/s/day], linear acceleration
 *   theta[6] = gamma     [m/s], single-instrument offset
 *   theta[7] = jitter    [m/s]
 *
 * Input data are the Butler et al. 51 Peg Lick velocities with columns
 *   JD, RV, RV_error.
 * EMPEROR centres one-instrument RVs and shifts BJD so the first epoch is 0.
 */
double logll_beta(double *ptr_one_chain, int nline_data,
                  double *data_NlineNdim, double beta_one)
{
    if (ndim_data < 3 || nline_data <= 0) {
        return -1.0e300;
    }

    const double period_days = ptr_one_chain[0];
    const double K = ptr_one_chain[1];
    const double lambda = ptr_one_chain[2];
    const double h = ptr_one_chain[3];
    const double k = ptr_one_chain[4];
    const double gamma_dot = ptr_one_chain[5];
    const double gamma = ptr_one_chain[6];
    const double jitter = ptr_one_chain[7];

    const double ecc = h * h + k * k;
    if (period_days <= 0.0 || K <= 0.0 || ecc >= 1.0 || jitter <= 0.0) {
        return -1.0e300;
    }

    double omega = 0.0;
    if (ecc >= 1.0e-6) {
        omega = atan2(h, k);
        if (omega < 0.0) omega += PI2;
    }

    const double t_ref = data_NlineNdim[0];
    double rv_mean = 0.0;
    for (int i = 0; i < nline_data; i++) {
        rv_mean += data_NlineNdim[i * ndim_data + 1];
    }
    rv_mean /= (double)nline_data;

    double logll = -0.5 * log(PI2) * (double)nline_data;

    for (int i = 0; i < nline_data; i++) {
        const double t = data_NlineNdim[i * ndim_data + 0] - t_ref;
        const double rv_obs = data_NlineNdim[i * ndim_data + 1] - rv_mean;
        const double rv_err = data_NlineNdim[i * ndim_data + 2];
        const double err2 = rv_err * rv_err + jitter * jitter;

        if (err2 <= 0.0 || !isfinite(err2)) {
            return -1.0e300;
        }

        double M = PI2 * t / period_days + lambda - omega;
        M = fmod(M, PI2);
        if (M < 0.0) M += PI2;

        const double E = newton_solver(M, ecc);
        const double true_anomaly = 2.0 * atan2(sqrt(1.0 + ecc) * sin(E / 2.0),
                                                sqrt(1.0 - ecc) * cos(E / 2.0));
        const double trend_rv = gamma_dot * t;
        const double model_rv = K * (cos(true_anomaly + omega) + ecc * cos(omega))
                              + trend_rv + gamma;
        const double residual = rv_obs - model_rv;

        logll += -0.5 * (residual * residual / err2 + log(err2));
    }

    return beta_one * logll;
}


double newton_solver(double M, double e)
{
    double E = M;
    const double tolerance = 1.0e-12;

    for (int iter = 0; iter < 64; iter++) {
        const double f = E - e * sin(E) - M;
        const double fp = 1.0 - e * cos(E);
        const double step = f / fp;
        E -= step;
        if (fabs(step) < tolerance) break;
    }

    return E;
}
