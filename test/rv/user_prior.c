#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DEST_SIZE 100
#define PI 3.14159265358979323846
#define PI2 6.28318530717958647692

extern int N_parm;
extern double sigma_scale_min;
extern double sigma_scale_max;
extern double init_gp_ratio;
extern char *results_dir;

double *sigma_parm_min = NULL;
double *sigma_parm_max = NULL;
static double *para_min = NULL;
static double *para_max = NULL;

double log_prior(double *ptr_one_chain);
int calc_sigma_scale_boundary(double sigma_scale_min, double sigma_scale_max,
                              double *sigma_parm_min, double *sigma_parm_max);
double r8_unif_ab(double a, double b);
double r8_normal_ab(double a, double b);
char *read_onepara(char *path, char *para_name);
int save_debug_para_boundary(char *p_s, double p, double p_min, double p_max);

static double bounce_inside(double para, double p_min, double p_max)
{
    if (para > p_max) {
        para = p_max - (para - p_max);
        para = fmax(para, p_min);
    }
    if (para < p_min) {
        para = p_min + (p_min - para);
        para = fmin(para, p_max);
    }
    return para;
}

static double normal_cdf(double x, double mu, double sd)
{
    return 0.5 * (1.0 + erf((x - mu) / (sd * sqrt(2.0))));
}

static double log_uniform_prior(double x, double p_min, double p_max)
{
    if (x < p_min || x > p_max || p_max <= p_min) return -INFINITY;
    return -log(p_max - p_min);
}

static double log_trunc_normal_prior(double x, double p_min, double p_max,
                                     double mu, double sd)
{
    if (x < p_min || x > p_max || sd <= 0.0 || p_max <= p_min) {
        return -INFINITY;
    }

    const double z_norm = normal_cdf(p_max, mu, sd) - normal_cdf(p_min, mu, sd);
    if (z_norm <= 0.0 || !isfinite(z_norm)) return -INFINITY;

    const double z = (x - mu) / sd;
    return -0.5 * z * z - log(sd * sqrt(2.0 * PI)) - log(z_norm);
}

void read_parm_range(char *path)
{
    char para_name[DEST_SIZE];
    char dummy[DEST_SIZE] = {0};

    if (para_min == NULL) para_min = (double *)malloc(sizeof(double) * N_parm);
    if (para_max == NULL) para_max = (double *)malloc(sizeof(double) * N_parm);

    for (int i = 0; i < N_parm; i++) {
        snprintf(para_name, sizeof para_name, "para%d_max", i);
        char *para_line = read_onepara(path, para_name);
        if (para_line == NULL) exit(2);
        sscanf(para_line, "%[^:]:%lf", dummy, &para_max[i]);
        free(para_line);

        snprintf(para_name, sizeof para_name, "para%d_min", i);
        para_line = read_onepara(path, para_name);
        if (para_line == NULL) exit(2);
        sscanf(para_line, "%[^:]:%lf", dummy, &para_min[i]);
        free(para_line);
    }
}

void init_parm_set(int seed, double *chain_parm)
{
    srand48(seed);
    read_parm_range("input.ini");

    for (int i = 0; i < N_parm; i++) {
        chain_parm[i] = r8_unif_ab(para_min[i], para_max[i]);
    }

    if (N_parm == 8) {
        while (chain_parm[3] * chain_parm[3] + chain_parm[4] * chain_parm[4] > 1.0) {
            chain_parm[3] = r8_unif_ab(para_min[3], para_max[3]);
            chain_parm[4] = r8_unif_ab(para_min[4], para_max[4]);
        }
    }

    if (sigma_parm_min == NULL) sigma_parm_min = (double *)malloc(sizeof(double) * N_parm);
    if (sigma_parm_max == NULL) sigma_parm_max = (double *)malloc(sizeof(double) * N_parm);
    calc_sigma_scale_boundary(sigma_scale_min, sigma_scale_max,
                              sigma_parm_min, sigma_parm_max);
}

int calc_sigma_scale_boundary(double sigma_scale_min, double sigma_scale_max,
                              double *sigma_parm_min, double *sigma_parm_max)
{
    int sigma_bound_bad = 0;

    for (int i = 0; i < N_parm; i++) {
        sigma_parm_min[i] = sigma_scale_min * (para_max[i] - para_min[i]);
        sigma_parm_max[i] = sigma_scale_max * (para_max[i] - para_min[i]);

        if (sigma_parm_min[i] < 0.0 || sigma_parm_max[i] < 0.0 ||
            sigma_parm_max[i] < sigma_parm_min[i]) {
            sigma_bound_bad++;
        }
    }

    if (sigma_bound_bad > 0) {
        printf("ERR: bad sigma_parm_min/max %d TIMES!\n", sigma_bound_bad);
        return 1;
    }
    return 0;
}

int init_gaussian_proposal(double *ptr_sigma_prop, double init_gp_ratio)
{
    for (int i = 0; i < N_parm; i++) {
        ptr_sigma_prop[i] = (para_max[i] - para_min[i]) * init_gp_ratio;
    }
    return 0;
}

int para_boundary(double *ptr_one_chain_new)
{
    for (int i = 0; i < N_parm; i++) {
        if (N_parm == 8 && i == 2) {
            const double width = para_max[i] - para_min[i];
            while (ptr_one_chain_new[i] < para_min[i]) ptr_one_chain_new[i] += width;
            while (ptr_one_chain_new[i] > para_max[i]) ptr_one_chain_new[i] -= width;
        } else {
            ptr_one_chain_new[i] = bounce_inside(ptr_one_chain_new[i],
                                                 para_min[i], para_max[i]);
        }
    }
    return 0;
}

/*
 * Priors for the final 51 Peg model: P,K,lambda,h,k,
 * gamma_dot and gamma are uniform; e=h^2+k^2 is constrained
 * to [0,1]; jitter uses the EMPEROR truncated Normal(mu=5, sd=5).
 */
double log_prior(double *theta)
{
    if (N_parm != 8) return -INFINITY;

    const double ecc = theta[3] * theta[3] + theta[4] * theta[4];
    if (ecc < 0.0 || ecc > 1.0) return -INFINITY;

    double lp = 0.0;
    lp += log_uniform_prior(theta[0], para_min[0], para_max[0]);
    lp += log_uniform_prior(theta[1], para_min[1], para_max[1]);
    lp += log_uniform_prior(theta[2], para_min[2], para_max[2]);
    lp += log_uniform_prior(theta[3], para_min[3], para_max[3]);
    lp += log_uniform_prior(theta[4], para_min[4], para_max[4]);
    lp += log_uniform_prior(theta[5], para_min[5], para_max[5]);
    lp += log_uniform_prior(theta[6], para_min[6], para_max[6]);
    lp += log_trunc_normal_prior(theta[7], para_min[7], para_max[7], 5.0, 5.0);

    return lp;
}

int save_debug_para_boundary(char *p_s, double p, double p_min, double p_max)
{
    FILE *out;
    char fname[100];

    snprintf(fname, sizeof fname, "%s%s%s%s", results_dir, "/", p_s,
             ".debug_para_boundary");
    if ((out = fopen(fname, "a")) == NULL) {
        fprintf(stderr, "Can't create output file!\n");
        exit(3);
    }

    fprintf(out, "%f  %f  %f\n", p, p_min, p_max);

    if (fclose(out) != 0) fprintf(stderr, "Error in closing file!\n");
    return 0;
}
