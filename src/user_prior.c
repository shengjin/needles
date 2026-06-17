
//////////////////////////////////
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DEST_SIZE 100

// NOTE: the description of all parameters starts from line 714

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////// Part 0
////// General part, no need to change in most cases
//////////////////////////////////////////////////////////////////////////////
//
// Number of parameters
extern int N_parm;
//
// Gaussian proposal array
double *sigma_parm_min;
double *sigma_parm_max;
//
// Range to scale gaussian proposal
extern double sigma_scale_min;
extern double sigma_scale_max;
//
// initial gaussian proposal ratio to the paramter range
extern double init_gp_ratio;
//
// Define results directory
extern char *results_dir; //dir to save outputs
//
/////////////////////////////////////
//NOTE: change log_prior in mpi_batch/init if func_prototype is changed
double log_prior(double *ptr_one_chain);
//
////////////////////
//calc sigma scale boundary (function in user_prior.c)
int calc_sigma_scale_boundary(double sigma_scale_min, double sigma_scale_max, double *sigma_parm_min, double *sigma_parm_max);
//
//
// useful random number genetors
double r8_normal_ab ( double a, double b );
double r8_unif_ab ( double a, double b );
double r8_logunif_ab ( double a, double b );
double normal_pdf(double x, double mean, double sd);
//
//
///////////////////
// prototype no need to chaneg, contents need change
void read_parm_range(char *path);
//
//
///////////////////////////////////////////////
// check if modification is needed or not
///////////////////
int save_debug_para_boundary(char *p_s, double p, double p_min, double p_max);
//
//

double bounce_inside(double para, double para_min, double para_max)
{
//
    if (para > para_max)
    {
        para = para_max - (para - para_max);
        para = fmax(para, para_min);
    }
    if (para < para_min)
    {
        para = para_min + (para_min - para);
        para = fmin(para, para_max);
    }
    //
    return para;
}


//////////////////////////////////////////////////////////////////////////////
////// Part 1
////// define how many parameters in the model
////// their values are read from input files by read_parm_range(char *path)
//
double para0_min;
double para1_min;
double para2_min;
double para3_min;
double para4_min;
double para5_min;
double para6_min;
double para7_min;
double para8_min;
double para9_min;
double para10_min;
double para11_min;
double para12_min;
//
double para0_max;
double para1_max;
double para2_max;
double para3_max;
double para4_max;
double para5_max;
double para6_max;
double para7_max;
double para8_max;
double para9_max;
double para10_max;
double para11_max;
double para12_max;
//
/////// Part 2: read prior range
///////   NOTE: Input Min and Max in default.
///////   NOTE: Add your own control para (mean, std, etc) in case of needed.
////////////////////////////////////////////
//
void read_parm_range(char *path)
{
    char *para_line;
    char *para_name;
    //
    char *read_onepara(char *path, char *para_name);
    //
    char dummy[DEST_SIZE] = {0};
    //
    //
    para_name = "para0_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para0_max);
    para_name = "para0_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para0_min);
    //
    para_name = "para1_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para1_max);
    para_name = "para1_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para1_min);
    //
    para_name = "para2_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para2_max);
    para_name = "para2_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para2_min);
    //
    para_name = "para3_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para3_max);
    para_name = "para3_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para3_min);
    //
    para_name = "para4_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para4_max);
    para_name = "para4_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para4_min);
    //
    para_name = "para5_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para5_max);
    para_name = "para5_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para5_min);
    //
    para_name = "para6_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para6_max);
    para_name = "para6_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para6_min);
    //
    para_name = "para7_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para7_max);
    para_name = "para7_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para7_min);
    //
    para_name = "para8_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para8_max);
    para_name = "para8_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para8_min);
    //
    para_name = "para9_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para9_max);
    para_name = "para9_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para9_min);
    //
    para_name = "para10_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para10_max);
    para_name = "para10_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para10_min);
    //
    para_name = "para11_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para11_max);
    para_name = "para11_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para11_min);
    //
    para_name = "para12_max";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para12_max);
    para_name = "para12_min";
    para_line = read_onepara(path, para_name);
    sscanf(para_line, "%[^:]:%lf", dummy, &para12_min);
    //

    //
    free(para_line);
    para_line = NULL;
}

///////////////////////////////////////////////////////////////////////////////
///////   Part 3: initialize prior, part 1 
///////   NOTE: uniform initialization in default.
///////   NOTE: Specify your own init function in case of needed.
////////////////////////////////////////////
//
// set init parameter at N_ITER = 0
double para0_init(double para0_min,double para0_max)
{
    double para0;
    para0 = r8_unif_ab(para0_min, para0_max);
    return para0;
}
//
double para1_init(double para1_min,double para1_max)
{
    double para1;
    para1 = r8_unif_ab(para1_min, para1_max);
    return para1;
}
//
double para2_init(double para2_min,double para2_max)
{
    double para2;
    para2 = r8_unif_ab(para2_min, para2_max);
    return para2;
}
//
double para3_init(double para3_min,double para3_max)
{
    double para3;
    para3 = r8_unif_ab(para3_min, para3_max);
    return para3;
}
//
double para4_init(double para4_min,double para4_max)
{
    double para4;
    para4 = r8_unif_ab(para4_min, para4_max);
    return para4;
}
//
double para5_init(double para5_min,double para5_max)
{
    double para5;
    para5 = r8_unif_ab(para5_min, para5_max);
    return para5;
}
//
double para6_init(double para6_min,double para6_max)
{
    double para6;
    para6 = r8_unif_ab(para6_min, para6_max);
    return para6;
}
//
double para7_init(double para7_min,double para7_max)
{
    double para7;
    para7 = r8_unif_ab(para7_min, para7_max);
    return para7;
}
//
double para8_init(double para8_min,double para8_max)
{
    double para8;
    para8 = r8_unif_ab(para8_min, para8_max);
    return para8;
}
//
double para9_init(double para9_min,double para9_max)
{
    double para9;
    para9 = r8_unif_ab(para9_min, para9_max);
    return para9;
}
//
double para10_init(double para10_min,double para10_max)
{
    double para10;
    para10 = r8_unif_ab(para10_min, para10_max);
    return para10;
}
//
double para11_init(double para11_min,double para11_max)
{
    double para11;
    para11 = r8_unif_ab(para11_min, para11_max);
    return para11;
}
//
double para12_init(double para12_min,double para12_max)
{
    double para12;
    para12 = r8_unif_ab(para12_min, para12_max);
    return para12;
}
//

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
///////   Part 4 initialize, part 2
///////   NOTE: uniform initialization in default.
///////   NOTE: Specify your own init function in case of needed.
//
// set init parameter at N_ITER = 0
void init_parm_set(int seed, double* chain_parm)
{
    // rand seed
    srand48(seed);
    //
    char *input_file; 
    input_file = "input.ini";
    // read the range
    read_parm_range(input_file);
    //
    //
    chain_parm[0] = para0_init(para0_min, para0_max);
    chain_parm[1] = para1_init(para1_min, para1_max);
    chain_parm[2] = para2_init(para2_min, para2_max);
    chain_parm[3] = para3_init(para3_min, para3_max);
    chain_parm[4] = para4_init(para4_min, para4_max);
    chain_parm[5] = para5_init(para5_min, para5_max);
    chain_parm[6] = para6_init(para6_min, para6_max);
    chain_parm[7] = para7_init(para7_min, para7_max);
    chain_parm[8] = para8_init(para8_min, para8_max);
    chain_parm[9] = para9_init(para9_min, para9_max);
    chain_parm[10] = para10_init(para10_min, para10_max);
    chain_parm[11] = para11_init(para11_min, para11_max);
    chain_parm[12] = para12_init(para12_min, para12_max);
    //
    ////// set sigma tuning range of parms
    //
    sigma_parm_min = (double *) malloc(sizeof(double)*N_parm);
    sigma_parm_max = (double *) malloc(sizeof(double)*N_parm);
    //
    // initilize
    for (int i=0; i<N_parm; i++)
    {
       sigma_parm_min[i] = 0;
       sigma_parm_max[i] = 0;
    }
    //
    calc_sigma_scale_boundary(sigma_scale_min, sigma_scale_max, sigma_parm_min, sigma_parm_max);
}
//
void init_parm_manual(double *init_loci, double* chain_parm)
{
    for (int i=0; i<N_parm; i++)
    {
        chain_parm[i] = init_loci[i];
    }
    //
    char *input_file;
    input_file = "input.ini";
    // read the range
    read_parm_range(input_file);
    //
       ////// set sigma tuning range of parms
    //
    sigma_parm_min = (double *) malloc(sizeof(double)*N_parm);
    sigma_parm_max = (double *) malloc(sizeof(double)*N_parm);
    calc_sigma_scale_boundary(sigma_scale_min, sigma_scale_max, sigma_parm_min, sigma_parm_max);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///////   Part 5: calc the scaling boundary of each para
////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////
int calc_sigma_scale_boundary(double sigma_scale_min, double sigma_scale_max, double *sigma_parm_min, double *sigma_parm_max)
{
    //
    int sigma_bound_ok = 0;
    
    sigma_parm_min[0]= sigma_scale_min * (para0_max - para0_min);
    sigma_parm_min[1]= sigma_scale_min * (para1_max - para1_min);
    sigma_parm_min[2]= sigma_scale_min * (para2_max - para2_min);
    sigma_parm_min[3]= sigma_scale_min * (para3_max - para3_min);
    sigma_parm_min[4]= sigma_scale_min * (para4_max - para4_min);
    sigma_parm_min[5]= sigma_scale_min * (para5_max - para5_min);
    sigma_parm_min[6]= sigma_scale_min * (para6_max - para6_min);
    sigma_parm_min[7]= sigma_scale_min * (para7_max - para7_min);
    sigma_parm_min[8]= sigma_scale_min * (para8_max - para8_min);
    sigma_parm_min[9]= sigma_scale_min * (para9_max - para9_min);
    sigma_parm_min[10]= sigma_scale_min * (para10_max - para10_min);
    sigma_parm_min[11]= sigma_scale_min * (para11_max - para11_min);
    sigma_parm_min[12]= sigma_scale_min * (para12_max - para12_min);
    sigma_parm_max[0]= sigma_scale_max * (para0_max - para0_min);
    sigma_parm_max[1]= sigma_scale_max * (para1_max - para1_min);
    sigma_parm_max[2]= sigma_scale_max * (para2_max - para2_min);
    sigma_parm_max[3]= sigma_scale_max * (para3_max - para3_min);
    sigma_parm_max[4]= sigma_scale_max * (para4_max - para4_min);
    sigma_parm_max[5]= sigma_scale_max * (para5_max - para5_min);
    sigma_parm_max[6]= sigma_scale_max * (para6_max - para6_min);
    sigma_parm_max[7]= sigma_scale_max * (para7_max - para7_min);
    sigma_parm_max[8]= sigma_scale_max * (para8_max - para8_min);
    sigma_parm_max[9]= sigma_scale_max * (para9_max - para9_min);
    sigma_parm_max[10]= sigma_scale_max * (para10_max - para10_min);
    sigma_parm_max[11]= sigma_scale_max * (para11_max - para11_min);
    sigma_parm_max[12]= sigma_scale_max * (para12_max - para12_min);

    //
    for (int i=0; i<N_parm; i++)
    {
        // each of the following case is a error, sum how many bad sigma
        if (sigma_parm_min[i] < 0)
        {
            sigma_bound_ok++;
        }
        //
        if (sigma_parm_max[i] < 0)
        {
            sigma_bound_ok++;
        }
        if ((sigma_parm_max[i] - sigma_parm_min[i]) < 0)
        {
            sigma_bound_ok++;
        }
    }
    //
    ///////
    //
    if (sigma_bound_ok > 0)
    {
        printf("ERR: bad sigma_parm_min/max TIMES!\n");
        return 1;
    }
    else
    {
        return 0;
    }
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///////   Part 6: init gaussian proposal
///////
////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////
//
// set init gaussian proposal for the sampling
int init_gaussian_proposal(double *ptr_sigma_prop, double init_gp_ratio)
{
    //
    // PingPong outilers inside the allowed prior range
        //
    double para0_gp;
    double para1_gp;
    double para2_gp;
    double para3_gp;
    double para4_gp;
    double para5_gp;
    double para6_gp;
    double para7_gp;
    double para8_gp;
    double para9_gp;
    double para10_gp;
    double para11_gp;
    double para12_gp;
    //
    para0_gp = (para0_max - para0_min) * init_gp_ratio;
    para1_gp = (para1_max - para1_min) * init_gp_ratio;
    para2_gp = (para2_max - para2_min) * init_gp_ratio;
    para3_gp = (para3_max - para3_min) * init_gp_ratio;
    para4_gp = (para4_max - para4_min) * init_gp_ratio;
    para5_gp = (para5_max - para5_min) * init_gp_ratio;
    para6_gp = (para6_max - para6_min) * init_gp_ratio;
    para7_gp = (para7_max - para7_min) * init_gp_ratio;
    para8_gp = (para8_max - para8_min) * init_gp_ratio;
    para9_gp = (para9_max - para9_min) * init_gp_ratio;
    para10_gp = (para10_max - para10_min) * init_gp_ratio;
    para11_gp = (para11_max - para11_min) * init_gp_ratio;
    para12_gp = (para12_max - para12_min) * init_gp_ratio;
    //
    *(ptr_sigma_prop+0) = para0_gp;
    *(ptr_sigma_prop+1) = para1_gp;
    *(ptr_sigma_prop+2) = para2_gp;
    *(ptr_sigma_prop+3) = para3_gp;
    *(ptr_sigma_prop+4) = para4_gp;
    *(ptr_sigma_prop+5) = para5_gp;
    *(ptr_sigma_prop+6) = para6_gp;
    *(ptr_sigma_prop+7) = para7_gp;
    *(ptr_sigma_prop+8) = para8_gp;
    *(ptr_sigma_prop+9) = para9_gp;
    *(ptr_sigma_prop+10) = para10_gp;
    *(ptr_sigma_prop+11) = para11_gp;
    *(ptr_sigma_prop+12) = para12_gp;

    //
    return 0; 
    //
}
//
int init_gaussian_proposal_manual(double *ptr_sigma_prop, double *init_step, double *scale_step_beta, int my_rank)
{
    *(ptr_sigma_prop+0) = (para0_max - para0_min) * init_step[0] * scale_step_beta[my_rank];
    *(ptr_sigma_prop+1) = (para1_max - para1_min) * init_step[1] * scale_step_beta[my_rank];
    *(ptr_sigma_prop+2) = (para2_max - para2_min) * init_step[2] * scale_step_beta[my_rank];
    *(ptr_sigma_prop+3) = (para3_max - para3_min) * init_step[3] * scale_step_beta[my_rank];
    *(ptr_sigma_prop+4) = (para4_max - para4_min) * init_step[4] * scale_step_beta[my_rank];
    *(ptr_sigma_prop+5) = (para5_max - para5_min) * init_step[5] * scale_step_beta[my_rank];
    *(ptr_sigma_prop+6) = (para6_max - para6_min) * init_step[6] * scale_step_beta[my_rank];
    *(ptr_sigma_prop+7) = (para7_max - para7_min) * init_step[7] * scale_step_beta[my_rank];
    *(ptr_sigma_prop+8) = (para8_max - para8_min) * init_step[8] * scale_step_beta[my_rank];
    *(ptr_sigma_prop+9) = (para9_max - para9_min) * init_step[9] * scale_step_beta[my_rank];
    *(ptr_sigma_prop+10) = (para10_max - para10_min) * init_step[10] * scale_step_beta[my_rank];
    *(ptr_sigma_prop+11) = (para11_max - para11_min) * init_step[11] * scale_step_beta[my_rank];
    *(ptr_sigma_prop+12) = (para12_max - para12_min) * init_step[12] * scale_step_beta[my_rank];

    //
    return 0;
    //
}
   
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///////   Part 6: check if a proposed point is within the range
///////
////////////////////////////////////////////
////////////////////////////////////////////
//
// TODO: add a new function for this ugly long boring function.
//
int para_boundary(double *ptr_one_chain_new)
{
    //
    // PingPong outilers inside the allowed prior range
    //
    double para0;
    double para1;
    double para2;
    double para3;
    double para4;
    double para5;
    double para6;
    double para7;
    double para8;
    double para9;
    double para10;
    double para11;
    double para12;
    //
    para0 = *(ptr_one_chain_new+0);
    para1 = *(ptr_one_chain_new+1);
    para2 = *(ptr_one_chain_new+2);
    para3 = *(ptr_one_chain_new+3);
    para4 = *(ptr_one_chain_new+4);
    para5 = *(ptr_one_chain_new+5);
    para6 = *(ptr_one_chain_new+6);
    para7 = *(ptr_one_chain_new+7);
    para8 = *(ptr_one_chain_new+8);
    para9 = *(ptr_one_chain_new+9);
    para10 = *(ptr_one_chain_new+10);
    para11 = *(ptr_one_chain_new+11);
    para12 = *(ptr_one_chain_new+12);
    //
    *(ptr_one_chain_new+0) = bounce_inside(para0, para0_min, para0_max);
    //
    *(ptr_one_chain_new+1) = bounce_inside(para1, para1_min, para1_max);
    //
    *(ptr_one_chain_new+2) = bounce_inside(para2, para2_min, para2_max);
    //
    *(ptr_one_chain_new+3) = bounce_inside(para3, para3_min, para3_max);
    //
    *(ptr_one_chain_new+4) = bounce_inside(para4, para4_min, para4_max);
    //
    *(ptr_one_chain_new+5) = bounce_inside(para5, para5_min, para5_max);
    //
    *(ptr_one_chain_new+6) = bounce_inside(para6, para6_min, para6_max);
    //
    *(ptr_one_chain_new+7) = bounce_inside(para7, para7_min, para7_max);
    //
    *(ptr_one_chain_new+8) = bounce_inside(para8, para8_min, para8_max);
    //
    *(ptr_one_chain_new+9) = bounce_inside(para9, para9_min, para9_max);
    //
    *(ptr_one_chain_new+10) = bounce_inside(para10, para10_min, para10_max);
    //
    *(ptr_one_chain_new+11) = bounce_inside(para11, para11_min, para11_max);
    //
    *(ptr_one_chain_new+12) = bounce_inside(para12, para12_min, para12_max);
    //

    //
    return 0;
}
       
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///////   Part 7: log prior
///////   NOTE: Write your own prior function of all parameters
///////
////////////////////////////////////////////
////////////////////////////////////////////
//
// NOTE: para_boundary ensures: ((a<=a_max) && (a>=a_min))
//
////////////////////////////////////////////
//

/////////////////////////
double prior_para0(double para0_min, double para0_max)
{
    return 1/(para0_max-para0_min);
}
//
double prior_para1(double para1_min, double para1_max)
{
    return 1/(para1_max-para1_min);
}
//
double prior_para2(double para2_min, double para2_max)
{
    return 1/(para2_max-para2_min);
}
//
double prior_para3(double para3_min, double para3_max)
{
    return 1/(para3_max-para3_min);
}
//
double prior_para4(double para4_min, double para4_max)
{
    return 1/(para4_max-para4_min);
}
//
double prior_para5(double para5_min, double para5_max)
{
    return 1/(para5_max-para5_min);
}
//
double prior_para6(double para6_min, double para6_max)
{
    return 1/(para6_max-para6_min);
}
//
double prior_para7(double para7_min, double para7_max)
{
    return 1/(para7_max-para7_min);
}
//
double prior_para8(double para8_min, double para8_max)
{
    return 1/(para8_max-para8_min);
}
//
double prior_para9(double para9_min, double para9_max)
{
    return 1/(para9_max-para9_min);
}
//
double prior_para10(double para10_min, double para10_max)
{
    return 1/(para10_max-para10_min);
}
//
double prior_para11(double para11_min, double para11_max)
{
    return 1/(para11_max-para11_min);
}
//
double prior_para12(double para12_min, double para12_max)
{
    return 1/(para12_max-para12_min);
}
//

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/////// Part 8: Combine all prior distributions 
///////
////////////////////////////////////////////
//
double log_prior(double *ptr_one_chain)
{
    //
    double para0;
    double para1;
    double para2;
    double para3;
    double para4;
    double para5;
    double para6;
    double para7;
    double para8;
    double para9;
    double para10;
    double para11;
    double para12;
    //

    //
    double log_prior;
    //
    para0 = *(ptr_one_chain+0);
    para1 = *(ptr_one_chain+1);
    para2 = *(ptr_one_chain+2);
    para3 = *(ptr_one_chain+3);
    para4 = *(ptr_one_chain+4);
    para5 = *(ptr_one_chain+5);
    para6 = *(ptr_one_chain+6);
    para7 = *(ptr_one_chain+7);
    para8 = *(ptr_one_chain+8);
    para9 = *(ptr_one_chain+9);
    para10 = *(ptr_one_chain+10);
    para11 = *(ptr_one_chain+11);
    para12 = *(ptr_one_chain+12);
    //

    log_prior = 0;
                /*
                log(prior_para0(para0_min, para0_max)) +
                log(prior_para1(para1_min, para1_max)) +
                log(prior_para2(para2_min, para2_max)) +
                log(prior_para3(para3_min, para3_max)) +
                log(prior_para4(para4_min, para4_max)) +
                log(prior_para5(para5_min, para5_max)) +
                log(prior_para6(para6_min, para6_max)) +
                log(prior_para7(para7_min, para7_max)) +
                log(prior_para8(para8_min, para8_max)) +
                log(prior_para9(para9_min, para9_max)) +
                log(prior_para10(para10_min, para10_max)) +
                log(prior_para11(para11_min, para11_max)) +
                log(prior_para12(para12_min, para12_max)) +
                +0;    //
		*/

    //
    return log_prior;
    // 
}

////////////////////////////////////////////
////////////////////////////////////////////
/////  Save para debug
/////  NO need to change
////////////////////////////////////////////
////////////////////////////////////////////
//
//
int save_debug_para_boundary(char *p_s, double p, double p_min, double p_max)
{
    FILE *out;
    //
    // set fnames
    char fname[100];

    snprintf(fname, sizeof fname, "%s%s%s%s", results_dir, "/", p_s, ".debug_para_boundary");
    //printf("%s\n", fname);
    //
    // open files
    if ((out = fopen(fname, "a")) == NULL)
    {
        fprintf(stderr, "Can't create output file!\n");
        exit(3);
    }
    //
    fprintf(out, "%f", p);
    fprintf(out, "  ");
    fprintf(out, "%f", p_min);
    fprintf(out, "  ");
    fprintf(out, "%f", p_max);
    fprintf(out, "\n");
    //
    // close files
    if (fclose(out) != 0)
        fprintf(stderr, "Error in closing file!\n");
    //
    return 0;
}
