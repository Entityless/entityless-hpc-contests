#include "../c_header/c_public_var.h"
#include "../c_header/c_public_const.h"
#include "mpi.h"


struct weno_init_param param_weno;

double* hj;
double* c_u;
double* c_u1;
double* c_cur_u;
double* c_last_u;
int sz_c_u;
int iter_step;

volatile long host_flag[FLAG_SIZE];
volatile long slave_flag[FLAG_SIZE];
volatile long local_cc[FLAG_SIZE];
volatile long flag_to_wait;

void c_param_init_(int* my_id,
                   int* nx, int* ny, int *nz,
                   int *npx, int *npy, int *npz,
                   double *dt, double *end_time, double *hx,
                   double *u, double *u1,
                   int* iperiodic, double* tt, int* istep,
                   int* k_step_show, int *k_step_save,
                   int* LAP)
{
    param_weno.my_id = *my_id;

    printf("process %d enter c_param_init_\n", param_weno.my_id);fflush(stdout);

    param_weno.nx = *nx;
    param_weno.ny = *ny;
    param_weno.nz = *nz;
    param_weno.npx = *npx;
    param_weno.npy = *npy;
    param_weno.npz = *npz;

    param_weno.dt = *dt;
    param_weno.end_time = *end_time;
    param_weno.hx = *hx;

    param_weno.u = u;
    param_weno.u1 = u1;
    param_weno.lap = *LAP;

    param_weno.iperiodic = iperiodic;

    hj = malloc(sizeof(double) * (1 + param_weno.nx));

    sz_c_u = (param_weno.nx + 2 * param_weno.lap) * (param_weno.ny + 2 * param_weno.lap) * (param_weno.nz + 2 * param_weno.lap);
    c_u1 = malloc(sizeof(double) * sz_c_u);
    c_u = malloc(sizeof(double) * sz_c_u);

    printf("process %d malloc, sz = %d\n", param_weno.my_id, sz_c_u);

    param_weno.tt = tt;
    param_weno.istep = istep;

    param_weno.k_step_show = *k_step_show;
    param_weno.k_step_save = *k_step_save;


    memcpy(c_u, u, sizeof(double) * sz_c_u);
    memcpy(c_u1, u, sizeof(double) * sz_c_u);

    c_cur_u = c_u;
    c_last_u = c_u1;

    //touch_arr();

    printf("process %d quit c_param_init_\n", param_weno.my_id);fflush(stdout);
}

void c_finalize_()
{
    free(hj);
}

void wait_slave_flag()
{
    //if(my_rank == 0)
    //    printf("wait flag %ld\n", flag_to_wait);
    while(slave_flag[0] < flag_to_wait);
    //if(my_rank == 0)
    //    printf("get flag %ld\n", flag_to_wait);
    flag_to_wait++;
}
