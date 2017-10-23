#ifndef C_PUBLIC_VAR_H
#define C_PUBLIC_VAR_H

#include <stdio.h>
#include <memory.h>
#include <stdlib.h>

#include "c_public_const.h"
#include "slave_kernel.h"

extern struct weno_init_param param_weno;

void c_param_init_(int *my_id,
                   int* nx, int* ny, int *nz,
                   int *npx, int *npy, int *npz,
                   double *dt, double *end_time, double *hx,
                   double *u, double *u1, int *iperiodic, double *tt,
                    int *istep, int *k_step_show, int *k_step_save,
                    int *DAT);

void c_finalize_();

extern double* hj;
extern double* c_u;
extern double* c_u1;
extern double* c_cur_u;
extern double* c_last_u;
extern int sz_c_u;
extern int iter_step;

//#define f(dim1, dim2, dim3) f[dim1 - f_s1 + (dim2 - f_s2) * f_dim2 + (dim3 - f_s3) * f_dim3]
//#define fx(dim1, dim2, dim3) fx[dim1 - fx_s1 + (dim2 - fx_s2) * fx_dim2 + (dim3 - fx_s3) * fx_dim3]

#define u_index(dim1, dim2, dim3) (dim1 - f_s1 + (dim2 - f_s2) * f_dim2 + (dim3 - f_s3) * f_dim3)
#define u1_index(dim1, dim2, dim3) (dim1 - fx_s1 + (dim2 - fx_s2) * fx_dim2 + (dim3 - fx_s3) * fx_dim3)

//#define hj(dim1) hj[dim1]

extern volatile long host_flag[FLAG_SIZE];
extern volatile long slave_flag[FLAG_SIZE];
extern volatile long local_cc[FLAG_SIZE];
extern volatile long flag_to_wait;

void wait_slave_flag();

static inline unsigned long rpcc()
{
    unsigned long time;
    asm("rtc %0": "=r" (time) : );
    return time;
}

#endif
