#ifndef C_PUBLIC_VAR_H
#define C_PUBLIC_VAR_H

#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include "mpi.h"

#include "c_public_const.h"
#include "slave_kernel.h"
#include <assert.h>
#include <memory.h>

//extern struct weno_init_param param_weno;


void c_param_init_(struct dcomplex* _u0, struct dcomplex* _u1, struct dcomplex* _u2, double* _twiddle, int* _dims, int* _niter,
                   int* _fftblockpad_default, int* _fftblock, int* _maxdim, int* _layout_type,
                   struct dcomplex* _u, int* _me, int* _np1, int* _np2, int* _ntdivnp,
                   int* _np, int* _transblockpad, int* _transblock, int* _fftblockpad,
                   int* _nx, int* _ny, int* _nz,
                   int* _xstart,int *_ystart,int *_zstart,
                   int* _xend,int *_yend,int *_zend, struct dcomplex* _sums);

void c_finalize_();

extern volatile long host_flag[FLAG_SIZE];
extern volatile long slave_flag[FLAG_SIZE];
extern volatile long local_cc[FLAG_SIZE];
extern volatile long flag_to_wait;

////FFF
//ptr
extern struct dcomplex* f_u0_, *f_u1_, *f_u2_;//f
extern int* dims_;
extern double* f_twiddle_;//f
extern struct dcomplex* f_u_;//f
extern int* xstart_, *ystart_, *zstart_;
extern int* xend_, *yend_, *zend_;
extern struct dcomplex* f_sums_;//f
//val
extern int niter_;
extern int fftblockpad_default_, maxdim_, layout_type_;
extern int fftblock_, fftblockpad_;
extern int me_, np1_, np2_, np_, ntdivnp_;
extern int transblock_, transblockpad_;
extern int nx_, ny_, nz_;
extern double ntotal_f_;
//easy access
#define dims(a, b) dims_[a - 1 + (b - 1) * 3]
//local pre-allocated ptr
extern my_fc* scratch_;
//comm
extern MPI_Comm commslice1_, commslice2_;
//c_local
//extern my_fc* c_u0_, *c_u1_, *c_u2_;
//extern my_ft* c_twiddle_;
//extern dfc* c_u_;
extern dfc* c_sums_;
extern double* c_u0_real_, *c_u0_imag_;
extern double* c_u1_real_, *c_u1_imag_;
extern double* c_u2_real_, *c_u2_imag_;
extern double* c_u_real_, *c_u_imag_;
extern double* c_twiddle_;
//support
extern int chk_cnt_;
extern int* chk_xs_;
extern int* chk_ys_;
extern int* chk_zs_;
extern dfc* buf_chk_;
extern int* chk_core_cnt_;
extern int* chk_core_stp_;


////mpe benchmark
extern unsigned long mpe_cc_cur_[100];
extern unsigned long mpe_cc_total_[100];
#define MPE_EVV 0
#define MPE_CHK 1
#define MPE_FT1 2
#define MPE_FT2 3
#define MPE_FT3 4
#define MPE_T2F 5
#define MPE_T2G 6
#define MPE_T2L 7


void wait_slave_flag();

static inline unsigned long rpcc()
{
    unsigned long time;
    asm("rtc %0": "=r" (time) : );
    return time;
}

static inline double ccts(unsigned long cc)
{
    return cc * 1.0 / CCPS;
}

#endif
