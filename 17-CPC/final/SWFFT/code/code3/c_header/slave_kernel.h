#ifndef SLAVE_KERNEL_H
#define SLAVE_KERNEL_H

#include "c_public_const.h"

// (o゜▽゜)o☆
struct fft_init_param
{
    int my_id;
    long* host_flag, *slave_flag;

    int iter;

    //ptr on mpe
    double* twiddle;
    //double* u_real, *u_imag;
    double* u_real_left, *u_imag_left;
    double* u_real_right, *u_imag_right;
    double last_ureal, last_uimag;


    dfc* buf_chk;
    int* chk_cnt, *chk_stp, *chk_xs, *chk_ys, *chk_zs;

    int nx, ny, nz;
    int npc;
};

void cpe_athread_daemon(void *_param);

#define EXIT_FLAG       12450
//#define ALLIN_FFT_FLAG  2333
//#define EVV_FFT_FLAG 23232
//#define SG_FFT_FLAG 289346
#define CP_EVV_FFT_FLAG 54612
//#define COUTLE_FFT_FLAG 43287
#define DB_FFT_FLAG 438712
//#define TP_LOCAL_FLAG  243653
#define TP_FINISH_FLAG 435432


#define FLAG_SIZE        32
#define MPI_RANK         1
#define KERNEL_ACTION    2
#define GROUP_SIZE       3
#define REMAIN_POINT     4
#define IN_PTR_REAL      5
#define OUT_PTR_REAL     6
#define IN_PTR_IMAG      7
#define OUT_PTR_IMAG     8
#define IN_STRIDE        9
#define OUT_STRIDE       10
#define FFT_D1           11
#define FFT_D2           12
#define FFT_D3           13
#define SUP_PTR_REAL     14
#define SUP_PTR_IMAG     15
#define FFT_D1_FNS       16//for finish
#define FFT_D2_FNS       17
#define FFT_D3_FNS       18

//step 40, 3 for rem, 39 for step rem, and extra 6 for safety
#define SLAVE_SAFE_PAD 48


#define CPE_TOTAL_SYNC 0x0000FFFF

#endif
