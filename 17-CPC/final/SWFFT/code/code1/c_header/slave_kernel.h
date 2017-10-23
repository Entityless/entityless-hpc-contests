#ifndef SLAVE_KERNEL_H
#define SLAVE_KERNEL_H

#include "c_public_const.h"

struct fft_init_param
{
    int my_id;
    long* host_flag, *slave_flag;

    int iter;

    //ptr on mpe
    double* u0_real, *u0_imag;//用于读入，包括evv的初始值
    double* u1_real, *u1_imag;//用于写出
    double* u2_real, *u2_imag;//用于读入
    double* twiddle;
    double* u_real, *u_imag;


    dfc* buf_chk;
    int* chk_cnt, *chk_stp, *chk_xs, *chk_ys, *chk_zs;

    int nx, ny, nz;
    int npc;

};

void cpe_athread_daemon(void *_param);

#define EXIT_FLAG       12450
#define ALLIN_FFT_FLAG  2333

#define FLAG_SIZE        32
#define MPI_RANK         1
#define KERNEL_ACTION    2
#define GROUP_SIZE       3
#define REMAIN_POINT     4
#define IN_PTR           5
#define OUT_PTR          6
#define IN_STRIDE        7
#define OUT_STRIDE       8
#define REQUIRE_IN       9
#define REQUIRE_OUT      10

//step 40, 3 for rem, 39 for step rem, and extra 6 for safety
#define SLAVE_SAFE_PAD 48


#define CPE_TOTAL_SYNC 0x0000FFFF

#endif
