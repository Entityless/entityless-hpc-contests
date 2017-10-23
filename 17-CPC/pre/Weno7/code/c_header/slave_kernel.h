#ifndef SLAVE_KERNEL_H
#define SLAVE_KERNEL_H

struct weno_init_param
{
    int my_id;

    int nx;
    int ny;
    int nz;
    int npx;
    int npy;
    int npz;
    double dt;
    double end_time;
    double hx;

    double* u;
    double* u1;

    int* iperiodic;

    double* tt;
    int* istep;

    int k_step_show;
    int k_step_save;

    long* host_flag, *slave_flag;

    int lap;
};

void cpe_athread_daemon(void *_param);

#define EXIT_FLAG       12450
#define WENO_INNER_FLAG 23232

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
