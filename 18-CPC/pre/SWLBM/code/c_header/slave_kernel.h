#ifndef SLAVE_KERNEL_H
#define SLAVE_KERNEL_H

#include "c_public_const.h"

struct lbm_init_param
{
    int my_id;
    long* host_flag, *slave_flag;

    int iter;

    int x_sec, y_sec, Z;

    int* walls;
    int* flags;

    // #ifndef Real
    // #define Real float
    // #endif

    float nu, omega, CSmago;

    //不会提供nodes，因为其本身就是轮换数组
    //所以会通过flag进行指针传递
//    lbm_flag_type* walls;//not 0 for wallsdd
//    lbm_flag_type* flags;//0 for fluid, 3 for bounce

//    int *halo_xs, *halo_ys;
//    int *inner_xs, *inner_ys;

    //如果后面使用了insane优化，那么还会有这些
//    int *insane_halo_xs, *insane_halo_ys;
//    int *insane_inner_xs, *insane_inner_ys;
};

void cpe_athread_daemon(void *_param);

#define EXIT_FLAG        12450
#define HANDSHAKE_FLAG   19260817  //时间之握
#define STD_LBM_FLAG     114514

#define FLAG_SIZE        32
#define MPI_RANK         1
#define KERNEL_ACTION    2 // work and exit

//STD_LBM
#define STD_POINT_CNT        3
#define STD_XS_PTR           4
#define STD_YS_PTR           5
#define STD_CURRENT_HEAD     6
#define STD_OTHER_HEAD       7
#define INSANE_POINT_CNT     8
#define INSANE_XS_PTR        9
#define INSANE_YS_PTR        10

//std lbm
//#define GROUP_SIZE       3
//#define REMAIN_POINT     4
//#define IN_PTR           5
//#define OUT_PTR          6
//#define IN_STRIDE        7
//#define OUT_STRIDE       8
//#define REQUIRE_IN       9
//#define REQUIRE_OUT      10

//step 40, 3 for rem, 39 for step rem, and extra 6 for safety
#define SLAVE_SAFE_PAD 48

#define STEP_SIZE 50
#define INSANE_SIZE 20

#define BLOCK_SIZE(block_id, total_blocks, n) (((n) / (total_blocks)) + (((n) % (total_blocks) > (block_id)) ? 1 : 0))
#define BLOCK_LOW(block_id, total_blocks, n) (((n) / (total_blocks)) * (block_id) + (((n) % (total_blocks) > (block_id)) ? (block_id) : (n) % (total_blocks)))

#define CPE_TOTAL_SYNC 0x0000FFFF

#endif
