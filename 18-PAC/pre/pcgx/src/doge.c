// #include <math.h>
#define _GNU_SOURCE
#include <sched.h>
#include <stdio.h>
#include <unistd.h>
#include <omp.h>
#include <mpi.h>
#include <memory.h>
#include <stddef.h>

//用了2年半的祖传函数
static inline int BLOCK_SIZE(int block_id, int total_blocks, int n)
{
    return (n / total_blocks) + ((n % total_blocks > block_id) ? 1 : 0);
}

static inline int BLOCK_LOW(int block_id, int total_blocks, int n)
{
    return (n / total_blocks) * block_id + ((n % total_blocks > block_id) ? block_id : n % total_blocks);
}

int hch_doge_my_rank, hch_doge_comm_sz;

//下划线表示这个变量来自于F90代码
double *B_, *R_, *S_, *Q_, *WORK0_, *WORK1_, *X_;
const double *A0_, *AN_, *AE_, *ANE_, *RCALCT_B_;
MPI_Comm *dist_comm_;

#define realt float

double local_sum;

//喵喵喵
int full_len, inner_len;//总数组的长度
int i_inner_len, j_inner_len;//;, block_fat, block_vec_cnt, block_vec_remain;
#define VEC_RED_LEN 32 //4 cache line
//还具有提高精度的作用
realt local_sum_red[VEC_RED_LEN] __attribute__ ((aligned(64)));
realt *post_A0;//预处理的A0
realt *post_mask;//预处理的mask

//喵喵喵喵喵喵喵？喵喵喵喵、喵喵喵喵喵喵。
realt *CQ, *CX, *CR;
realt *CS0;//和S_完全相同，只是因为我不知道要怎么办才暂时设立的。
//初始化后只读，A
realt *A00, *AN0, *AE0, *ANE0, *B0;

int nx_block_, ny_block_, max_blocks_tropic_, nblocks_tropic_, solv_ncheck_;

//数组的缩减在后面再进行

int *ibs_, *ies_, *jbs_, *jes_, ib, ie, jb, je;

extern void hch_timer_start(int num);
extern void hch_timer_stop(int num);

//numa opt

//manual pin core

void bind_to_core(int cid)
{
    // cpu_set_t mask;
    // cpu_set_t get;
    // int i, j, k;
    // int num = sysconf(_SC_NPROCESSORS_CONF);
    // //printf("system has %d processor(s)\n", num);

    // int tid = cid;

    // CPU_ZERO_S(sizeof(cpu_set_t), &mask);
    // CPU_SET_S(tid, sizeof(cpu_set_t), &mask);

    // if (pthread_setaffinity_np(pthread_self(), sizeof(mask), &mask) < 0) {
    //     fprintf(stderr, "set thread affinity failed\n");
    // }

    // CPU_ZERO_S(sizeof(cpu_set_t), &get);
    // if (pthread_getaffinity_np(pthread_self(), sizeof(get), &get) < 0) {
    //     fprintf(stderr, "get thread affinity failed\n");
    // }

    // for (j = 0; j < num; j++) {
    //     if (CPU_ISSET_S(j, sizeof(cpu_set_t), &get)) {
    //         // printf("thread %d is running in processor %d\n", (int)pthread_self(), j);
    //     }
    // }
    kmp_affinity_mask_t mask;
    int err;
    kmp_create_affinity_mask(&mask);
    printf("rank = %d, mask created\n", hch_doge_my_rank);
    err = kmp_set_affinity_mask_proc (cid, &mask);
    printf("rank = %d, set mask err = %d\n", hch_doge_my_rank, err);
    err = kmp_set_affinity(mask);
    printf("rank = %d, set affinity err = %d\n", hch_doge_my_rank, err);

}

#define HCH_CPU_CNT 2

int cpu_slot_start[HCH_CPU_CNT];
int cpu_slot_size[HCH_CPU_CNT];
int my_core;
int my_cpu;
int total_core;
int core_per_cpu;
int* core_proc_start;
int* core_proc_size;

int numa_init_flag = 0;
//how
void numa_init()
{
    if(numa_init_flag == 0)
    {
        //the content

        int nproc = hch_doge_comm_sz;
        total_core = sysconf(_SC_NPROCESSORS_CONF);
        core_per_cpu = total_core / HCH_CPU_CNT;
        if(hch_doge_my_rank == 0)
            printf("total core = %d, core per cpu = %d\n", total_core, core_per_cpu);

        core_proc_start = (int*)malloc(sizeof(int) * total_core);
        core_proc_size = (int*)malloc(sizeof(int) * total_core);

        int i, local_core;
        for(i = 0; i < HCH_CPU_CNT; i++)
        {
            cpu_slot_start[i] = BLOCK_LOW(i, HCH_CPU_CNT, nproc);
            cpu_slot_size [i] = BLOCK_SIZE(i, HCH_CPU_CNT, nproc);

            if(hch_doge_my_rank == 0)
                printf("i = %d, cpu_slot_start = %d, cpu_slot_size = %d\n",
                        i, cpu_slot_start[i], cpu_slot_size[i]);

            //locate my cpu
            if(hch_doge_my_rank >= cpu_slot_start[i] && hch_doge_my_rank < cpu_slot_start[i] + cpu_slot_size[i])
                my_cpu = i;

            //initial core data
            for(local_core = 0; local_core < core_per_cpu; local_core++)
            {
                int actual_cid = local_core + i * core_per_cpu;

                core_proc_start[actual_cid] = cpu_slot_start[i] + BLOCK_LOW(local_core, core_per_cpu, cpu_slot_size[i]);
                core_proc_size[actual_cid] = BLOCK_SIZE(local_core, core_per_cpu, cpu_slot_size[i]);

                if(hch_doge_my_rank == 0)
                    printf("actual_cid = %d, core_proc_start = %d, core_proc_size = %d\n",
                            actual_cid, core_proc_start[actual_cid], core_proc_size[actual_cid]);

                if(hch_doge_my_rank >= core_proc_start[actual_cid] && hch_doge_my_rank < core_proc_start[actual_cid] + core_proc_size[actual_cid])
                    my_core = actual_cid;
            }
            //init proc data
        }

        numa_init_flag = 1;
    }
}

void pcg_iter_init_(int *_nx_block, int *_ny_block, int *_max_blocks_tropic, int *_nblocks_tropic, int *_solv_ncheck,// int* local_block
                    double *_X, double *_B, double *_R, realt *_CR, double *_S, double *_Q, double *_WORK0, double *_WORK1,
                    const double *_A0, const double *_AN, const double *_AE, const double *_ANE, const double *_RCALCT_B, MPI_Comm* _dist_comm)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &hch_doge_my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &hch_doge_comm_sz);

    // numa_init();

    if(hch_doge_my_rank == 999)
    {
        printf("pcg_iter_init_ rank = %d, %d, %d, %d, %d, %d\n", hch_doge_my_rank, *_nx_block, *_ny_block, *_max_blocks_tropic, *_nblocks_tropic, *_solv_ncheck);
    }

    //val copy
    nx_block_ = *_nx_block;
    ny_block_ = *_ny_block;
    max_blocks_tropic_ = *_max_blocks_tropic;
    nblocks_tropic_ = *_nblocks_tropic;
    solv_ncheck_ = *_solv_ncheck;
    RCALCT_B_ = _RCALCT_B;

    //var ptr copy
    B_ = _B;
    R_ = _R;
    CR = _CR;
    S_ = _S;
    Q_ = _Q;
    WORK0_ = _WORK0;
    WORK1_ = _WORK1;

    //const ptr copy
    X_ = _X;
    A0_ = _A0;
    AE_ = _AE;
    AN_ = _AN;
    ANE_ = _ANE;

    ibs_ = malloc(sizeof(int) * nblocks_tropic_);
    ies_ = malloc(sizeof(int) * nblocks_tropic_);
    jbs_ = malloc(sizeof(int) * nblocks_tropic_);
    jes_ = malloc(sizeof(int) * nblocks_tropic_);

    full_len = nblocks_tropic_ * nx_block_ * ny_block_;
    CS0 = (realt*)_mm_malloc(sizeof(realt) * full_len, 2*1024*1024);
    A00 = (realt*)_mm_malloc(sizeof(realt) * full_len, 2*1024*1024);
    AN0 = (realt*)_mm_malloc(sizeof(realt) * full_len, 2*1024*1024);
    AE0 = (realt*)_mm_malloc(sizeof(realt) * full_len, 2*1024*1024);
    ANE0 = (realt*)_mm_malloc(sizeof(realt) * full_len, 2*1024*1024);
    B0 = (realt*)_mm_malloc(sizeof(realt) * full_len, 2*1024*1024);

    post_A0 = (realt*)_mm_malloc(sizeof(realt) * full_len, 2*1024*1024);
    CX = (realt*)_mm_malloc(sizeof(realt) * full_len, 2*1024*1024);
    // CR = (realt*)_mm_malloc(sizeof(realt) * full_len, 2*1024*1024);
    post_mask = (realt*)_mm_malloc(sizeof(realt) * full_len, 2*1024*1024);
    CQ = (realt*)_mm_malloc(sizeof(realt) * full_len, 2*1024*1024);

    dist_comm_ = _dist_comm;
}

void pcg_iter_finalize_()
{
    free(ibs_);
    free(ies_);
    free(jbs_);
    free(jes_);
    _mm_free(post_A0);
    _mm_free(post_mask);
}

void ijbe_init_(int* pos, int* ib, int* ie, int* jb, int* je)
{
    ibs_[*pos] = *ib;
    ies_[*pos] = *ie;
    jbs_[*pos] = *jb;
    jes_[*pos] = *je;
}

//第一次跑16那个循环。预判断并保存==0的索引；预计算A0的倒数从而减少除法开销。
//WORK里面是0的地方，总是0，也就是说只需要赋值一次。
void loop_16_first_()
{
    //the inner initial
    //proved that *bs and *es won't change in this code
    i_inner_len = ies_[0] - ibs_[0] + 1;
    j_inner_len = jes_[0] - jbs_[0] + 1;
    ib = ibs_[0];
    ie = ies_[0];
    jb = jbs_[0];
    je = jes_[0];

    inner_len = i_inner_len * j_inner_len * nblocks_tropic_;

    int iblock, d1, d2, d3;
    for(iblock = 0; iblock < nblocks_tropic_; iblock++)
    {
        int idx_d3 = ny_block_ * nx_block_ * iblock;
        for(d2 = jb - 1; d2 <= je - 1; d2++)
        // for(d2 = 0; d2 < ny_block_; d2++)
        {
            int idx_d2 = idx_d3 + nx_block_ * d2;
            for(d1 = ib - 1; d1 <= ie - 1; d1++)
            // for(d1 = 0; d1 < nx_block_; d1++)
            {
                int idx = idx_d2 + d1;


            }
        }
        for(d2 = 0; d2 < ny_block_; d2++)
        {
            int idx_d2 = idx_d3 + nx_block_ * d2;
            for(d1 = 0; d1 < nx_block_; d1++)
            {
                int idx = idx_d2 + d1;

                if(A0_[idx] != 0.0)
                {
                    post_A0[idx] = 1.0 / A0_[idx];
                }
                else
                {
                    post_A0[idx] = 0.0;
                }

                CS0[idx] = S_[idx];
                A00[idx] = A0_[idx];
                AN0[idx] = AN_[idx];
                AE0[idx] = AE_[idx];
                ANE0[idx] = ANE_[idx];
                B0[idx] = B_[idx];
                post_mask[idx] = RCALCT_B_[idx];
                CX[idx] = X_[idx];
                CR[idx] = R_[idx];
            }
        }
    }
}

//还可以继续优化
void loop_16_()
{
    int iblock, d1, d2, d3;
    memset(local_sum_red, 0, sizeof(realt) * VEC_RED_LEN);
    local_sum = 0.0;
    int vec_cnt = i_inner_len / VEC_RED_LEN;

    for(iblock = 0; iblock < nblocks_tropic_; iblock++)
    {
        int idx_d3 = ny_block_ * nx_block_ * iblock;
        // for(d2 = 0; d2 < ny_block_; d2++)
        for(d2 = jb - 1; d2 <= je - 1; d2++)
        {
            int idx_d2 = idx_d3 + nx_block_ * d2;
            // for(d1 = 0; d1 < nx_block_; d1++)
            for(d1 = ib - 1; d1 <= ie - 1; d1++)
            {
                int idx = idx_d2 + d1;

                CQ[idx] = CR[idx] * post_A0[idx];
            }

            for(d1 = 0; d1 < vec_cnt * VEC_RED_LEN; d1 += VEC_RED_LEN)
            {
                int idx = idx_d2 + d1 + ib - 1;
                int d0;
                for(d0 = 0; d0 < VEC_RED_LEN; d0++)
                {
                    local_sum_red[d0] += CQ[idx + d0] * CR[idx + d0] * post_mask[idx + d0];
                }
            }
            for(d1 = vec_cnt * VEC_RED_LEN; d1 < i_inner_len; d1++)
            {
                int idx = idx_d2 + d1 + ib - 1;
                local_sum_red[d1 - vec_cnt * VEC_RED_LEN] += CQ[idx] * CR[idx] * post_mask[idx];
            }
        }

    }
    for(d1 = 0; d1 < VEC_RED_LEN; d1++)
    {
        local_sum += local_sum_red[d1];
    }
}

void regular_masked_global_sum_dbl0_(double* ret)
{
    MPI_Allreduce(&local_sum, ret, 1, MPI_DOUBLE, MPI_SUM, *dist_comm_);
}

void loop_17_(double* eta0, double* eta1)
{
    int iblock, d1, d2, d3;
    local_sum = 0.0;
    memset(local_sum_red, 0, sizeof(realt) * VEC_RED_LEN);
    const double eta_div = (*eta1 / *eta0);

    for(iblock = 0; iblock < nblocks_tropic_; iblock++)
    {
        int idx_d3 = ny_block_ * nx_block_ * iblock;
        for(d2 = jb - 1; d2 <= je - 1; d2++)
        {
            int idx_d2 = idx_d3 + nx_block_ * d2;
            for(d1 = ib - 1; d1 <= ie - 1; d1++)
            {
                int idx = idx_d2 + d1;

                CS0[idx] = CQ[idx] + CS0[idx] * eta_div;
            }
        }

        for(d2 = jb - 1; d2 <= je - 1; d2++)
        {
            int idx_d2 = idx_d3 + nx_block_ * d2;

            int idx00 = idx_d2 + ib - 1;//idx i j
            //p -- plus 1, m -- minus 1
            int idx0p = idx00 + nx_block_;
            int idx0m = idx00 - nx_block_;
            int idxm0 = idx00 - 1;
            int idxmp = idxm0 + nx_block_;
            int idxmm = idxm0 - nx_block_;
            int idxp0 = idx00 + 1;
            int idxpp = idxp0 + nx_block_;
            int idxpm = idxp0 - nx_block_;

            #define CONV_UNROLL 16
            //这个数值一定要小于VEC_RED_LEN

            int unroll_cnt = i_inner_len / CONV_UNROLL;
            // unroll_cnt = 0;

            realt A00_red[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt AN0_00_red[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt AN0_0m_red[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt AE0_00_red[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt AE0_m0_red[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt ANE0_00_red[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt ANE0_0m_red[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt ANE0_m0_red[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt ANE0_mm_red[CONV_UNROLL] __attribute__ ((aligned(64)));

            realt CS_00[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt CS_0p[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt CS_0m[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt CS_p0[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt CS_m0[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt CS_pp[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt CS_pm[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt CS_mp[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt CS_mm[CONV_UNROLL] __attribute__ ((aligned(64)));

            realt conv_q[CONV_UNROLL] __attribute__ ((aligned(64)));

            for(d1 = 0; d1 < unroll_cnt * CONV_UNROLL; d1 += CONV_UNROLL)
            {
                // int idx_plus = d1 * CONV_UNROLL;//大错特错
                int d0;

                //read
                for(d0 = 0; d0 < CONV_UNROLL; d0++)
                {
                    A00_red    [d0] = A00 [idx00 + d1 + d0];
                    AN0_00_red [d0] = AN0 [idx00 + d1 + d0];
                    AN0_0m_red [d0] = AN0 [idx0m + d1 + d0];
                    AE0_00_red [d0] = AE0 [idx00 + d1 + d0];
                    AE0_m0_red [d0] = AE0 [idxm0 + d1 + d0];
                    ANE0_00_red[d0] = ANE0[idx00 + d1 + d0];
                    ANE0_0m_red[d0] = ANE0[idx0m + d1 + d0];
                    ANE0_m0_red[d0] = ANE0[idxm0 + d1 + d0];
                    ANE0_mm_red[d0] = ANE0[idxmm + d1 + d0];

                    CS_00[d0]  = CS0[idx00 + d1 + d0];
                    CS_0p[d0]  = CS0[idx0p + d1 + d0];
                    CS_0m[d0]  = CS0[idx0m + d1 + d0];
                    CS_p0[d0]  = CS0[idxp0 + d1 + d0];
                    CS_m0[d0]  = CS0[idxm0 + d1 + d0];
                    CS_pp[d0]  = CS0[idxpp + d1 + d0];
                    CS_pm[d0]  = CS0[idxpm + d1 + d0];
                    CS_mp[d0]  = CS0[idxmp + d1 + d0];
                    CS_mm[d0]  = CS0[idxmm + d1 + d0];
                }

                //conv
                for(d0 = 0; d0 < CONV_UNROLL; d0++)
                {
                    conv_q[d0] = A00_red    [d0] * CS_00[d0] +
                                 AN0_00_red [d0] * CS_0p[d0] +
                                 AN0_0m_red [d0] * CS_0m[d0] +
                                 AE0_00_red [d0] * CS_p0[d0] +
                                 AE0_m0_red [d0] * CS_m0[d0] +
                                 ANE0_00_red[d0] * CS_pp[d0] +
                                 ANE0_0m_red[d0] * CS_pm[d0] +
                                 ANE0_m0_red[d0] * CS_mp[d0] +
                                 ANE0_mm_red[d0] * CS_mm[d0]; 
                }

                //store and local_sum
                for(d0 = 0; d0 < CONV_UNROLL; d0++)
                {
                    int idx = idx_d2 + (d1 + ib - 1) + d0;
                    CQ[idx] = conv_q[d0];
                    local_sum_red[d0] += CQ[idx] * CS_00[d0] * post_mask[idx];
                }
            }


            for(d1 = unroll_cnt * CONV_UNROLL; d1 < i_inner_len; d1++)
            {
                int idx = idx_d2 + (d1 + ib - 1);
                CQ[idx] = A00 [idx00 + d1] * CS0[idx00 + d1] +
                                 AN0 [idx00 + d1] * CS0[idx0p + d1] +
                                 AN0 [idx0m + d1] * CS0[idx0m + d1] +
                                 AE0 [idx00 + d1] * CS0[idxp0 + d1] +
                                 AE0 [idxm0 + d1] * CS0[idxm0 + d1] +
                                 ANE0[idx00 + d1] * CS0[idxpp + d1] +
                                 ANE0[idx0m + d1] * CS0[idxpm + d1] +
                                 ANE0[idxm0 + d1] * CS0[idxmp + d1] +
                                 ANE0[idxmm + d1] * CS0[idxmm + d1];  //AX, X全是参数传入

                //由于CONV_UNROLL一定小于VEC_RED_LEN，可以使用下面的方法
                //无法在测试环境中使用
                // local_sum_red[d1 - unroll_cnt * CONV_UNROLL] += CQ[idx] * CS0[idx00 + d1] * post_mask[idx];
                local_sum_red[0] += CQ[idx] * CS0[idx00 + d1] * post_mask[idx];
            }
        }
    }

    for(d1 = 0; d1 < VEC_RED_LEN; d1++)
    {
        local_sum += local_sum_red[d1];
    }
}

void loop_18_(double* eta1)
{
    int iblock, d1, d2, d3;

    for(iblock = 0; iblock < nblocks_tropic_; iblock++)
    {
        int idx_d3 = ny_block_ * nx_block_ * iblock;
        for(d2 = jb - 1; d2 <= je - 1; d2++)
        {
            int idx_d2 = idx_d3 + nx_block_ * d2;
            for(d1 = ib - 1; d1 <= ie - 1; d1++)
            {
                int idx = idx_d2 + d1;

                CX[idx] += *eta1 * CS0[idx];
                CR[idx] -= *eta1 * CQ[idx];
            }
        }
    }
}

//基本不算热点了
void c_pcg_btrop_check_(double* rr)
{

    local_sum = 0.0;
    memset(local_sum_red, 0, sizeof(realt) * VEC_RED_LEN);

    int iblock, d1, d2, d3;
    for(iblock = 0; iblock < nblocks_tropic_; iblock++)
    {
        int idx_d3 = ny_block_ * nx_block_ * iblock;

        for(d2 = jbs_[iblock] - 1; d2 <= jes_[iblock] - 1; d2++)
        {

            int idx_d2 = idx_d3 + nx_block_ * d2;

            int idx00 = idx_d2 + ib - 1;//idx i j
            //p -- plus 1, m -- minus 1
            int idx0p = idx00 + nx_block_;
            int idx0m = idx00 - nx_block_;
            int idxm0 = idx00 - 1;
            int idxmp = idxm0 + nx_block_;
            int idxmm = idxm0 - nx_block_;
            int idxp0 = idx00 + 1;
            int idxpp = idxp0 + nx_block_;
            int idxpm = idxp0 - nx_block_;

            #define CONV_UNROLL 16
            //这个数值一定要小于VEC_RED_LEN

            int unroll_cnt = i_inner_len / CONV_UNROLL;
            // unroll_cnt = 0;

            realt A00_red[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt AN0_00_red[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt AN0_0m_red[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt AE0_00_red[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt AE0_m0_red[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt ANE0_00_red[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt ANE0_0m_red[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt ANE0_m0_red[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt ANE0_mm_red[CONV_UNROLL] __attribute__ ((aligned(64)));

            realt CS_00[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt CS_0p[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt CS_0m[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt CS_p0[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt CS_m0[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt CS_pp[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt CS_pm[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt CS_mp[CONV_UNROLL] __attribute__ ((aligned(64)));
            realt CS_mm[CONV_UNROLL] __attribute__ ((aligned(64)));

            realt conv_q[CONV_UNROLL] __attribute__ ((aligned(64)));

            for(d1 = 0; d1 < unroll_cnt * CONV_UNROLL; d1 += CONV_UNROLL)
            {
                // int idx_plus = d1 * CONV_UNROLL;//大错特错
                int d0;

                //read
                for(d0 = 0; d0 < CONV_UNROLL; d0++)
                {
                    A00_red    [d0] = A00 [idx00 + d1 + d0];
                    AN0_00_red [d0] = AN0 [idx00 + d1 + d0];
                    AN0_0m_red [d0] = AN0 [idx0m + d1 + d0];
                    AE0_00_red [d0] = AE0 [idx00 + d1 + d0];
                    AE0_m0_red [d0] = AE0 [idxm0 + d1 + d0];
                    ANE0_00_red[d0] = ANE0[idx00 + d1 + d0];
                    ANE0_0m_red[d0] = ANE0[idx0m + d1 + d0];
                    ANE0_m0_red[d0] = ANE0[idxm0 + d1 + d0];
                    ANE0_mm_red[d0] = ANE0[idxmm + d1 + d0];

                    CS_00[d0]  = CX[idx00 + d1 + d0];
                    CS_0p[d0]  = CX[idx0p + d1 + d0];
                    CS_0m[d0]  = CX[idx0m + d1 + d0];
                    CS_p0[d0]  = CX[idxp0 + d1 + d0];
                    CS_m0[d0]  = CX[idxm0 + d1 + d0];
                    CS_pp[d0]  = CX[idxpp + d1 + d0];
                    CS_pm[d0]  = CX[idxpm + d1 + d0];
                    CS_mp[d0]  = CX[idxmp + d1 + d0];
                    CS_mm[d0]  = CX[idxmm + d1 + d0];
                }

                //conv
                for(d0 = 0; d0 < CONV_UNROLL; d0++)
                {
                    conv_q[d0] = A00_red    [d0] * CS_00[d0] +
                                 AN0_00_red [d0] * CS_0p[d0] +
                                 AN0_0m_red [d0] * CS_0m[d0] +
                                 AE0_00_red [d0] * CS_p0[d0] +
                                 AE0_m0_red [d0] * CS_m0[d0] +
                                 ANE0_00_red[d0] * CS_pp[d0] +
                                 ANE0_0m_red[d0] * CS_pm[d0] +
                                 ANE0_m0_red[d0] * CS_mp[d0] +
                                 ANE0_mm_red[d0] * CS_mm[d0]; 
                }

                //store and local_sum
                for(d0 = 0; d0 < CONV_UNROLL; d0++)
                {
                    int idx = idx_d2 + (d1 + ib - 1) + d0;
                    realt resu = B0[idx] - conv_q[d0];
                    local_sum_red[d0] += resu * resu;
                }
            }


            for(d1 = unroll_cnt * CONV_UNROLL; d1 < i_inner_len; d1++)
            {
                int idx = idx_d2 + (d1 + ib - 1);
                realt resu = A00 [idx00 + d1] * CX[idx00 + d1] +
                                 AN0 [idx00 + d1] * CX[idx0p + d1] +
                                 AN0 [idx0m + d1] * CX[idx0m + d1] +
                                 AE0 [idx00 + d1] * CX[idxp0 + d1] +
                                 AE0 [idxm0 + d1] * CX[idxm0 + d1] +
                                 ANE0[idx00 + d1] * CX[idxpp + d1] +
                                 ANE0[idx0m + d1] * CX[idxpm + d1] +
                                 ANE0[idxm0 + d1] * CX[idxmp + d1] +
                                 ANE0[idxmm + d1] * CX[idxmm + d1];  //AX, X全是参数传入

                resu= B0[idx] - resu;

                local_sum_red[0] += resu * resu;
            }
        }
    }

    for(d1 = 0; d1 < VEC_RED_LEN; d1++)
    {
        local_sum += local_sum_red[d1];
    }

    regular_masked_global_sum_dbl0_(rr);
}

//在check之前，从C代码本地将结果拷贝回去
//只需要R和X
void c_pcg_restore_result_(double* eta1)
{
    // return;
    int iblock, d1, d2, d3;

    for(iblock = 0; iblock < nblocks_tropic_; iblock++)
    {
        int idx_d3 = ny_block_ * nx_block_ * iblock;
        for(d2 = jb - 1; d2 <= je - 1; d2++)
        {
            int idx_d2 = idx_d3 + nx_block_ * d2;
            for(d1 = ib - 1; d1 <= ie - 1; d1++)
            {
                int idx = idx_d2 + d1;

                R_[idx] = CR[idx];
                X_[idx] = CX[idx];
            }
        }
    }
}

void c_pcg_restore_x_(double* eta1)
{
    // return;
    int iblock, d1, d2, d3;

    for(iblock = 0; iblock < nblocks_tropic_; iblock++)
    {
        int idx_d3 = ny_block_ * nx_block_ * iblock;
        for(d2 = jb - 1; d2 <= je - 1; d2++)
        {
            int idx_d2 = idx_d3 + nx_block_ * d2;
            for(d1 = ib - 1; d1 <= ie - 1; d1++)
            {
                int idx = idx_d2 + d1;

                // R_[idx] = CR[idx];
                X_[idx] = CX[idx];
            }
        }
    }
}