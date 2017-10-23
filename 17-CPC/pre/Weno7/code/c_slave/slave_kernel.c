#include "../c_header/slave_kernel.h"
#include "../c_header/c_public_const.h"

#include "slave.h"
#include "dma.h"
#include "simd.h"
#include <stdio.h>
#include <memory.h>

//param
__thread_local struct weno_init_param param;

//control
__thread_local int my_core;
__thread_local volatile long local_flag[FLAG_SIZE];
__thread_local volatile long out_flag[FLAG_SIZE];
__thread_local long slave_flag_to_wait;
__thread_local volatile intv8 bcast_tmp256;
__thread_local volatile unsigned long get_reply, put_reply;
__thread_local volatile unsigned long get_reply_god, put_reply_god;
__thread_local volatile unsigned long get_reply_target, put_reply_target;
__thread_local volatile long stcc, edcc;
__thread_local volatile long stcc2, edcc2;
__thread_local volatile long stcc3, edcc3;

//weno
__thread_local double* cpe_hj;
__thread_local double* cpe_u_in;
__thread_local double* cpe_u_in_next;
__thread_local double* cpe_u_out;
__thread_local int max_step20_u_out;
__thread_local int max_step40_u_out;

#define SWAP_DB_PTR(ptra, ptrb) {double* tmp; tmp = ptra; ptra = ptrb; ptrb = tmp;}

//align
void* _nice_ptr_head(char* ptr, int alignment, int right_shft)
{
    long ptrl = (long)ptr;
    int rem = ptrl % alignment;
    if(rem > 0)
    {
        return ptr + (alignment - rem + right_shft);
    }
    else
    {
        return ptr + right_shft;
    }
}

__thread_local double* _alloc_cpe_hj;
__thread_local double* _alloc_cpe_u_in_1;
__thread_local double* _alloc_cpe_u_in_2;
__thread_local double* _alloc_cpe_u_out;


//debug
__thread_local int single_print_flag;
__thread_local int cpe_cur_k, cpe_cur_j;

//optimize hj, vectorize
#define _MM_MA4_PD(out, a0, a1, a2, a3, b0, b1, b2, b3) \
out = a0 * b0;                                          \
out = a1 * b1 + out;                                    \
out = a2 * b2 + out;                                    \
out = a3 * b3 + out;

#define _MM_VM4_PD(out, l_vec_ptr, vec1, vec2, vec3, vec4)                                    \
{                                                                                             \
    doublev4 brod1 = *(l_vec_ptr    );                                                        \
    doublev4 brod2 = *(l_vec_ptr + 1);                                                        \
    doublev4 brod3 = *(l_vec_ptr + 2);                                                        \
    doublev4 brod4 = *(l_vec_ptr + 3);                                                        \
    _MM_MA4_PD(out, brod1, brod2, brod3, brod4, vec1, vec2, vec3, vec4);                      \
}

#define _MM_VM4_PD_2(out, out_2, l_vec_ptr, vec1, vec2, vec3, vec4, vec1_2, vec2_2, vec3_2, vec4_2) \
{                                                                                                   \
    doublev4 brod1 = *(l_vec_ptr    );                                                              \
    doublev4 brod2 = *(l_vec_ptr + 1);                                                              \
    doublev4 brod3 = *(l_vec_ptr + 2);                                                              \
    doublev4 brod4 = *(l_vec_ptr + 3);                                                              \
    _MM_MA4_PD(out, brod1, brod2, brod3, brod4, vec1, vec2, vec3, vec4);                            \
    _MM_MA4_PD(out_2, brod1, brod2, brod3, brod4, vec1_2, vec2_2, vec3_2, vec4_2);                  \
}

__thread_local const double S0_mat_row[16] __attribute__((aligned(32))) =
{0.2308680555555555181e1,-0.8176041666666666430e1,  0.9759375000000000355e1,-0.3892013888888888662e1,
        -0.8217708333333332504e1, 0.2973645833333333144e2, -0.3631979166666666003e2, 0.1480104166666666643e2,
        0.9926041666666666430e1,-0.3669479166666666003e2,  0.4661145833333333144e2,-0.1984270833333333073e2,
        -0.4017013888888889106e1, 0.1513437500000000036e2, -0.2005104166666666288e2, 0.8933680555555556069e1};

__thread_local const double S1_mat_row[16] __attribute__((aligned(32))) =
{0.1100347222222222143e1, -0.3384374999999999911e1,   0.3301041666666666874e1, -0.1017013888888888884e1,
        -0.3342708333333332948e1,  0.1161145833333333144e2,  -0.1219479166666666536e2,  0.3926041666666666874e1,
        0.3301041666666666874e1, -0.1231979166666666536e2,   0.1423645833333333144e2, -0.5217708333333332504e1,
        -0.1058680555555555403e1,  0.4092708333333333393e1,  -0.5342708333333332504e1,  0.2308680555555555181e1};



__thread_local const double S2_mat_row[16] __attribute__((aligned(32))) =
{0.1225347222222222143e1, -0.3176041666666666874e1,  0.3009374999999999911e1,  -0.1058680555555555403e1,
        -0.3051041666666666874e1,  0.1098645833333333321e2, -0.1231979166666666536e2,   0.4384374999999999467e1,
        0.2842708333333332948e1, -0.1219479166666666536e2,  0.1486145833333333144e2,  -0.5509374999999999467e1,
        -0.1017013888888888884e1,  0.4384374999999999467e1, -0.5551041666666665542e1,   0.2183680555555555625e1};



__thread_local const double S3_mat_row[16] __attribute__((aligned(32))) =
{0.8933680555555556069e1, -0.2005104166666666288e2,   0.1513437500000000036e2, -0.4017013888888889106e1,
        -0.1984270833333333073e2,  0.4661145833333333144e2,  -0.3669479166666666003e2,  0.9926041666666666430e1,
        0.1480104166666666643e2, -0.3631979166666666003e2,   0.2973645833333333144e2, -0.8217708333333332504e1,
        -0.3892013888888888662e1,  0.9759375000000000355e1,  -0.8176041666666666430e1,  0.2308680555555555181e1};

__thread_local double e_mat[16] __attribute__((aligned(32)));

__thread_local const double c[4] __attribute__((aligned(32))) = {1.e0/35.e0, 12.e0/35.e0, 18.e0/35.e0, 4.e0/35.e0};
__thread_local double res[4] __attribute__((aligned(32)));
__thread_local double S[4] __attribute__((aligned(32)));
__thread_local double q[16] __attribute__((aligned(32)));

__inline doublev4 __attribute__((__always_inline__)) batch4_vmv4x4_rowver(const double* __attribute__((aligned(32))) mat,
    doublev4 r_row1, doublev4 r_row2, doublev4 r_row3, doublev4 r_row4)
{


    doublev4 row_res[4];
    double value1, value2, value3, value4;
    int i;
    for(i=0; i<4; i++)
    {
        // value1 = *(mat + i * 4);
        // value2 = *(mat + i * 4 + 1);
        // value3 = *(mat + i * 4 + 2);
        // value4 = *(mat + i * 4 + 3);

        // doublev4 brod1 = simd_set_doublev4(value1, value1, value1, value1);
        // doublev4 brod2 = simd_set_doublev4(value2, value2, value2, value2);
        // doublev4 brod3 = simd_set_doublev4(value3, value3, value3, value3);
        // doublev4 brod4 = simd_set_doublev4(value4, value4, value4, value4);

        // doublev4 brod1 = *(mat + i * 4);
        // doublev4 brod2 = *(mat + i * 4 + 1);
        // doublev4 brod3 = *(mat + i * 4 + 2);
        // doublev4 brod4 = *(mat + i * 4 + 3);
        doublev4 data;
        simd_load(data, mat + i * 4);
        doublev4 brod1 = simd_vshff(data, data, 0x00);
        doublev4 brod2 = simd_vshff(data, data, 0x55);
        doublev4 brod3 = simd_vshff(data, data, 0xAA);
        doublev4 brod4 = simd_vshff(data, data, 0xFF);
        _MM_MA4_PD(row_res[i], brod1, brod2, brod3, brod4, r_row1, r_row2, r_row3, r_row4);
    }

    doublev4 res;
    _MM_MA4_PD(res, row_res[0], row_res[1], row_res[2], row_res[3], r_row1, r_row2, r_row3, r_row4);

    return res;
}


#define BATCH4_VMV4X8_ROWVER(res, res_2, mat,                                                                  \
    r_row1, r_row2, r_row3, r_row4, r_row5, r_row6, r_row7, r_row8)                                            \
{                                                                                                              \
    doublev4 row_res[4];                                                                                       \
    doublev4 row_res_2[4];                                                                                     \
    int i;                                                                                                     \
    for(i=0; i<4; i++)                                                                                         \
    {                                                                                                          \
        doublev4 brod1 = *(mat + i * 4);                                                                       \
        doublev4 brod2 = *(mat + i * 4 + 1);                                                                   \
        doublev4 brod3 = *(mat + i * 4 + 2);                                                                   \
        doublev4 brod4 = *(mat + i * 4 + 3);                                                                   \
        _MM_MA4_PD(row_res[i],   brod1, brod2, brod3, brod4, r_row1, r_row2, r_row3, r_row4);                  \
        _MM_MA4_PD(row_res_2[i], brod1, brod2, brod3, brod4, r_row5, r_row6, r_row7, r_row8);                  \
    }                                                                                                          \
    _MM_MA4_PD(res, row_res[0], row_res[1], row_res[2], row_res[3], r_row1, r_row2, r_row3, r_row4);           \
    _MM_MA4_PD(res_2, row_res_2[0], row_res_2[1], row_res_2[2], row_res_2[3], r_row5, r_row6, r_row7, r_row8); \
}

int BLOCK_SIZE(int block_id, int total_blocks, int n)
{
    return (n / total_blocks) + ((n % total_blocks > block_id) ? 1 : 0);
}

int BLOCK_LOW(int block_id, int total_blocks, int n)
{
    return (n / total_blocks) * block_id + ((n % total_blocks > block_id) ? block_id : n % total_blocks);
}

void cpe0_bcast_256(int my_core, volatile intv8* bcast_data)
{
    bcast_tmp256 = *bcast_data;

    if (my_core / 8 == 0)
    {
        __builtin_sw64_putC(bcast_tmp256, 0x0000000F);
    } else {
        bcast_tmp256 = __builtin_sw64_getxC(bcast_tmp256);
    }

    if (my_core % 8 == 0)
    {
        asm volatile ("nop":::"memory");
        __builtin_sw64_putR(bcast_tmp256, 0x0000000F);
        *bcast_data = bcast_tmp256;
    } else {
        asm volatile ("nop":::"memory");
        *bcast_data = __builtin_sw64_getxR(bcast_tmp256);
    }
}

void standard_wait_flag(int my_core)
{
    if(my_core == 0)
    {
        while(1)
        {
            get_reply = 0;
            athread_get(
                PE_MODE, param.host_flag, &local_flag[0],
                sizeof(long) * FLAG_SIZE, (void*)(&get_reply), 0, 0, 0
            );
            while(get_reply != 1);
            asm volatile ("#nop":::"memory");
            if(local_flag[0] >= slave_flag_to_wait)
                break;
        }

        slave_flag_to_wait++;
    }

    cpe0_bcast_256(my_core, &local_flag[0]);
    cpe0_bcast_256(my_core, &local_flag[4]);
    cpe0_bcast_256(my_core, &local_flag[8]);
    cpe0_bcast_256(my_core, &local_flag[12]);
}

void standard_write_flag(int my_core)
{
    out_flag[0] = out_flag[0] + 1;

    if(my_core == 0)
    {
        put_reply = 0;
        athread_put(
            PE_MODE, &out_flag[0], param.slave_flag,
            sizeof(long) * FLAG_SIZE, (void*)(&put_reply), 0, 0
        );
        while(put_reply != 1);

    }
    athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
}

void seq_write_data(void* mpe_data_ptr, void* cpe_data_ptr, int size, int my_core, int write_flag)
{
    int wcid;
    for (wcid = 63; wcid >= 0; wcid--)
    {
        if (my_core == wcid && write_flag != 0)
        {
            put_reply = 0;
            athread_put(
                PE_MODE, cpe_data_ptr, mpe_data_ptr,
                size, (void*)(&put_reply), 0, 0
            );
            while (put_reply != 1);
        }
        //athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
    }
}

void async_write_data(void* mpe_data_ptr, void* cpe_data_ptr, int size, int my_core, int write_flag)
{
    int wcid;
    //for (wcid = 63; wcid >= 0; wcid--)
    //{
    //    if (my_core == wcid && write_flag != 0)
    //    {
            //put_reply = 0;
            athread_put(
                PE_MODE, cpe_data_ptr, mpe_data_ptr,
                size, (void*)(&put_reply_god), 0, 0
            );
            put_reply_target++;
            //while (put_reply != 1);
    //    }
    //    //athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
    //}
}

inline static void wait_all_async_write_data()
{
    while(put_reply_god != put_reply_target);
}

inline static void wait_all_async_get_data()
{
    while(get_reply_god != get_reply_target);
}



void weno_hj_compute(int i_start, int i_end)
{
    int i;

    double e11,e12,e13,e14,e21,e22,e23,e24,e31,e32,e33,e34,e41,e42,e43,e44;
    double hx,hx_1,hx_6,hx_12,hx_60,hx_420,hx_840,hx_2520;

    double S0,S1,S2,S3,S10,S11,S12,S13,S20,S21,S22,S23,S30,S31,S32,S33,a0,a1,a2,a3,am,q0,q1,q2,q3,ep;
    double S0S0, Sl0_2;
    double S1S1, Sl1_2;
    double S2S2, Sl2_2;
    double S3S3, Sl3_2;
    double S0S3_2;
    double S1S2_2;
    double above_val, below_val;

    const double C0=1.e0/35.e0, C1=12.e0/35.e0, C2=18.e0/35.e0,  C3=4.e0/35.e0;

    const double a11=-2.e0/6.e0, a12=9.e0/6.e0, a13=-18.e0/6.e0, a14=11.e0/6.e0,
                 a21=1.e0/6.e0,                 a23=3.e0/6.e0,   a24=2.e0/6.e0 ,
                 a31=-2.e0/6.e0, a32=-3.e0/6.e0,                 a34=-1.e0/6.e0,
                 a41=-11.e0/6.e0,a42=18.e0/6.e0,a43=-9.e0/6.e0,  a44=2.e0/6.e0 ;

    const double b12=4.e0, b13=-5.e0, b14=2.e0,
                 b22= -2.e0,

                 b41=2.e0, b42=-5.e0, b43=4.e0,
                 c12=3.e0,

                 d12=13.e0/12.e0, d13=1043.e0/960.e0, d14=1.e0/12.e0;

    hx = param.hx;

    hx_1=1.e0/hx;
    hx_6=1.e0/(6.e0*hx);
    hx_12=1.e0/(12.e0*hx);
    hx_60=1.e0/(60.e0*hx);
    hx_420=1.e0/(420.e0*hx);
    hx_840=1.e0/(840.e0*hx);
    hx_2520=1.e0/(2520.e0*hx);

    e11=-3.e0/12.e0*hx_1;e12=13.e0/12.e0*hx_1;e13=-23.e0/12.e0*hx_1;e14=25.e0/12.e0*hx_1;
    e21=1.e0/12.e0*hx_1;e22=-5.e0/12.e0*hx_1;e23=13.e0/12.e0*hx_1;e24=3.e0/12.e0*hx_1;
    e31=-1.e0/12.e0*hx_1;e32=7.e0/12.e0*hx_1;e33=7.e0/12.e0*hx_1;e34=-1.e0/12.e0*hx_1;
    e41=3.e0/12.e0*hx_1;e42=13.e0/12.e0*hx_1;e43=-5.e0/12.e0*hx_1;e44=1.e0/12.e0*hx_1;

    e_mat[0]=-3.e0/12.e0*hx_1; e_mat[1]=13.e0/12.e0*hx_1;  e_mat[2]=-23.e0/12.e0*hx_1;  e_mat[3]=25.e0/12.e0*hx_1;
    e_mat[4]=1.e0/12.e0*hx_1;  e_mat[5]=-5.e0/12.e0*hx_1;  e_mat[6]=13.e0/12.e0*hx_1;   e_mat[7]=3.e0/12.e0*hx_1;
    e_mat[8]=-1.e0/12.e0*hx_1; e_mat[9]=7.e0/12.e0*hx_1;   e_mat[10]=7.e0/12.e0*hx_1;   e_mat[11]=-1.e0/12.e0*hx_1;
    e_mat[12]=3.e0/12.e0*hx_1; e_mat[13]=13.e0/12.e0*hx_1; e_mat[14]=-5.e0/12.e0*hx_1;  e_mat[15]=1.e0/12.e0*hx_1;

    ep=1.e-16;

    double fm3, fm2, fm1, f0, fp1, fp2, fp3;

    int total_len = i_end - i_start + 1;
    int hahablocks = total_len / 4;

    for(i = i_start; i <= i_end - 8; i += 8)
    {
        doublev4 r_row1  = simd_set_doublev4(cpe_u_in[i +  0], cpe_u_in[i +  1], cpe_u_in[i +  2], cpe_u_in[i +  3]);
        doublev4 r_row4  = simd_set_doublev4(cpe_u_in[i +  3], cpe_u_in[i +  4], cpe_u_in[i +  5], cpe_u_in[i +  6]);
        doublev4 r_row7  = simd_set_doublev4(cpe_u_in[i +  6], cpe_u_in[i +  7], cpe_u_in[i +  8], cpe_u_in[i +  9]);
        doublev4 r_row10 = simd_set_doublev4(cpe_u_in[i +  9], cpe_u_in[i + 10], cpe_u_in[i + 11], cpe_u_in[i + 12]);
        doublev4 r_row11 = simd_set_doublev4(cpe_u_in[i + 10], cpe_u_in[i + 11], cpe_u_in[i + 12], cpe_u_in[i + 13]);


        // doublev4 r_row2 = simd_vshff(r_row1, r_row4, 0x61);
        // doublev4 r_row3 = simd_vshff(r_row1, r_row4, 0xB6);
        // doublev4 r_row5 = simd_vshff(r_row4, r_row7, 0x61);
        // doublev4 r_row6 = simd_vshff(r_row4, r_row7, 0xB6);
        // doublev4 r_row8 = simd_vshff(r_row7, r_row10, 0x61);
        // doublev4 r_row9 = simd_vshff(r_row7, r_row10, 0xB6);

        doublev4 r_row2 = simd_vshff(r_row4,   r_row1,  0x49);
        doublev4 r_row3 = simd_vshff(r_row4,   r_row1,  0x9E);
        doublev4 r_row5 = simd_vshff(r_row7,   r_row4,  0x49);
        doublev4 r_row6 = simd_vshff(r_row7,   r_row4,  0x9E);
        doublev4 r_row8 = simd_vshff(r_row10,  r_row7,  0x49);
        doublev4 r_row9 = simd_vshff(r_row10,  r_row7,  0x9E);

        // doublev4 S0_batch4   = batch4_vmv4x4_rowver(S0_mat_row, r_row1, r_row2,  r_row3,  r_row4);
        // doublev4 S0_batch4_2 = batch4_vmv4x4_rowver(S0_mat_row, r_row5, r_row6,  r_row7,  r_row8);

        // doublev4 S1_batch4   = batch4_vmv4x4_rowver(S1_mat_row, r_row2, r_row3,  r_row4,  r_row5);
        // doublev4 S1_batch4_2 = batch4_vmv4x4_rowver(S1_mat_row, r_row6, r_row7,  r_row8,  r_row9);

        // doublev4 S2_batch4   = batch4_vmv4x4_rowver(S2_mat_row, r_row3, r_row4,  r_row5,  r_row6);
        // doublev4 S2_batch4_2 = batch4_vmv4x4_rowver(S2_mat_row, r_row7, r_row8,  r_row9, r_row10);

        // doublev4 S3_batch4   = batch4_vmv4x4_rowver(S3_mat_row, r_row4, r_row5,  r_row6,  r_row7);
        // doublev4 S3_batch4_2 = batch4_vmv4x4_rowver(S3_mat_row, r_row8, r_row9, r_row10, r_row11);
        doublev4 S0_batch4  ;
        doublev4 S0_batch4_2;
        doublev4 S1_batch4  ;
        doublev4 S1_batch4_2;
        doublev4 S2_batch4  ;
        doublev4 S2_batch4_2;
        doublev4 S3_batch4  ;
        doublev4 S3_batch4_2;

        BATCH4_VMV4X8_ROWVER(S0_batch4, S0_batch4_2, S0_mat_row, r_row1, r_row2,  r_row3,  r_row4, r_row5, r_row6,  r_row7,  r_row8);
        BATCH4_VMV4X8_ROWVER(S1_batch4, S1_batch4_2, S1_mat_row, r_row2, r_row3,  r_row4,  r_row5, r_row6, r_row7,  r_row8,  r_row9);
        BATCH4_VMV4X8_ROWVER(S2_batch4, S2_batch4_2, S2_mat_row, r_row3, r_row4,  r_row5,  r_row6, r_row7, r_row8,  r_row9, r_row10);
        BATCH4_VMV4X8_ROWVER(S3_batch4, S3_batch4_2, S3_mat_row, r_row4, r_row5,  r_row6,  r_row7, r_row8, r_row9, r_row10, r_row11);

        doublev4 ep_vec = simd_set_doublev4(ep, ep, ep, ep);
        // doublev4 ep_vec = ep;
        S0_batch4 = S0_batch4 + ep_vec;
        S1_batch4 = S1_batch4 + ep_vec;
        S2_batch4 = S2_batch4 + ep_vec;
        S3_batch4 = S3_batch4 + ep_vec;

        S0_batch4_2 = S0_batch4_2 + ep_vec;
        S1_batch4_2 = S1_batch4_2 + ep_vec;
        S2_batch4_2 = S2_batch4_2 + ep_vec;
        S3_batch4_2 = S3_batch4_2 + ep_vec;

        doublev4 a_sum, a_sum_2;
        S0_batch4 = S0_batch4 * S0_batch4;
        S1_batch4 = S1_batch4 * S1_batch4;
        S2_batch4 = S2_batch4 * S2_batch4;
        S3_batch4 = S3_batch4 * S3_batch4;

        S0_batch4_2 = S0_batch4_2 * S0_batch4_2;
        S1_batch4_2 = S1_batch4_2 * S1_batch4_2;
        S2_batch4_2 = S2_batch4_2 * S2_batch4_2;
        S3_batch4_2 = S3_batch4_2 * S3_batch4_2;

        doublev4 r_c;
        // r_c = simd_set_doublev4(c[0], c[0], c[0], c[0]);
        r_c = c[0];
        S0_batch4 = r_c / S0_batch4;
        S0_batch4_2 = r_c / S0_batch4_2;

        // r_c = simd_set_doublev4(c[1], c[1], c[1], c[1]);
        r_c = c[1];
        S1_batch4 = r_c / S1_batch4;
        S1_batch4_2 = r_c / S1_batch4_2;

        r_c = simd_set_doublev4(c[2], c[2], c[2], c[2]);
        // r_c = c[2];
        S2_batch4 = r_c / S2_batch4;
        S2_batch4_2 = r_c / S2_batch4_2;

        // r_c = simd_set_doublev4(c[3], c[3], c[3], c[3]);
        r_c = c[3];
        S3_batch4 = r_c / S3_batch4;
        S3_batch4_2 = r_c / S3_batch4_2;

        a_sum = S0_batch4 + S1_batch4 + S2_batch4 + S3_batch4;
        a_sum_2 = S0_batch4_2 + S1_batch4_2 + S2_batch4_2 + S3_batch4_2;

        doublev4 q0_batch4, q1_batch4, q2_batch4, q3_batch4;
        doublev4 q0_batch4_2, q1_batch4_2, q2_batch4_2, q3_batch4_2;

        // _MM_VM4_PD(q0_batch4, e_mat     , r_row1, r_row2, r_row3, r_row4);
        // _MM_VM4_PD(q0_batch4_2, e_mat     , r_row5, r_row6,  r_row7,  r_row8);
        _MM_VM4_PD_2(q0_batch4, q0_batch4_2, e_mat, r_row1, r_row2, r_row3, r_row4, r_row5, r_row6,  r_row7,  r_row8);

        // _MM_VM4_PD(q1_batch4, e_mat + 4 , r_row2, r_row3, r_row4, r_row5);
        // _MM_VM4_PD(q1_batch4_2, e_mat + 4 , r_row6, r_row7,  r_row8,  r_row9);
        _MM_VM4_PD_2(q1_batch4, q1_batch4_2, e_mat + 4, r_row2, r_row3, r_row4, r_row5, r_row6,  r_row7,  r_row8, r_row9);

        // _MM_VM4_PD(q2_batch4, e_mat + 8 , r_row3, r_row4, r_row5, r_row6);
        // _MM_VM4_PD(q2_batch4_2, e_mat + 8 , r_row7, r_row8,  r_row9, r_row10);
        _MM_VM4_PD_2(q2_batch4, q2_batch4_2, e_mat + 8, r_row3, r_row4, r_row5, r_row6,  r_row7,  r_row8, r_row9, r_row10);

        // _MM_VM4_PD(q3_batch4, e_mat + 12, r_row4, r_row5, r_row6, r_row7);
        // _MM_VM4_PD(q3_batch4_2, e_mat + 12, r_row8, r_row9, r_row10, r_row11);
        _MM_VM4_PD_2(q3_batch4, q3_batch4_2, e_mat + 12, r_row4, r_row5, r_row6,  r_row7,  r_row8, r_row9, r_row10, r_row11);

        doublev4 out_vec, out_vec_2;
        _MM_MA4_PD(out_vec, S0_batch4, S1_batch4, S2_batch4, S3_batch4, q0_batch4, q1_batch4, q2_batch4, q3_batch4);
        _MM_MA4_PD(out_vec_2, S0_batch4_2, S1_batch4_2, S2_batch4_2, S3_batch4_2, q0_batch4_2, q1_batch4_2, q2_batch4_2, q3_batch4_2);

        out_vec = out_vec / a_sum;
        out_vec_2 = out_vec_2 / a_sum_2;

        simd_storeu(out_vec, cpe_hj + i - 1);
        simd_storeu(out_vec_2, cpe_hj + i - 1 + 4);


    }

    for(;i <= i_end - 4; i += 4)
    {
        doublev4 r_row1 = simd_set_doublev4(cpe_u_in[i + 0], cpe_u_in[i + 1], cpe_u_in[i + 2], cpe_u_in[i + 3]);
        doublev4 r_row2 = simd_set_doublev4(cpe_u_in[i + 1], cpe_u_in[i + 2], cpe_u_in[i + 3], cpe_u_in[i + 4]);
        doublev4 r_row3 = simd_set_doublev4(cpe_u_in[i + 2], cpe_u_in[i + 3], cpe_u_in[i + 4], cpe_u_in[i + 5]);
        doublev4 r_row4 = simd_set_doublev4(cpe_u_in[i + 3], cpe_u_in[i + 4], cpe_u_in[i + 5], cpe_u_in[i + 6]);
        doublev4 r_row5 = simd_set_doublev4(cpe_u_in[i + 4], cpe_u_in[i + 5], cpe_u_in[i + 6], cpe_u_in[i + 7]);
        doublev4 r_row6 = simd_set_doublev4(cpe_u_in[i + 5], cpe_u_in[i + 6], cpe_u_in[i + 7], cpe_u_in[i + 8]);
        doublev4 r_row7 = simd_set_doublev4(cpe_u_in[i + 6], cpe_u_in[i + 7], cpe_u_in[i + 8], cpe_u_in[i + 9]);

        doublev4 S0_batch4 = batch4_vmv4x4_rowver(S0_mat_row, r_row1, r_row2, r_row3, r_row4);
        doublev4 S1_batch4 = batch4_vmv4x4_rowver(S1_mat_row, r_row2, r_row3, r_row4, r_row5);
        doublev4 S2_batch4 = batch4_vmv4x4_rowver(S2_mat_row, r_row3, r_row4, r_row5, r_row6);
        doublev4 S3_batch4 = batch4_vmv4x4_rowver(S3_mat_row, r_row4, r_row5, r_row6, r_row7);

        doublev4 ep_vec = simd_set_doublev4(ep, ep, ep, ep);
        S0_batch4 = S0_batch4 + ep_vec;
        S1_batch4 = S1_batch4 + ep_vec;
        S2_batch4 = S2_batch4 + ep_vec;
        S3_batch4 = S3_batch4 + ep_vec;
        doublev4 a_sum;
        S0_batch4 = S0_batch4 * S0_batch4;
        S1_batch4 = S1_batch4 * S1_batch4;
        S2_batch4 = S2_batch4 * S2_batch4;
        S3_batch4 = S3_batch4 * S3_batch4;

        doublev4 r_c;
        r_c = simd_set_doublev4(c[0], c[0], c[0], c[0]);
        S0_batch4 = r_c / S0_batch4;
        r_c = simd_set_doublev4(c[1], c[1], c[1], c[1]);
        S1_batch4 = r_c / S1_batch4;
        r_c = simd_set_doublev4(c[2], c[2], c[2], c[2]);
        S2_batch4 = r_c / S2_batch4;
        r_c = simd_set_doublev4(c[3], c[3], c[3], c[3]);
        S3_batch4 = r_c / S3_batch4;

        a_sum = S0_batch4 + S1_batch4 + S2_batch4 + S3_batch4;
        doublev4 q0_batch4, q1_batch4, q2_batch4, q3_batch4;

        _MM_VM4_PD(q0_batch4, e_mat     , r_row1, r_row2, r_row3, r_row4);
        _MM_VM4_PD(q1_batch4, e_mat + 4 , r_row2, r_row3, r_row4, r_row5);
        _MM_VM4_PD(q2_batch4, e_mat + 8 , r_row3, r_row4, r_row5, r_row6);
        _MM_VM4_PD(q3_batch4, e_mat + 12, r_row4, r_row5, r_row6, r_row7);
        doublev4 out_vec;
        _MM_MA4_PD(out_vec, S0_batch4, S1_batch4, S2_batch4, S3_batch4, q0_batch4, q1_batch4, q2_batch4, q3_batch4);

        out_vec = out_vec / a_sum;
        simd_storeu(out_vec, cpe_hj + i - 1);

    }

    for(; i <= i_end; ++i)
    {
        fm3 = cpe_u_in[i - 0];
        fm2 = cpe_u_in[i + 1];
        fm1 = cpe_u_in[i + 2];
        f0  = cpe_u_in[i + 3];
        fp1 = cpe_u_in[i + 4];
        fp2 = cpe_u_in[i + 5];
        fp3 = cpe_u_in[i + 6];

        S10=a11*fm3+a12*fm2+a13*fm1 +a14*f0;
        S11=a21*fm2 -   fm1+a23*f0+a24*fp1;
        S12=a31*fm1+a32*f0  +    fp1+a34*fp2;
        S13=a41*f0  +a42*fp1+a43*fp2+a44*fp3;
        S20=-fm3+b12*fm2+b13*fm1+b14*f0;
        S21=             fm1+b22*f0  +fp1;
        S22=             f0  +b22*fp1+fp2;
        S23=b41*f0+b42*fp1+b43*fp2-fp3;
        S30=-fm3+c12*(fm2-fm1) +f0;
        S31=-fm2+c12*(fm1-f0)   +fp1;
        S32=-fm1+c12*(f0-fp1)   +fp2;
        S33=-f0  +c12*(fp1-fp2) +fp3;

        S0=S10*S10+d12*S20*S20  +d13*S30*S30 +d14*S10*S30;
        S1=S11*S11+d12*S21*S21  +d13*S31*S31 +d14*S11*S31;
        S2=S12*S12+d12*S22*S22  +d13*S32*S32 +d14*S12*S32;
        S3=S13*S13+d12*S23*S23  +d13*S33*S33 +d14*S13*S33;

        S0 += ep;
        S1 += ep;
        S2 += ep;
        S3 += ep;
        S0S0 = S0 * S0;
        S1S1 = S1 * S1;
        S2S2 = S2 * S2;
        S3S3 = S3 * S3;
        S0S3_2 = S0S0 * S3S3;
        S1S2_2 = S1S1 * S2S2;
        Sl0_2 = S1S2_2 * S3S3;
        Sl1_2 = S0S3_2 * S2S2;
        Sl2_2 = S0S3_2 * S1S1;
        Sl3_2 = S1S2_2 * S0S0;

        below_val = (C0 * Sl0_2 + C1 * Sl1_2 + C2 * Sl2_2 + C3 * Sl3_2);

        q0=e11*fm3+e12*fm2+e13*fm1 +e14*f0;
        q1=e21*fm2+e22*fm1+e23*f0   +e24*fp1;
        q2=e31*fm1+e32*f0  +e33*fp1 +e34*fp2;
        q3=e41*f0  +e42*fp1+e43*fp2 +e44*fp3;
        q0 *= C0;
        q1 *= C1;
        q2 *= C2;
        q3 *= C3;

        above_val = (q0 * Sl0_2 + q1 * Sl1_2 + q2 * Sl2_2 + q3 * Sl3_2);

        cpe_hj[i - 1] = above_val / below_val;

        //a0=C0/((ep+S0)*(ep+S0));
        //a1=C1/((ep+S1)*(ep+S1));
        //a2=C2/((ep+S2)*(ep+S2));
        //a3=C3/((ep+S3)*(ep+S3));

        //am=a0+a1+a2+a3;
        //q0=e11*fm3+e12*fm2+e13*fm1 +e14*f0;
        //q1=e21*fm2+e22*fm1+e23*f0   +e24*fp1;
        //q2=e31*fm1+e32*f0  +e33*fp1 +e34*fp2;
        //q3=e41*f0  +e42*fp1+e43*fp2 +e44*fp3;
        //cpe_hj[i - 1]=(a0*q0+a1*q1+a2*q2+a3*q3)/am;
    }
}

void weno_u_out_compute()
{
    double dt = param.dt;
    int i, step;

    doublev4 dtv4 = dt;
    doublev4 hji01, hji02, hji03, hji04, hji05;
    doublev4 hji11, hji12, hji13, hji14, hji15;
    doublev4 hji01d, hji02d, hji03d, hji04d, hji05d;
    doublev4 hji11d, hji12d, hji13d, hji14d, hji15d;//d means double
    doublev4 uo1, uo2, uo3, uo4, uo5;
    doublev4 uo1d, uo2d, uo3d, uo4d, uo5d;

    double* tmphj0, *tmphj1;//, *tmpuo;
    doublev4* tmpuo4;
    doublev4* tmpui4;

    tmphj0 = cpe_hj + 1;
    tmphj1 = cpe_hj;//alighed
    //tmpuo = cpe_u_out;
    tmpuo4 = (doublev4*)cpe_u_out;
    tmpui4 = (doublev4*)(cpe_u_in + param.lap);

    for(step = 0; step < max_step40_u_out; step++)
    {

        // hji01 = simd_set_doublev4(tmphj0[   0], tmphj0[   1], tmphj0[   2], tmphj0[   3]);
        // hji02 = simd_set_doublev4(tmphj0[   4], tmphj0[   5], tmphj0[   6], tmphj0[   7]);
        // hji03 = simd_set_doublev4(tmphj0[   8], tmphj0[   9], tmphj0[  10], tmphj0[  11]);
        // hji04 = simd_set_doublev4(tmphj0[  12], tmphj0[  13], tmphj0[  14], tmphj0[  15]);
        // hji05 = simd_set_doublev4(tmphj0[  16], tmphj0[  17], tmphj0[  18], tmphj0[  19]);

        // hji01d = simd_set_doublev4(tmphj0[   0 + 20], tmphj0[   1 + 20], tmphj0[   2 + 20], tmphj0[   3 + 20]);
        // hji02d = simd_set_doublev4(tmphj0[   4 + 20], tmphj0[   5 + 20], tmphj0[   6 + 20], tmphj0[   7 + 20]);
        // hji03d = simd_set_doublev4(tmphj0[   8 + 20], tmphj0[   9 + 20], tmphj0[  10 + 20], tmphj0[  11 + 20]);
        // hji04d = simd_set_doublev4(tmphj0[  12 + 20], tmphj0[  13 + 20], tmphj0[  14 + 20], tmphj0[  15 + 20]);
        hji05d = simd_set_doublev4(tmphj0[  16 + 20], tmphj0[  17 + 20], tmphj0[  18 + 20], tmphj0[  19 + 20]);

        //unaligned exception
        //simd_loadu(hji01, tmphj1 +  1);
        //simd_loadu(hji02, tmphj1 +  5);
        //simd_loadu(hji03, tmphj1 +  9);
        //simd_loadu(hji04, tmphj1 + 13);
        //simd_loadu(hji05, tmphj1 + 17);

        simd_load(hji11, tmphj1 +  0);
        simd_load(hji12, tmphj1 +  4);
        simd_load(hji13, tmphj1 +  8);
        simd_load(hji14, tmphj1 + 12);
        simd_load(hji15, tmphj1 + 16);

        simd_load(hji11d, tmphj1 +  0 + 20);
        simd_load(hji12d, tmphj1 +  4 + 20);
        simd_load(hji13d, tmphj1 +  8 + 20);
        simd_load(hji14d, tmphj1 + 12 + 20);
        simd_load(hji15d, tmphj1 + 16 + 20);

//
        doublev4 temp1 = simd_vshff(hji12, hji11, 0x4E);
        doublev4 temp2 = simd_vshff(hji13, hji12, 0x4E);
        doublev4 temp3 = simd_vshff(hji14, hji13, 0x4E);
        doublev4 temp4 = simd_vshff(hji15, hji14, 0x4E);
        doublev4 temp5 = simd_vshff(hji11d, hji15, 0x4E);

        hji01 = simd_vshff(temp1, hji11, 0x99);
        hji02 = simd_vshff(temp2, hji12, 0x99);
        hji03 = simd_vshff(temp3, hji13, 0x99);
        hji04 = simd_vshff(temp4, hji14, 0x99);
        hji05 = simd_vshff(temp5, hji15, 0x99);

        doublev4 temp1d = simd_vshff(hji12d, hji11d, 0x4E);
        doublev4 temp2d = simd_vshff(hji13d, hji12d, 0x4E);
        doublev4 temp3d = simd_vshff(hji14d, hji13d, 0x4E);
        doublev4 temp4d = simd_vshff(hji15d, hji14d, 0x4E);

        hji01d = simd_vshff(temp1d, hji11d, 0x99);
        hji02d = simd_vshff(temp2d, hji12d, 0x99);
        hji03d = simd_vshff(temp3d, hji13d, 0x99);
        hji04d = simd_vshff(temp4d, hji14d, 0x99);
//

        uo1 = (hji01 - hji11);
        uo2 = (hji02 - hji12);
        uo3 = (hji03 - hji13);
        uo4 = (hji04 - hji14);
        uo5 = (hji05 - hji15);

        uo1d = (hji01d - hji11d);
        uo2d = (hji02d - hji12d);
        uo3d = (hji03d - hji13d);
        uo4d = (hji04d - hji14d);
        uo5d = (hji05d - hji15d);


        tmpuo4[0    ] = tmpui4[0    ] + uo1 * dtv4;
        tmpuo4[1    ] = tmpui4[1    ] + uo2 * dtv4;
        tmpuo4[2    ] = tmpui4[2    ] + uo3 * dtv4;
        tmpuo4[3    ] = tmpui4[3    ] + uo4 * dtv4;
        tmpuo4[4    ] = tmpui4[4    ] + uo5 * dtv4;
        tmpuo4[0 + 5] = tmpui4[0 + 5] + uo1d * dtv4;
        tmpuo4[1 + 5] = tmpui4[1 + 5] + uo2d * dtv4;
        tmpuo4[2 + 5] = tmpui4[2 + 5] + uo3d * dtv4;
        tmpuo4[3 + 5] = tmpui4[3 + 5] + uo4d * dtv4;
        tmpuo4[4 + 5] = tmpui4[4 + 5] + uo5d * dtv4;

        tmpuo4 += 10;
        tmpui4 += 10;
        tmphj1 += 40;
        tmphj0 += 40;
    }
}



/*
void weno_u_out_compute(int i_start, int i_end)
{
    double dt = param.dt;
    int i;

    for(i = i_start; i <= i_end; i++)
    {
        if(param.my_id == 0 && i == 40 && cpe_cur_j == 40 && cpe_cur_k == 40)
        {
            printf("uo ori = %.8f, hji = %.8f, hji1 = %.8f, dt = %.8f, ", cpe_u_out[i - 1], cpe_hj[i], cpe_hj[i-1], dt);
        }
        cpe_u_out[i - 1] += (cpe_hj[i]-cpe_hj[i-1]) * dt;

        if(param.my_id == 0 && i == 40 && cpe_cur_j == 40 && cpe_cur_k == 40)
        {
            printf("uo fin = %.8f\n", cpe_u_out[i - 1]);
        }
    }
}*/

void weno_slave_step()
{
    int ig, loop_length, jkid, in_stride, out_stride, my_start, my_len, jk_trueid, jk_trueid_out, j, k;
    int require_out, require_in;
    double* ptr_in, *ptr_out;
    int remain_length = local_flag[REMAIN_POINT];
    in_stride = local_flag[IN_STRIDE];
    out_stride = local_flag[OUT_STRIDE];
    require_in = local_flag[REQUIRE_IN];
    require_out = local_flag[REQUIRE_OUT];

    ptr_in = (double*)(local_flag[IN_PTR]);
    ptr_out = (double*)(local_flag[OUT_PTR]);

    loop_length = local_flag[GROUP_SIZE];

    my_start = BLOCK_LOW(my_core, 64, loop_length * 64 + remain_length);
    my_len = BLOCK_SIZE(my_core, 64, loop_length * 64 + remain_length);

    if(my_core == 0 && param.my_id == 0)
    {
        //printf("ftw = %ld, loop_length = %d, nx = %d, in_stride = %d, out_stride = %d\n"
        //       , slave_flag_to_wait, loop_length, in_stride, out_stride);
    }

    single_print_flag = 0;

    //fetch for the first iteration
    if(my_len > 0)//为了防止世界被破坏
    {
        jkid = my_start + 0;
        j = jkid % param.ny;
        k = jkid / param.ny;
        jk_trueid = (j + param.lap) + (param.ny + 2 * param.lap) * (k + param.lap);

        get_reply = 0;
        athread_get(
                    PE_MODE, ptr_in + in_stride * jk_trueid, cpe_u_in,
                    sizeof(double) * require_in, (void*)(&get_reply), 0, 0, 0
                    );
        while(get_reply != 1);
    }

    for(ig = 0; ig < loop_length + 1; ig++)
    {
        if(ig != 0)
        {
            standard_wait_flag(my_core);
        }

        //if(my_core == 0 && param.my_id == 0)
        //    printf("ig = %d, step 0.777\n", ig);

        if(my_len == ig)//happens only when ig == loop_length
        {
            //seq_write_data(ptr_out + (out_stride * jkid + param.lap/* + 1 + param.lap*/), cpe_u_out/* + (1 + param.lap)*/,
            //               sizeof(double) * param.nx, my_core, 0
            //               );
            //standard_write_flag(my_core);
            ////////RPCC(edcc);
            out_flag[PF_WRITE_DATA] += edcc - stcc;
            stcc = edcc;
            break;
        }

        //this will decide if there is "next loop"

        wait_all_async_get_data();
        ////////RPCC(edcc);
        out_flag[PF_MEMCPY_U] += edcc - stcc;
        stcc = edcc;


        if(ig + 1 < my_len)
        {
            jkid = my_start + ig + 1;
            j = jkid % param.ny;
            k = jkid / param.ny;
            jk_trueid = (j + param.lap) + (param.ny + 2 * param.lap) * (k + param.lap);

            //get_reply = 0;
            //athread_get(
            //            PE_MODE, ptr_in + in_stride * jk_trueid, cpe_u_in,
            //            sizeof(double) * in_stride, (void*)(&get_reply), 0, 0, 0
            //            );
            //while(get_reply != 1);


            athread_get(
                    PE_MODE, ptr_in + in_stride * jk_trueid, cpe_u_in_next,
                    sizeof(double) * require_in, (void*)(&get_reply_god), 0, 0, 0
                    );
            get_reply_target++;
            //wait_all_async_get_data();
        }

        jkid = my_start + ig;//dbg
        j = jkid % param.ny;
        k = jkid / param.ny;
        cpe_cur_k = k + 1;
        cpe_cur_j = j + 1;//dbg
        jk_trueid_out = (j + param.lap) + (param.ny + 2 * param.lap) * (k + param.lap);
        //jk_trueid = jk_trueid_out;

        //wait_all_async_get_data();//this means the last loop's cpe_u_in_next is ok, which is this loop's cpe_u_in
        //if(ig == 0)
        //{
        //    get_reply = 0;
        //    athread_get(
        //                PE_MODE, ptr_in + in_stride * jk_trueid, cpe_u_in,
        //                sizeof(double) * in_stride, (void*)(&get_reply), 0, 0, 0
        //                );
        //    while(get_reply != 1);
        //}

        ////////RPCC(edcc);
        out_flag[PF_GET_DATA] += edcc - stcc;
        stcc = edcc;

        //memcpy(cpe_u_out, cpe_u_in + param.lap, param.nx * sizeof(double));


        //if(my_core == 0 && param.my_id == 0)
        //    printf("ig = %d, step 1.233\n", ig);

        //weno_hj_compute(1 + param.lap, param.nx - param.lap + 1);
        weno_hj_compute(1, param.nx + 1);

        ////////RPCC(edcc);
        out_flag[PF_HJ_COMPUTE] += edcc - stcc;
        stcc = edcc;

        //if(my_core == 0 && param.my_id == 0)
        //    printf("ig = %d, step 2.456\n", ig);

        //weno_u_out_compute(1 + param.lap, param.nx - param.lap);

        wait_all_async_write_data();//in case of overwrite data!

        ////////RPCC(edcc);
        out_flag[PF_MEMCPY_U] += edcc - stcc;//this is wait flag
        stcc = edcc;

        weno_u_out_compute(1, param.nx);

        ////////RPCC(edcc);
        out_flag[PF_U_OUT_COMPUTE] += edcc - stcc;
        stcc = edcc;

        //if(my_core == 0 && param.my_id == 0)
        //    printf("ig = %d, step 2.789\n", ig);

        //seq_write_data(ptr_out + (out_stride * jk_trueid + param.lap/* + 1 + param.lap*/), cpe_u_out/* + (1 + param.lap)*/,
        //               sizeof(double) * param.nx, my_core, 1
        //               );

        async_write_data(ptr_out + (out_stride * jk_trueid_out/* + 1 + param.lap*/), cpe_u_out/* + (1 + param.lap)*/,
                         sizeof(double) * require_out, my_core, 1
                         );

        //if(my_core == 0 && param.my_id == 0)
        //    printf("ig = %d, step 3.789\n", ig);

        //standard_write_flag(my_core);

        ////////RPCC(edcc);
        out_flag[PF_WRITE_DATA] += edcc - stcc;
        stcc = edcc;

        SWAP_DB_PTR(cpe_u_in_next, cpe_u_in);//cpe_u_in of this loop has been squeezed out yet.

        //if(my_core == 0 && param.my_id == 0)
        //    printf("ig = %d, step 4.789\n", ig);
    }

    standard_write_flag(my_core);

    ////////RPCC(edcc);
    out_flag[PF_WRITE_DATA] += edcc - stcc;
    stcc = edcc;
}

void weno_slave_init()
{
    max_step20_u_out = param.nx / 20;
    if(param.nx % 20 > 0)
        max_step20_u_out++;
    max_step40_u_out = param.nx / 40;
    if(param.nx % 40 > 0)
        max_step40_u_out++;

    _alloc_cpe_hj = (double*)ldm_malloc(sizeof(double) * (param.nx + SLAVE_SAFE_PAD));
    _alloc_cpe_u_in_2 = (double*)ldm_malloc(sizeof(double) * (param.nx + 2 * param.lap + SLAVE_SAFE_PAD));
    _alloc_cpe_u_in_1 = (double*)ldm_malloc(sizeof(double) * (param.nx + 2 * param.lap + SLAVE_SAFE_PAD));
    _alloc_cpe_u_out = (double*)ldm_malloc(sizeof(double) * (param.nx + SLAVE_SAFE_PAD));

    cpe_hj = _nice_ptr_head((char*)_alloc_cpe_hj, 32, 0);
    cpe_u_in = _nice_ptr_head((char*)_alloc_cpe_u_in_1, 32, param.lap % 4);
    cpe_u_in_next = _nice_ptr_head((char*)_alloc_cpe_u_in_2, 32, param.lap % 4);
    cpe_u_out = _nice_ptr_head((char*)_alloc_cpe_u_out, 32, 0);

    if(my_core == 0 && param.my_id == 0)
    {
        printf("hj u u = %x %x %x %x\n", cpe_hj, cpe_u_in, cpe_u_out, cpe_u_in_next);
    }
}

void weno_slave_finalize()
{
    ldm_free(_alloc_cpe_hj, sizeof(double) * (param.nx + SLAVE_SAFE_PAD));
    ldm_free(_alloc_cpe_u_in_1, sizeof(double) * (param.nx + 2 * 4 + SLAVE_SAFE_PAD));
    ldm_free(_alloc_cpe_u_in_2, sizeof(double) * (param.nx + 2 * 4 + SLAVE_SAFE_PAD));
    ldm_free(_alloc_cpe_u_out, sizeof(double) * (param.nx + SLAVE_SAFE_PAD));
}

void cpe_athread_daemon(void *_param)
{
    int i, j, k;
    stcc = 0;
    edcc = 0;

    param = *((struct weno_init_param*) _param);

    my_core = athread_get_id(-1);

    //must do this
    if(my_core == 0 && param.my_id == 0)
        printf("nx = %d, ny = %d, nz = %d, id = %d, host_flag = %x, slave_flag = %x\n",
           param.nx, param.ny, param.nz, param.my_id, param.host_flag, param.slave_flag);

    for (i = 0; i < FLAG_SIZE; i++){ local_flag[i] = 0; out_flag[i] = 0;}
    slave_flag_to_wait = 1;

    weno_slave_init();

    get_reply_god = 0, put_reply_god = 0;//the Jesus do not like waiting for someone
    get_reply_target = 0, put_reply_target = 0;

    ////////RPCC(stcc);

    while(1)
    {
        standard_wait_flag(my_core);
        ////////RPCC(edcc);
        out_flag[PF_OUT_WAIT_FLAG] += edcc - stcc;
        stcc = edcc;

        if(local_flag[KERNEL_ACTION] == WENO_INNER_FLAG)
        {
            weno_slave_step();
        }

        if (local_flag[KERNEL_ACTION] == EXIT_FLAG)
        {
            standard_write_flag(my_core);
            standard_write_flag(my_core);
            break;
        }
    }

    weno_slave_finalize();
}
