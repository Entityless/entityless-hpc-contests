#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "slave.h"
#include "dma.h"
#include "simd.h"

#include "../c_header/c_public_const.h"
#include "../c_header/slave_kernel.h"

//#define PPG_INNER_DATA (PPG_BI_DATA * _kl * _jl)

// ??????
__thread_local struct cpe_init_param param;                 // 28 * (4B / 8B) + 4B
__thread_local volatile long local_flag[FLAG_SIZE];          // 16 * 8B
__thread_local volatile long out_flag[FLAG_SIZE];          // 16 * 8B
__thread_local long slave_flag_to_wait;
__thread_local volatile unsigned long get_reply, put_reply; // 2 * 8B
__thread_local volatile unsigned long get_reply_god, put_reply_god; // 2 * 8B
__thread_local volatile unsigned long get_reply_target, put_reply_target; // 2 * 8B
__thread_local volatile intv8 bcast_tmp256;                 // 8B

// ??? c_public_var.c ???Χ??????
__thread_local int   cpe_xy_ranges[4];                      // 4 * 4B
__thread_local float cpe_const_values[8];                   // 8 * 4B
__thread_local int   cpe_ixs, cpe_ixl, cpe_iys, cpe_iyl;    // 4 * 4B
__thread_local float cpe_deltt5, cpe_deltt, cpe_cong;       // 3 * 4B
__thread_local float cpe_al31, cpe_al21, cpe_al13;          // 3 * 4B
__thread_local float cpe_al23, cpe_al11;                    // 2 * 4B
__thread_local int   cpe_ix_size, cpe_dim4, cpe_dim3;       // 6 * 4B
__thread_local int   cpe_dim2, cpe_iy_size, cpe_ixiy_size;  // 3 * 4B

/* ---------- for propagat ---------- */
__thread_local int cpe_ppg_cnt = 0;                         // 4B
//__thread_local float cpe_ppg_packed_data[PPG_INNER_DATA];   // 12 * 25 * 22 * 4B
__thread_local float cpe_ppg_e[_jl * _kl];                  // 12 * 25 * 4B
//Huang Chenghuan 20170320
__thread_local float cpe_ppg_ee_now[_jl * _kl * 10];//10 for mean1
__thread_local int cpe_ppg_idxs[_jl * _kl * 5];
__thread_local float cpe_ppg_values[_jl * _kl * 6];

__thread_local int ppg_src_x[6];
__thread_local int ppg_src_y[6];
__thread_local int ppg_dst_x[6];
__thread_local int ppg_dst_y[6];
__thread_local int ppg_input_mask[9];
__thread_local int ppg_copy_count;


/* ---------- for implsch ---------- */
__thread_local int cpe_ips_cnt = 0;                         // 4B
// ???????
__thread_local int   cpe_ips_ikp   [_kl];                   // 25 * 4B
__thread_local int   cpe_ips_ikp1  [_kl];                   // 25 * 4B
__thread_local int   cpe_ips_ikm   [_kl];                   // 25 * 4B
__thread_local int   cpe_ips_ikm1  [_kl];                   // 25 * 4B
__thread_local float cpe_ips_wks17 [_kl];                   // 25 * 4B
__thread_local float cpe_ips_grolim[_kl];                   // 25 * 4B
__thread_local float cpe_ips_wp [2 * 2 * _kl];              // 2 * 2 * 25 * 4B
__thread_local float cpe_ips_wm [2 * 2 * _kl];              // 2 * 2 * 25 * 4B
__thread_local float cpe_ips_wkh[_kldp1];                   // 31 * 4B
__thread_local float cpe_ips_wk [_kldp1];                   // 31 * 4B
__thread_local float cpe_ips_sqrt_wk [_kldp1];              // 31 * 4B
__thread_local float cpe_ips_sqrt_div_wk [_kldp1];          // 31 * 4B
__thread_local float cpe_ips_dwk[_kldp1];                   // 31 * 4B
__thread_local int   cpe_ips_jp1[2 * _jl];                  // 2 * 12 * 4B
__thread_local int   cpe_ips_jp2[2 * _jl];                  // 2 * 12 * 4B
__thread_local int   cpe_ips_jm1[2 * _jl];                  // 2 * 12 * 4B
__thread_local int   cpe_ips_jm2[2 * _jl];                  // 2 * 12 * 4B
__thread_local int   cpe_ips_jpms[8 * _jl];
__thread_local float cpe_ips_cos_theta[_jlp1];              // 13 * 4B
__thread_local float cpe_ips_sin_theta[_jlp1];              // 13 * 4B
// ??ε???????????????
__thread_local volatile int   cpe_ips_ia, cpe_ips_ic;                // 2 * 4B
__thread_local float cpe_ips_wf [_kldp1];                   // 31 * 4B
__thread_local float cpe_ips_ssbos [_kldp1];                // 31 * 4B
__thread_local float cpe_ips_ccg[_kldp1];                   // 31 * 4B
__thread_local float cpe_ips_packed_data[PPG_SV_DATA];      // 8 * 4B
__thread_local float cpe_ips_e[_kl * _jl];                  // 12 * 25 * 4B
// ???д???????
__thread_local float cpe_ips_ee[_kl * _jl];                 // 12 * 25 * 4B
__thread_local float cpe_pexx_output_buf[NUM_PEXX];         // 3 * 4B
// ?????????????д??????м????
struct c_implsch_local_tmp
{
    float asi, awk, ae, ark, awf, ks0, enh;                 // 6 * 4B
    float fconst0[_kl];                                     // 25 * 4B
    float   se[_klp1 * _jl + 8];                                // 26 * 12 * 4B
    float  dse[_klp1 * _jl + 8];                                // 26 * 12 * 4B
    float sein[_kl   * _jl + 8];                                // 25 * 12 * 4B
    float sebo[_kl   * _jl + 8];                                // 25 * 12 * 4B
    float seds[_kl   * _jl + 8];                                // 25 * 12 * 4B
};

__thread_local struct c_implsch_local_tmp thread_local_buff;

//__thread_local unsigned int cpe_jse_idxs[_kl * _jl * 2 * 8];
__thread_local unsigned int cpe_jse_cnts[_klp1];

//mean1 getem
__thread_local float cpe_mean1_dwf [_kldp1];
__thread_local float cpe_mean1_wf [_kldp1];
__thread_local float *cpe_mean1_next_dwf;
__thread_local float *cpe_mean1_next_wf ;
//__thread_local float cpe_mean1_ees[_jl * _kl * 5];
__thread_local float* cpe_mean1_ees;
__thread_local float* cpe_mean1_ees_next;
__thread_local float cpe_mean1_em[_jl * _kl]  __attribute__ ((aligned (32)));
__thread_local float cpe_mean1_data_in[GETEM_IN_FLOAT];
__thread_local float cpe_mean1_data_next[GETEM_IN_FLOAT];
__thread_local int cpe_mean1_index_in[GETEM_IN_INT];
__thread_local int cpe_mean1_index_next[GETEM_IN_INT];
__thread_local float cpe_mean1_data_out[GETEM_OUT_DATA];
__thread_local int   cpe_mean1_ia, cpe_mean1_ic, cpe_mean1_nsp;
__thread_local int   cpe_mean1_next_ia, cpe_mean1_next_ic, cpe_mean1_next_nsp;
__thread_local int cpe_ppg_ia, cpe_ppg_ic, ppgid, cpe_ppg_nsp;             // 2 * 4B
__thread_local int wild_my_core;

__thread_local char strlocal[100];
__thread_local volatile long stcc, edcc;
__thread_local volatile long stcc2, edcc2;
__thread_local volatile long stcc3, edcc3;

__thread_local float cpe_nomarine_es[_jl * _kl];

__thread_local const floatv4 onev4={1.0,1.0,1.0,1.0};
__thread_local float cpe_tmp_fv[8] __attribute__ ((aligned (32)));
__thread_local float cpe_tmp_fv1[8] __attribute__ ((aligned (32)));
__thread_local float cpe_tmp_fv2[8] __attribute__ ((aligned (32)));

__thread_local int cpe_ppg_kv[_klmjl] __attribute__ ((aligned (32)));
__thread_local int cpe_ppg_jv[_klmjl] __attribute__ ((aligned (32)));

//mean1 const

int BLOCK_SIZE(int block_id, int total_blocks, int n)
{
    return (n / total_blocks) + ((n % total_blocks > block_id) ? 1 : 0);
}

int BLOCK_LOW(int block_id, int total_blocks, int n)
{
    return (n / total_blocks) * block_id + ((n % total_blocks > block_id) ? block_id : n % total_blocks);
}

// ?? CPE0 ???????? CPE
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

// CPE 0 ???????????????? flag??????????????????????????????????????
void standard_wait_flag(int my_core)
{
    if(my_core == 0)
    {
        while(1)
        {
            get_reply = 0;
            athread_get(
                PE_MODE, (void*)param.host_flag, (void*)&local_flag[0],
                sizeof(long) * FLAG_SIZE, (void*)(&get_reply), 0, 0, 0
            );
            while(get_reply != 1);
            asm volatile ("#nop":::"memory");
            if(local_flag[0] >= slave_flag_to_wait)
                break;
        }
        slave_flag_to_wait++;
    }

    int i;
    for(i = 0; i < FLAG_SIZE / 4; i++)
    {
        cpe0_bcast_256(my_core, (void*)&local_flag[i * 4]);
    }
}

void standard_write_flag(int my_core)
{
    out_flag[0] = out_flag[0] + 1;
    athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
    if(my_core == 0)
    {
        put_reply = 0;
        athread_put(
            PE_MODE, (void*)&out_flag[0], (void*)param.slave_flag,
            sizeof(long) * FLAG_SIZE, (void*)(&put_reply), 0, 0
        );
        while(put_reply != 1);

    }
    athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
}

static __inline floatv4 __attribute__((__always_inline__)) bilinear_interpolation_qr(floatv4 aa, floatv4 bb, floatv4 cc, floatv4 dd, floatv4 q, floatv4 r)
{
    floatv4 phia = (onev4 - q) * (onev4 - r);
    floatv4 phib = q * (onev4 - r);
    floatv4 phic = (onev4 - q) * r;
    floatv4 phid = q * r;
    return (aa * phia + bb * phib + cc * phic + dd * phid);
}

void cpe_propagat_kernel(int cpe_ppg_ia, int cpe_ppg_ic, int my_core)
{
    int jk = 0;
    //float *packed_data_ptr;
    floatv4 e1s[4], e2s[4], e3s[4], e4s[4], es[4], exxyy;
    floatv4 q_4, r_4, q, r, fien, fien1;
    floatv4 temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8,
            temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16;
	 
    int ixx[4], jyy[4], iwk[4], iwk1[4], jth[4], jth1[4];
    int *idxs_ptr;
    float* values_ptr;
    int left_jk[4][4], right_iaic[4][4];
    int i, index[64];
    values_ptr = &cpe_ppg_values[0];
    idxs_ptr = &cpe_ppg_idxs[0];



    //packed_data_ptr = cpe_ppg_packed_data;
    for (jk = 0; jk < _jl * _kl; jk+=4)
    {
        q_4  = simd_set_floatv4(values_ptr[VALUES_Q_4_OFF]
                               ,values_ptr[VALUES_Q_4_OFF+6]
                               ,values_ptr[VALUES_Q_4_OFF+12]
                               ,values_ptr[VALUES_Q_4_OFF+18]);
        r_4  = simd_set_floatv4(values_ptr[VALUES_R_4_OFF]
                                ,values_ptr[VALUES_R_4_OFF+6]
                                ,values_ptr[VALUES_R_4_OFF+12]
                                ,values_ptr[VALUES_R_4_OFF+18]);
        q    = simd_set_floatv4(values_ptr[VALUES_Q_OFF]
                                ,values_ptr[VALUES_Q_OFF+6]
                                ,values_ptr[VALUES_Q_OFF+12]
                                ,values_ptr[VALUES_Q_OFF+18]);
        r    = simd_set_floatv4(values_ptr[VALUES_R_OFF]
                                ,values_ptr[VALUES_R_OFF+6]
                                ,values_ptr[VALUES_R_OFF+12]
                                ,values_ptr[VALUES_R_OFF+18]); 
        fien  = simd_set_floatv4(values_ptr[VALUES_FIEN_OFF]
                                ,values_ptr[VALUES_FIEN_OFF+6]
                                ,values_ptr[VALUES_FIEN_OFF+12]
                                ,values_ptr[VALUES_FIEN_OFF+18]); 
        fien1 = simd_set_floatv4(values_ptr[VALUES_FIEN1_OFF]
                                ,values_ptr[VALUES_FIEN1_OFF+6]
                                ,values_ptr[VALUES_FIEN1_OFF+12]
                                ,values_ptr[VALUES_FIEN1_OFF+18]); 
        
        int ixx[4]  =  { idxs_ptr[IDXS_IXX_OFF] - cpe_ppg_ia + 1       
                         ,idxs_ptr[IDXS_IXX_OFF+5] - cpe_ppg_ia + 1
                         ,idxs_ptr[IDXS_IXX_OFF+10] - cpe_ppg_ia + 1
                         ,idxs_ptr[IDXS_IXX_OFF+15] - cpe_ppg_ia + 1};                             
        int jyy[4]  =  { idxs_ptr[IDXS_JYY_OFF] - cpe_ppg_ic + 1
                         ,idxs_ptr[IDXS_JYY_OFF+5] - cpe_ppg_ic + 1
                         ,idxs_ptr[IDXS_JYY_OFF+10] - cpe_ppg_ic + 1
                         ,idxs_ptr[IDXS_JYY_OFF+15] - cpe_ppg_ic + 1};
        int jth[4]  =  { idxs_ptr[IDXS_JTH_OFF]
                         ,idxs_ptr[IDXS_JTH_OFF+5]
                         ,idxs_ptr[IDXS_JTH_OFF+10]
                         ,idxs_ptr[IDXS_JTH_OFF+15]};
        int iwk[4]  =  { idxs_ptr[IDXS_IWK_OFF]
                         ,idxs_ptr[IDXS_IWK_OFF+5]
                         ,idxs_ptr[IDXS_IWK_OFF+10]
                         ,idxs_ptr[IDXS_IWK_OFF+15]};
        int iwk1[4]  = { idxs_ptr[IDXS_IWK1_OFF]
                         ,idxs_ptr[IDXS_IWK1_OFF+5]
                         ,idxs_ptr[IDXS_IWK1_OFF+10]
                         ,idxs_ptr[IDXS_IWK1_OFF+15]};
                               


        for (i=0;i<4;i++){
            jth1[i] = jth[i] + 1;
            
            if (jth[i] == _jlp1){
                jth[i] = _jl;
                jth1[i] = 1;
            }
            if (jth1[i] == _jlp1) jth1[i] = 1;
            left_jk[0][i]    = (iwk[i]  - 1) + (jth[i]  - 1) * _kl;
            left_jk[1][i]    = (iwk1[i] - 1) + (jth[i]  - 1) * _kl;
            left_jk[2][i]    = (iwk[i]  - 1) + (jth1[i] - 1) * _kl;
            left_jk[3][i]   = (iwk1[i] - 1) + (jth1[i] - 1) * _kl;

            right_iaic[0][i] = (ixx[i]  + 3 * jyy[i] ) * _kl * _jl;
            right_iaic[1][i] = (ixx[i] + 1 + 3 * jyy[i] ) * _kl * _jl;
            right_iaic[2][i] = (ixx[i]  + 3 * jyy[i] + 3) * _kl * _jl;
            right_iaic[3][i] = (ixx[i] + 1 + 3 * jyy[i] +3) * _kl * _jl;
        }

        index[0] = left_jk[0][0] + right_iaic[0][0];
        index[1] = left_jk[0][1] + right_iaic[0][1];
        index[2] = left_jk[0][2] + right_iaic[0][2];
        index[3] = left_jk[0][3] + right_iaic[0][3];
        temp1 = simd_set_floatv4(cpe_ppg_ee_now[index[0]],
                                 cpe_ppg_ee_now[index[1]],
                                 cpe_ppg_ee_now[index[2]],
                                 cpe_ppg_ee_now[index[3]]); 
        e1s[0] = temp1 * fien;
        
        index[4] =  left_jk[1][0] + right_iaic[0][0];
        index[5] =  left_jk[1][1] + right_iaic[0][1];
        index[6] =  left_jk[1][2] + right_iaic[0][2];
        index[7] =  left_jk[1][3] + right_iaic[0][3];
        temp2 = simd_set_floatv4(cpe_ppg_ee_now[index[4]],
                                 cpe_ppg_ee_now[index[5]],
                                 cpe_ppg_ee_now[index[6]],
                                 cpe_ppg_ee_now[index[7]]); 
        e2s[0] = temp2 * fien1;
        
        index[8] =  left_jk[2][0] + right_iaic[0][0];
        index[9] =  left_jk[2][1] + right_iaic[0][1];
        index[10] =  left_jk[2][2] + right_iaic[0][2];
        index[11] =  left_jk[2][3] + right_iaic[0][3];
        temp3 = simd_set_floatv4(cpe_ppg_ee_now[index[8]],
                                 cpe_ppg_ee_now[index[9]],
                                 cpe_ppg_ee_now[index[10]],
                                 cpe_ppg_ee_now[index[11]]); 
        e3s[0] = temp3 * fien;
        
        index[12] =  left_jk[3][0] + right_iaic[0][0];
        index[13] =  left_jk[3][1] + right_iaic[0][1];
        index[14] =  left_jk[3][2] + right_iaic[0][2];
        index[15] =  left_jk[3][3] + right_iaic[0][3];
        temp4 = simd_set_floatv4(cpe_ppg_ee_now[index[12]],
                                 cpe_ppg_ee_now[index[13]],
                                 cpe_ppg_ee_now[index[14]],
                                 cpe_ppg_ee_now[index[15]]); 
        e4s[0] = temp4 * fien1;
        
        index[16] =  left_jk[0][0] + right_iaic[1][0];
        index[17] =  left_jk[0][1] + right_iaic[1][1];
        index[18] =  left_jk[0][2] + right_iaic[1][2];
        index[19] =  left_jk[0][3] + right_iaic[1][3];
        temp5 = simd_set_floatv4(cpe_ppg_ee_now[index[16]],
                                 cpe_ppg_ee_now[index[17]],
                                 cpe_ppg_ee_now[index[18]],
                                 cpe_ppg_ee_now[index[19]]); 
        e1s[1] = temp5 * fien;
        
        index[20] =  left_jk[1][0] + right_iaic[1][0];
        index[21] =  left_jk[1][1] + right_iaic[1][1];
        index[22] =  left_jk[1][2] + right_iaic[1][2];
        index[23] =  left_jk[1][3] + right_iaic[1][3];
        temp6 = simd_set_floatv4(cpe_ppg_ee_now[index[20]],
                                 cpe_ppg_ee_now[index[21]],
                                 cpe_ppg_ee_now[index[22]],
                                 cpe_ppg_ee_now[index[23]]); 
        e2s[1] = temp6 * fien1;
        
        
        index[24] =  left_jk[2][0] + right_iaic[1][0];
        index[25] =  left_jk[2][1] + right_iaic[1][1];
        index[26] =  left_jk[2][2] + right_iaic[1][2];
        index[27] =  left_jk[2][3] + right_iaic[1][3];
        temp7 = simd_set_floatv4(cpe_ppg_ee_now[index[24]],
                                 cpe_ppg_ee_now[index[25]],
                                 cpe_ppg_ee_now[index[26]],
                                 cpe_ppg_ee_now[index[27]]); 
        e3s[1] = temp7 * fien;
        
        
        index[28] =  left_jk[3][0] + right_iaic[1][0];
        index[29] =  left_jk[3][1] + right_iaic[1][1];
        index[30] =  left_jk[3][2] + right_iaic[1][2];
        index[31] =  left_jk[3][3] + right_iaic[1][3];
        temp8 = simd_set_floatv4(cpe_ppg_ee_now[index[28]],
                                 cpe_ppg_ee_now[index[29]],
                                 cpe_ppg_ee_now[index[30]],
                                 cpe_ppg_ee_now[index[31]]); 
        e4s[1] = temp8 * fien1;
        
        
        index[32] =  left_jk[0][0] + right_iaic[2][0];
        index[33] =  left_jk[0][1] + right_iaic[2][1];
        index[34] =  left_jk[0][2] + right_iaic[2][2];
        index[35] =  left_jk[0][3] + right_iaic[2][3];
        temp9 = simd_set_floatv4(cpe_ppg_ee_now[index[32]],
                                 cpe_ppg_ee_now[index[33]],
                                 cpe_ppg_ee_now[index[34]],
                                 cpe_ppg_ee_now[index[35]]); 
        e1s[2] = temp9 * fien;
        
        
        index[36] =  left_jk[1][0] + right_iaic[2][0];
        index[37] =  left_jk[1][1] + right_iaic[2][1];
        index[38] =  left_jk[1][2] + right_iaic[2][2];
        index[39] =  left_jk[1][3] + right_iaic[2][3];
        temp10 = simd_set_floatv4(cpe_ppg_ee_now[index[36]],
                                  cpe_ppg_ee_now[index[37]],
                                  cpe_ppg_ee_now[index[38]],
                                  cpe_ppg_ee_now[index[39]]); 
        e2s[2] = temp10 * fien1;
        
        index[40] =  left_jk[2][0] + right_iaic[2][0];
        index[41] =  left_jk[2][1] + right_iaic[2][1];
        index[42] =  left_jk[2][2] + right_iaic[2][2];
        index[43] =  left_jk[2][3] + right_iaic[2][3];
        temp11 = simd_set_floatv4(cpe_ppg_ee_now[index[40]],
                                 cpe_ppg_ee_now[index[41]],
                                 cpe_ppg_ee_now[index[42]],
                                 cpe_ppg_ee_now[index[43]]); 
        e3s[2] = temp11 * fien;
        
        
        index[44] =  left_jk[3][0] + right_iaic[2][0];
        index[45] =  left_jk[3][1] + right_iaic[2][1];
        index[46] =  left_jk[3][2] + right_iaic[2][2];
        index[47] =  left_jk[3][3] + right_iaic[2][3];
        temp12 = simd_set_floatv4(cpe_ppg_ee_now[index[44]],
                                 cpe_ppg_ee_now[index[45]],
                                 cpe_ppg_ee_now[index[46]],
                                 cpe_ppg_ee_now[index[47]]); 
        e4s[2] = temp12 * fien1;
        
        
        index[48] =  left_jk[0][0] + right_iaic[3][0];
        index[49] =  left_jk[0][1] + right_iaic[3][1];
        index[50] =  left_jk[0][2] + right_iaic[3][2];
        index[51] =  left_jk[0][3] + right_iaic[3][3];
        temp13 = simd_set_floatv4(cpe_ppg_ee_now[index[48]],
                                 cpe_ppg_ee_now[index[49]],
                                 cpe_ppg_ee_now[index[50]],
                                 cpe_ppg_ee_now[index[51]]); 
        e1s[3] = temp13 * fien;
        
        
        index[52] =  left_jk[1][0] + right_iaic[3][0];
        index[53] =  left_jk[1][1] + right_iaic[3][1];
        index[54] =  left_jk[1][2] + right_iaic[3][2];
        index[55] =  left_jk[1][3] + right_iaic[3][3];
        temp14 = simd_set_floatv4(cpe_ppg_ee_now[index[52]],
                                cpe_ppg_ee_now[index[53]],
                                cpe_ppg_ee_now[index[54]],
                                cpe_ppg_ee_now[index[55]]); 
        e2s[3] = temp14 * fien1;
        
        
        index[56] =  left_jk[2][0] + right_iaic[3][0];
        index[57] =  left_jk[2][1] + right_iaic[3][1];
        index[58] =  left_jk[2][2] + right_iaic[3][2];
        index[59] =  left_jk[2][3] + right_iaic[3][3];
        temp15 = simd_set_floatv4(cpe_ppg_ee_now[index[56]],
                                cpe_ppg_ee_now[index[57]],
                                cpe_ppg_ee_now[index[58]],
                                cpe_ppg_ee_now[index[59]]); 
        e3s[3] = temp15 * fien;
        
        
        index[60] =  left_jk[3][0] + right_iaic[3][0];
        index[61] =  left_jk[3][1] + right_iaic[3][1];
        index[62] =  left_jk[3][2] + right_iaic[3][2];
        index[63] =  left_jk[3][3] + right_iaic[3][3];
        temp16 = simd_set_floatv4(cpe_ppg_ee_now[index[60]],
                                cpe_ppg_ee_now[index[61]],
                                cpe_ppg_ee_now[index[62]],
                                cpe_ppg_ee_now[index[63]]); 
        e4s[3] = temp16 * fien1;
         
        for (i = 0; i < 4; i++)
            es[i] = bilinear_interpolation_qr(e1s[i], e2s[i], e3s[i], e4s[i], q_4, r_4);

        exxyy = bilinear_interpolation_qr(
            es[0], es[1], es[2], es[3],
            q, r
        );
        float _exxyy[4];
        simd_store(exxyy,&_exxyy[0]);

        for(i = 0; i < 4; i++)
        {
            cpe_ppg_e[cpe_ppg_jv[jk + i] + cpe_ppg_kv[jk + i]] = max(_exxyy[i], _small);
        }
        
        //packed_data_ptr += PPG_BI_DATA;
        values_ptr = values_ptr + 24;
        idxs_ptr = idxs_ptr + 20;
    }
    
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
        athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
    }
}

void two_step_write_data(void* mpe_data_ptr, void* cpe_data_ptr, int size, int my_core, int write_flag)
{
//    put_reply = 0;
//    athread_put(
//        PE_MODE, cpe_data_ptr, mpe_data_ptr,
//        size, (void*)(&put_reply), 0, 0
//    );
//    while (put_reply != 1);
//    athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);

    int wcid;
    for (wcid = 1; wcid >= 0; wcid--)
    {
        if (my_core % 2 == wcid && write_flag != 0)
        {
            put_reply = 0;
            athread_put(
                PE_MODE, cpe_data_ptr, mpe_data_ptr,
                size, (void*)(&put_reply), 0, 0
            );
            while (put_reply != 1);
        }
        athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
    }
}
/*
void cpe_setspec2()
{
    const float gama=3.3, sq3=1.7320508075688772935274463415059;
    float vx, vy, w, xj, xj0, arlfa, wsj, wkj, theta0, sinth, costh, wk0, wf0, ws0, wl, sigma, alpha;
    int j, k;
    int index_tmp;

    vx = cpe_ips_packed_data[WX_OFFSET];
    vy = cpe_ips_packed_data[WY_OFFSET];

    w = sqrt(vx * vx + vy * vy);

    if(w <= 0.0f)
        w = 0.9f;

    xj0 = 200.0 * 1000.0;
    xj = _g * xj0 / (w * w);
    arlfa = (0.076 * (powf(xj, -0.4f))) / _pi;
    wsj = 22. * (powf(xj, -0.33f)) * _g / w;
    wkj = (wsj * wsj) / _g;

    for(j = 1; j <= _jl; j++)
    {
        costh = cpe_ips_cos_theta[j - 1];
        sinth = cpe_ips_sin_theta[j - 1];
        for(k = 1; k <= _kl; k++)
        {
            wk0 = cpe_ips_wk[k - 1];
            //index_tmp = (ic - _iys) * (_kldp1 * ix_size) + (ia - _ixs) * _kldp1;
            wf0 = cpe_ips_wf[k - 1];
            ws0 = _zpi * wf0;
            wl = vx * costh + vy * sinth;

            if(wl > 0.)
            {
                if(ws0 <= wsj)
                {
                    sigma = 0.07;
                }
                else
                {
                    sigma = 0.09;
                }
                alpha = arlfa / powf(wk0, 4) * expf(-1.25*(wkj/wk0)*(wkj/wk0))
                        * powf(gama, exp(-0.5 * powf(((1.0 - ws0 / wsj) / sigma), 2)))
                        * powf((wl / w), 2);
            }
            else
                alpha = 0.0;

            index_tmp = k - 1 + (j - 1) * _kl;

            //if(my_rank == 0 && k == 10 && j == 10)
            //{
            //    printf("c:\n  k, j, ia, ic: %d, %d, %d, %d\n", k, j, ia, ic);

            //}

            cpe_ips_ee[index_tmp] = max(alpha, _small);
        }
    }
}*/

void copy_ips_permanent_arrays(int my_core)
{
    // NOTICE : ???????? get ?? while????? cpe_const_values[] ??????????????
    // ikp, ikp1, ikm, ikm1, wks17, grolim : kl
    get_reply = 0;
    athread_get(
        PE_MODE, param.host_ikp, &cpe_ips_ikp[0],
        sizeof(int) * _kl, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    get_reply = 0;
    athread_get(
        PE_MODE, param.host_ikp1, &cpe_ips_ikp1[0],
        sizeof(int) * _kl, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    get_reply = 0;
    athread_get(
        PE_MODE, param.host_ikm, &cpe_ips_ikm[0],
        sizeof(int) * _kl, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    get_reply = 0;
    athread_get(
        PE_MODE, param.host_ikm1, &cpe_ips_ikm1[0],
        sizeof(int) * _kl, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    get_reply = 0;
    athread_get(
        PE_MODE, param.host_wks17, &cpe_ips_wks17[0],
        sizeof(float) * _kl, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    get_reply = 0;
    athread_get(
        PE_MODE, param.host_grolim, &cpe_ips_grolim[0],
        sizeof(float) * _kl, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    // wp, wm : 2 * 2 * kl
    get_reply = 0;
    athread_get(
        PE_MODE, param.host_wp, &cpe_ips_wp[0],
        sizeof(float) * _kl * 2 * 2, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    get_reply = 0;
    athread_get(
        PE_MODE, param.host_wm, &cpe_ips_wm[0],
        sizeof(float) * _kl * 2 * 2, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    // wkh, wk, dwk : kldp1
    get_reply = 0;
    athread_get(
        PE_MODE, param.host_wkh, &cpe_ips_wkh[0],
        sizeof(float) * _kldp1, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    get_reply = 0;
    athread_get(
        PE_MODE, param.host_wk, &cpe_ips_wk[0],
        sizeof(float) * _kldp1, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    get_reply = 0;
    athread_get(
        PE_MODE, param.host_dwk, &cpe_ips_dwk[0],
        sizeof(float) * _kldp1, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    // jp1, jp2, jm1, jm2 : 2 * jl
    get_reply = 0;
    athread_get(
        PE_MODE, param.host_jp1, &cpe_ips_jp1[0],
        sizeof(int) * 2 * _jl, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    get_reply = 0;
    athread_get(
        PE_MODE, param.host_jp2, &cpe_ips_jp2[0],
        sizeof(int) * 2 * _jl, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    get_reply = 0;
    athread_get(
        PE_MODE, param.host_jm1, &cpe_ips_jm1[0],
        sizeof(int) * 2 * _jl, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    get_reply = 0;
    athread_get(
        PE_MODE, param.host_jm2, &cpe_ips_jm2[0],
        sizeof(int) * 2 * _jl, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    // cos_theta, sin_theta : jlp1
    get_reply = 0;
    athread_get(
        PE_MODE, param.host_cos_theta, &cpe_ips_cos_theta[0],
        sizeof(int) * _jlp1, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    get_reply = 0;
    athread_get(
        PE_MODE, param.host_sin_theta, &cpe_ips_sin_theta[0],
        sizeof(int) * _jlp1, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    // ????????????Χ???????
    get_reply = 0;
    athread_get(
        PE_MODE, param.host_packed_ranges, &cpe_xy_ranges[0],
        sizeof(int) * NUM_HPR, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    get_reply = 0;
    athread_get(
        PE_MODE, param.host_packed_consts, &cpe_const_values[0],
        sizeof(float) * NUM_HPD, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    cpe_ixs = cpe_xy_ranges[IXS_OFFSET];
    cpe_ixl = cpe_xy_ranges[IXL_OFFSET];
    cpe_iys = cpe_xy_ranges[IYS_OFFSET];
    cpe_iyl = cpe_xy_ranges[IYL_OFFSET];

    cpe_ix_size   = cpe_ixl - cpe_ixs + 1;
    cpe_iy_size   = cpe_iyl - cpe_iys + 1;
    cpe_ixiy_size = cpe_ix_size * cpe_iy_size;
    cpe_dim4      = _kl * _jl * cpe_ix_size;
    cpe_dim3      = _kl * _jl;
    cpe_dim2      = _kl;

    cpe_deltt5 = cpe_const_values[DELTT5_OFFSET];
    cpe_deltt  = cpe_const_values[ DELTT_OFFSET];
    cpe_cong   = cpe_const_values[  CONG_OFFSET];
    cpe_al31   = cpe_const_values[  AL31_OFFSET];
    cpe_al21   = cpe_const_values[  AL21_OFFSET];
    cpe_al11   = cpe_const_values[  AL11_OFFSET];
    cpe_al13   = cpe_const_values[  AL13_OFFSET];
    cpe_al23   = cpe_const_values[  AL23_OFFSET];

    // WARNING : ?????????????????????????? c_acce_adapter.c?ж??????????????????
    if (param.host_mpi_rank == 0 && my_core == 0)
    {
        printf("CPE recv range : %d %d %d %d\n", cpe_ixs, cpe_ixl, cpe_iys, cpe_iyl);
        printf("CPE recv const : %e %e %e %e %e %e %e %e \n", cpe_deltt5, cpe_deltt, cpe_cong, cpe_al11, cpe_al21, cpe_al31, cpe_al13, cpe_al23);
        printf("CPE recv ptr   : %x %x %x %x %x %x\n", param.host_pexx_output_buf, param.host_e, param.host_ee, param.host_wf, param.host_ccg, param.host_ssbos);
    }
}

void load_mpe_ips_e_wf_ccg()
{
    int wf_ccg_offset, e_offset;

    //cpe_ips_ic = ia_ptr[compute_unit];
    //cpe_ips_ia = ic_ptr[compute_unit];//????????????

    wf_ccg_offset = (cpe_ips_ic - cpe_iys) * (_kldp1 * cpe_ix_size) + (cpe_ips_ia - cpe_ixs) * _kldp1;
    e_offset = (cpe_ips_ic - cpe_iys) * cpe_dim4 + (cpe_ips_ia - cpe_ixs) * cpe_dim3;

    // if (compute_unit == 0) printf("CU0, ia ic = %d %d, offsets = %d %d\n", cpe_ips_ia, cpe_ips_ic, wf_ccg_offset, e_offset);

    get_reply = 0;
    athread_get(
        PE_MODE, param.host_wf + wf_ccg_offset, &cpe_ips_wf[0],
        sizeof(float) * _kldp1, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    get_reply = 0;
    athread_get(
        PE_MODE, param.host_ccg + wf_ccg_offset, &cpe_ips_ccg[0],
        sizeof(float) * _kldp1, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    get_reply = 0;
    athread_get(
        PE_MODE, param.host_ssbos + wf_ccg_offset, &cpe_ips_ssbos[0],
        sizeof(float) * _kldp1, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);
}

void load_mpe_ips_packed_data(int compute_unit, float* packed_data_ptr)
{
    int ips_packed_data_offset = compute_unit * PPG_SV_DATA;
    get_reply = 0;
    athread_get(
        PE_MODE, packed_data_ptr + ips_packed_data_offset, &cpe_ips_packed_data,
        sizeof(float) * PPG_SV_DATA, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);
}

void cpe_c_mean2()
{
    int j, k, k1;
    floatv4 dwkk, wfk, wfk1, wsk, wsk1, wkk, wkk1, swkk, swkk1, ekj, ekj1;
    floatv4 zpiv4 = _zpi;
    float _aev4[4], _awfv4[4], _asiv4[4], _awkv4[4], _arkv4[4];
    floatv4 aev4=0.0, awfv4=0.0, asiv4=0.0, awkv4=0.0, arkv4=0.0;
    floatv4 inner_factor1, inner_factor2;
    int e_index, e_index2, cond_k1, cond_k2;
    double tmpd;
    float tmpf;
    
    thread_local_buff.ae  = 0.0;
    thread_local_buff.asi = 0.0;
    thread_local_buff.awf = 0.0;
    thread_local_buff.awk = 0.0;
    thread_local_buff.ark = 0.0;

    for (k = 1; k <= _kld; k++)
    {
        k1   = k + 1;
        dwkk = cpe_ips_dwk[k - 1];
        wfk  = cpe_ips_wf[k - 1];
        wfk1 = cpe_ips_wf[k1 - 1];
        tmpd = _zpi * cpe_ips_wf[k - 1];
        //tmpd = 1.0 / tmpd;
        tmpf = tmpd;
        wsk  = tmpf;
        tmpd = _zpi * cpe_ips_wf[k];
        //tmpd = 1.0 / tmpd;
        tmpf = tmpd;
        wsk1  = tmpf;
        wkk  = cpe_ips_wk[k - 1];
        wkk1 = cpe_ips_wk[k1 - 1];
        swkk = cpe_ips_sqrt_div_wk[k - 1];
        swkk1 = cpe_ips_sqrt_div_wk[k1 - 1];

        if (k < _kl)
        {
            inner_factor1 = 1.0;
            inner_factor2 = 1.0;
            cond_k1 = k;
            cond_k2 = k1;
        }
        else
        {
            inner_factor1 = cpe_ips_wkh[k - _kl + 1 - 1];
            inner_factor2 = cpe_ips_wkh[k - _kl + 2 - 1];
            cond_k1 = _kl;
            cond_k2 = _kl;
        }

        
        for (j = 1; j <= _jl; j+=4)
        {
            e_index  = (cond_k1 - 1) * _jl + j - 1;
            e_index2  = (cond_k2 - 1) * _jl + j - 1;

            ekj = simd_set_floatv4(cpe_ips_e[e_index]
                               ,cpe_ips_e[e_index + 1]
                               ,cpe_ips_e[e_index + 2]
                               ,cpe_ips_e[e_index + 3]);
            ekj1 = simd_set_floatv4(cpe_ips_e[e_index2]
                               ,cpe_ips_e[e_index2 + 1]
                               ,cpe_ips_e[e_index2 + 2]
                               ,cpe_ips_e[e_index2 + 3]);
            ekj = ekj * inner_factor1;
            ekj1 = ekj1 * inner_factor2;

            aev4 = aev4 + (ekj + ekj1) * dwkk;
            awfv4 = awfv4 + (ekj * wfk + ekj1 * wfk1) * dwkk;
            asiv4 = asiv4 + (ekj / wsk + ekj1 / wsk1) * dwkk;
            awkv4 = awkv4 + (ekj * wkk + ekj1 * wkk1) * dwkk;
            arkv4 = arkv4 + (ekj * swkk + ekj1 * swkk1) * dwkk;
        }

    }
    simd_store(aev4,&_aev4[0]);
    simd_store(awfv4,&_awfv4[0]);
    simd_store(asiv4,&_asiv4[0]);
    simd_store(awkv4,&_awkv4[0]);
    simd_store(arkv4,&_arkv4[0]);
    thread_local_buff.ae  =  _aev4[0] + _aev4[1] + _aev4[2] + _aev4[3];
    thread_local_buff.awf =  _awfv4[0] + _awfv4[1] + _awfv4[2] + _awfv4[3];
    thread_local_buff.asi =  _asiv4[0] + _asiv4[1] + _asiv4[2] + _asiv4[3];
    thread_local_buff.awk =  _awkv4[0] + _awkv4[1] + _awkv4[2] + _awkv4[3];
    thread_local_buff.ark =  _arkv4[0] + _arkv4[1] + _arkv4[2] + _arkv4[3];

    thread_local_buff.asi =  thread_local_buff.ae / thread_local_buff.asi;
    thread_local_buff.awf = thread_local_buff.awf /  thread_local_buff.ae;
    thread_local_buff.awk = thread_local_buff.awk /  thread_local_buff.ae;
    thread_local_buff.ark =  thread_local_buff.ae / thread_local_buff.ark;
    thread_local_buff.ark = thread_local_buff.ark * thread_local_buff.ark;
}

void cpe_c_get_thresholds()
{
    int i, kpmt, kakt, ks1, ks;
    float vx, vy, ww, wkpmt, wakt;

    vx = cpe_ips_packed_data[WX_OFFSET];
    vy = cpe_ips_packed_data[WY_OFFSET];
    ww = vx * vx + vy * vy;
    if (ww == 0.0) ww = 0.5;
    wkpmt = _cksp * _gc2 / ww;
    kpmt = (int) (log10(wkpmt / cpe_ips_wk[0]) / _alog10pwk);
    kpmt += 1;

    wakt = _cksa * thread_local_buff.awk;
    kakt = (int) (log10(wakt / cpe_ips_wk[0]) / _alog10pwk);
    kakt += 1;

    ks1 = max(kpmt, kakt);
    ks  = min(ks1, _kl);

    thread_local_buff.ks0   = ks;

    for (i = 1; i <= ks; i++) thread_local_buff.fconst0[i - 1] = 1.0;
    for (i = ks + 1; i <= _kl; i++) thread_local_buff.fconst0[i - 1] = 0.0;
}

void cpe_c_snonlin()
{
    float xx;
    floatv4 cwks17, ffacp, ffacp1, fcen, up, up1, um, um1, sap, sam;;
    floatv4 wp11, wp12, wp21, wp22, wm11, wm12, wm21, wm22;
    floatv4 wp112, wp122, wp212, wp222, wm112, wm122, wm212, wm222;

    floatv4 ea1, ea2, ea3, ea4, ea5, ea6, ea7, ea8;
    floatv4 eij, eij2, zua, ead1, ead2, ad, adp, adm, delad, deladp, deladm;
    floatv4 two = 2.0, cpe_al31v4 = cpe_al31, cpe_al21v4 = cpe_al21;
    floatv4 one = 1.0;
    floatv4 cpe_al11v4 = cpe_al11, cpe_al13v4 = cpe_al13, cpe_al23v4 = cpe_al23;
    floatv4 cpe_al31v4_inv = _inv_al31, cpe_al21v4_inv = _inv_al21;
    floatv4 cpe_al11v4_inv = _inv_al11, cpe_al13v4_inv = _inv_al13, cpe_al23v4_inv = _inv_al23;
    floatv4 cpe_al31v4_inv_two = cpe_al31v4_inv * two;
    floatv4 cpe_al13v4_inv_fcen, cpe_al23v4_inv_fcen;
    floatv4 cpe_al31v4_inv_two_minus = -cpe_al31v4_inv_two;

    floatv4 ans[9][2];
    float   _ans[9][2][4];

    floatv4 tmp_fv4_1, tmp_fv4_2, tmp_fv4_3, tmp_fv4_4;
    floatv4 two_fcen;
    floatv4 wp11_f, wp12_f, wp21_f, wp22_f;

    int kh, k, mr, i, j, im, im1, kp, kp1, kp2, kp3;
    int j11[4], j12[4], j21[4], j22[4], jpjm_index[4], se_index;
    int jse11[4], jse12[4], jse21[4], jse22[4];

    xx = 0.75 * cpe_ips_packed_data[D_OFFSET] * thread_local_buff.awk;
    if (xx < 0.5) xx = 0.5;
    thread_local_buff.enh = 1.0 + (5.5 / xx) * (1.0 - 0.833 * xx) * exp(-1.25 * xx);
    kh = 0;

    int e_index;
    int jpms_index1, jpms_index2;

    for (k = 1; k <= _kl; k++)
    {
        wp11 = cpe_ips_wp[_wpwm_index(k, 1, 1)];
        wp12 = cpe_ips_wp[_wpwm_index(k, 1, 2)];
        wp21 = cpe_ips_wp[_wpwm_index(k, 2, 1)];
        wp22 = cpe_ips_wp[_wpwm_index(k, 2, 2)];
        wm11 = cpe_ips_wm[_wpwm_index(k, 1, 1)];
        wm12 = cpe_ips_wm[_wpwm_index(k, 1, 2)];
        wm21 = cpe_ips_wm[_wpwm_index(k, 2, 1)];
        wm22 = cpe_ips_wm[_wpwm_index(k, 2, 2)];

        wp112 = wp11 * wp11;
        wp122 = wp12 * wp12;
        wp212 = wp21 * wp21;
        wp222 = wp22 * wp22;
        wm112 = wm11 * wm11;
        wm122 = wm12 * wm12;
        wm212 = wm21 * wm21;
        wm222 = wm22 * wm22;

        ffacp  = 1.0;
        ffacp1 = 1.0;

        im     = cpe_ips_ikm [k - 1];
        im1    = cpe_ips_ikm1[k - 1];
        kp     = cpe_ips_ikp [k - 1];
        kp1    = cpe_ips_ikp1[k - 1];
        kp2    = cpe_ips_ikp [k - 1];
        kp3    = cpe_ips_ikp1[k - 1];
        cwks17 = cpe_cong * cpe_ips_wks17[k - 1];

        if (kp >= _kl)
        {
            kh  = kh + 1;
            kp2 = _kl + 1;
            if (kp == _kl) kp2 = _kl;
            kp  = _kl;
            kp1 = _kl;
            kp3 = _kl + 1;
            ffacp  = cpe_ips_wkh[kh - 1];
            ffacp1 = cpe_ips_wkh[kh + 1 - 1];
        }

        wp11_f = wp11 * ffacp;
        wp12_f = wp12 * ffacp;
        wp21_f = wp21 * ffacp1;
        wp22_f = wp22 * ffacp1;

        fcen = thread_local_buff.fconst0[k - 1] * thread_local_buff.enh * cwks17;
        cpe_al13v4_inv_fcen = cpe_al13v4_inv * fcen; 
        cpe_al23v4_inv_fcen = cpe_al23v4_inv * fcen;
        two_fcen = fcen * two;

        kp = _jl * (kp - 1);
        kp1 = _jl * (kp1 - 1);
        kp2 = _jl * (kp2 - 1);
        kp3 = _jl * (kp3 - 1);
        im = _jl * (im - 1);
        im1 = _jl * (im1 - 1);

        jpms_index1 = 0;
        jpms_index2 = 0;

        for (mr = 1; mr <= 2; mr++)
        {
            for (j = 1; j <= _jl; j+=4)
            {
                //RPCC(stcc2);
                e_index = (k - 1) * _jl + j - 1;
                eij = simd_set_floatv4( cpe_ips_e[e_index + 0]
                                       ,cpe_ips_e[e_index + 1]
                                       ,cpe_ips_e[e_index + 2]
                                       ,cpe_ips_e[e_index + 3]);
                //if (eij < 1e-20) continue;
                /*
                for (i = 0; i < 4; i++)
                {
                    jpjm_index[i] = _jpjm_index(mr, j+i);
                    j11[i] = (cpe_ips_jp1[jpjm_index[i]] - 1);
                    j12[i] = (cpe_ips_jp2[jpjm_index[i]] - 1);
                    j21[i] = (cpe_ips_jm1[jpjm_index[i]] - 1);
                    j22[i] = (cpe_ips_jm2[jpjm_index[i]] - 1);
                    jse11[i] = (cpe_ips_jp1[jpjm_index[i]] - 1);
                    jse12[i] = (cpe_ips_jp2[jpjm_index[i]] - 1);
                    jse21[i] = (cpe_ips_jm1[jpjm_index[i]] - 1);
                    jse22[i] = (cpe_ips_jm2[jpjm_index[i]] - 1);
                }*/

                for (i = 0; i < 4; i++)
                {
                    //jpjm_index[i] = _jpjm_index(mr, j+i);

                    jse11[i] = j11[i] = cpe_ips_jpms[jpms_index1];
                    jse12[i] = j12[i] = cpe_ips_jpms[jpms_index1 + 1];
                    jse21[i] = j21[i] = cpe_ips_jpms[jpms_index1 + 2];
                    jse22[i] = j22[i] = cpe_ips_jpms[jpms_index1 + 3];

                    jpms_index1 += 4;
                    
                }

                ea1 = simd_set_floatv4( cpe_ips_e[j11[0] + kp]
                                       ,cpe_ips_e[j11[1] + kp]
                                       ,cpe_ips_e[j11[2] + kp]
                                       ,cpe_ips_e[j11[3] + kp]);
                ea2 = simd_set_floatv4( cpe_ips_e[j12[0] + kp]
                                       ,cpe_ips_e[j12[1] + kp]
                                       ,cpe_ips_e[j12[2] + kp]
                                       ,cpe_ips_e[j12[3] + kp]);
                ea3 = simd_set_floatv4( cpe_ips_e[j11[0] + kp1]
                                       ,cpe_ips_e[j11[1] + kp1]
                                       ,cpe_ips_e[j11[2] + kp1]
                                       ,cpe_ips_e[j11[3] + kp1]);
                ea4 = simd_set_floatv4( cpe_ips_e[j12[0] + kp1]
                                       ,cpe_ips_e[j12[1] + kp1]
                                       ,cpe_ips_e[j12[2] + kp1]
                                       ,cpe_ips_e[j12[3] + kp1]);
                ea5 = simd_set_floatv4( cpe_ips_e[j21[0] + im]
                                       ,cpe_ips_e[j21[1] + im]
                                       ,cpe_ips_e[j21[2] + im]
                                       ,cpe_ips_e[j21[3] + im]);
                ea6 = simd_set_floatv4( cpe_ips_e[j22[0] + im]
                                       ,cpe_ips_e[j22[1] + im]
                                       ,cpe_ips_e[j22[2] + im]
                                       ,cpe_ips_e[j22[3] + im]);
                ea7 = simd_set_floatv4( cpe_ips_e[j21[0] + im1]
                                       ,cpe_ips_e[j21[1] + im1]
                                       ,cpe_ips_e[j21[2] + im1]
                                       ,cpe_ips_e[j21[3] + im1]);
                ea8 = simd_set_floatv4( cpe_ips_e[j22[0] + im1]
                                       ,cpe_ips_e[j22[1] + im1]
                                       ,cpe_ips_e[j22[2] + im1]
                                       ,cpe_ips_e[j22[3] + im1]);

                //RPCC(edcc2);
                //out_flag[IPF_NONLIN_1] += edcc2 - stcc2;
                //stcc2 = edcc2;

                tmp_fv4_1 = wp11_f * ea1;
                up = tmp_fv4_1 + wp12_f * ea2;
                tmp_fv4_2 = wp21_f * ea3 + up;
                sap = tmp_fv4_2 + wp22_f * ea4;

                tmp_fv4_3 = wm11 * ea5;
                um = tmp_fv4_3 + wm12 * ea6;
                tmp_fv4_4 = wm21 * ea7 + um;
                sam = tmp_fv4_4 + wm22 * ea8;

                eij2 = eij * eij;
                zua  = eij * cpe_al31v4_inv_two;
                tmp_fv4_1 = sap * cpe_al11v4_inv;
                ead1 = tmp_fv4_1 + sam * cpe_al21v4_inv;

                tmp_fv4_2 = sap * sam;
                ead2 = tmp_fv4_2 * cpe_al31v4_inv_two_minus;

                tmp_fv4_3 = eij2 * ead1;
                ad = tmp_fv4_3 + ead2 * eij;
                adp  = ad * cpe_al13v4_inv_fcen;
                adm  = ad * cpe_al23v4_inv_fcen;

                tmp_fv4_4 = eij * two;
                delad = ead1 * tmp_fv4_4 + ead2;

                tmp_fv4_1 = zua * sam;
                tmp_fv4_2 = eij2 * cpe_al11v4_inv - tmp_fv4_1;
                deladp = tmp_fv4_2 * cpe_al13v4_inv_fcen;
                tmp_fv4_3 = zua * sap;
                tmp_fv4_4 = eij2 * cpe_al21v4_inv - tmp_fv4_3;
                deladm = tmp_fv4_4 * cpe_al23v4_inv_fcen;

                //RPCC(edcc2);
                //out_flag[IPF_NONLIN_2] += edcc2 - stcc2;
                //stcc2 = edcc2;


                ans[0][0] = two_fcen * ad;    simd_store(ans[0][0],(void*)&_ans[0][0]);
                ans[0][1] = two_fcen * delad; simd_store(ans[0][1],(void*)&_ans[0][1]);

                ans[1][0] = adp * wp11;       simd_store(ans[1][0],(void*)&_ans[1][0]);
                ans[1][1] = deladp * wp112;   simd_store(ans[1][1],(void*)&_ans[1][1]);

                ans[2][0] = adp * wp12;       simd_store(ans[2][0],(void*)&_ans[2][0]);
                ans[2][1] = deladp * wp122;   simd_store(ans[2][1],(void*)&_ans[2][1]);

                ans[3][0] = adp * wp21;       simd_store(ans[3][0],(void*)&_ans[3][0]);
                ans[3][1] = deladp * wp212;   simd_store(ans[3][1],(void*)&_ans[3][1]);

                ans[4][0] = adp * wp22;       simd_store(ans[4][0],(void*)&_ans[4][0]);
                ans[4][1] = deladp * wp222;   simd_store(ans[4][1],(void*)&_ans[4][1]);

                ans[5][0] = adm * wm11;       simd_store(ans[5][0],(void*)&_ans[5][0]);
                ans[5][1] = deladm * wm112;   simd_store(ans[5][1],(void*)&_ans[5][1]);

                ans[6][0] = adm * wm12;       simd_store(ans[6][0],(void*)&_ans[6][0]);
                ans[6][1] = deladm * wm122;   simd_store(ans[6][1],(void*)&_ans[6][1]);

                ans[7][0] = adm * wm21;       simd_store(ans[7][0],(void*)&_ans[7][0]);
                ans[7][1] = deladm * wm212;   simd_store(ans[7][1],(void*)&_ans[7][1]);

                ans[8][0] = adm * wm22;       simd_store(ans[8][0],(void*)&_ans[8][0]);
                ans[8][1] = deladm * wm222;   simd_store(ans[8][1],(void*)&_ans[8][1]);

                //RPCC(edcc2);
                //out_flag[IPF_NONLIN_3] += edcc2 - stcc2;
                //stcc2 = edcc2;

                for (i = 0; i < 4; i++)
                {
                    se_index = (j + i - 1) + (k - 1) * _jl;
                    thread_local_buff.se[se_index]  -= _ans[0][0][i];
                    thread_local_buff.dse[se_index] -= _ans[0][1][i];

                    se_index = kp2 + jse11[i];
                    thread_local_buff.se[se_index]  += _ans[1][0][i];
                    thread_local_buff.dse[se_index] += _ans[1][1][i];

                    se_index = kp2 + jse12[i];
                    thread_local_buff.se[se_index]  += _ans[2][0][i];
                    thread_local_buff.dse[se_index] += _ans[2][1][i];

                    se_index = kp3 + jse11[i];
                    thread_local_buff.se[se_index]  += _ans[3][0][i];
                    thread_local_buff.dse[se_index] += _ans[3][1][i];

                    se_index = kp3 + jse12[i];
                    thread_local_buff.se[se_index]  += _ans[4][0][i];
                    thread_local_buff.dse[se_index] += _ans[4][1][i];

                    se_index = im +jse21[i];
                    thread_local_buff.se[se_index]  += _ans[5][0][i];
                    thread_local_buff.dse[se_index] += _ans[5][1][i];

                    se_index = im + jse22[i];
                    thread_local_buff.se[se_index]  += _ans[6][0][i];
                    thread_local_buff.dse[se_index] += _ans[6][1][i];

                    se_index = im1 + jse21[i];
                    thread_local_buff.se[se_index]  += _ans[7][0][i];
                    thread_local_buff.dse[se_index] += _ans[7][1][i];

                    se_index = im1 + jse22[i];
                    thread_local_buff.se[se_index]  += _ans[8][0][i];
                    thread_local_buff.dse[se_index] += _ans[8][1][i];
                }

                //RPCC(edcc2);
                //out_flag[IPF_NONLIN_4] += edcc2 - stcc2;
            }
        }
    }
}

void cpe_c_sinput()
{
    int j, k, ks, sein_index, se_index, e_index;
    floatv4 vx, vy, cd, costh, sinth, wl, wlstar, wk0, wf0, ws0, bett, beta;
    float cdv1;
    floatv4 zpi_v4, esv4;
    floatv4 v_28, zero_v4;

    ks = thread_local_buff.ks0;
    vx = cpe_ips_packed_data[WX_OFFSET];
    vy = cpe_ips_packed_data[WY_OFFSET];
    cdv1 = (0.8 + 0.065 * cpe_ips_packed_data[W_OFFSET]) * 0.001;
    cdv1 = sqrt(cdv1);
    cd = cdv1;
    v_28 = 28.0;
    zero_v4 = 0.0;
    //transpose

    bett = _beta10;
    zpi_v4 = _zpi;

    for (k = 1; k <= ks; k++)
    {
        wk0 = cpe_ips_wk[k - 1];
        wf0 = cpe_ips_wf[k - 1];
        ws0 = zpi_v4 * wf0;

        for (j = 1; j <= _jl; j += 4)
        {
            sein_index = _sexx_index(k, j);
            se_index = e_index = (k - 1) * _jl + j - 1;
            esv4 = simd_set_floatv4(cpe_ips_e[e_index]
                               ,cpe_ips_e[e_index + 1]
                               ,cpe_ips_e[e_index + 2]
                               ,cpe_ips_e[e_index + 3]);

            sinth  = simd_set_floatv4(cpe_ips_sin_theta[j - 1],cpe_ips_sin_theta[j],
                                        cpe_ips_sin_theta[j + 1],cpe_ips_sin_theta[j + 2]);

            costh  = simd_set_floatv4(cpe_ips_cos_theta[j - 1],cpe_ips_cos_theta[j],
                                        cpe_ips_cos_theta[j + 1],cpe_ips_cos_theta[j + 2]);

            wl     = vx * costh + vy * sinth;
            wlstar = wl * cd;
            bett = _beta10;
            beta = bett * (wk0 * v_28 * wlstar - ws0);
            beta = __builtin_sw64_sllt(beta, zero_v4, beta);
            esv4 = esv4 * beta;
            simd_store(beta, &cpe_tmp_fv[0]);
            simd_store(esv4, &cpe_tmp_fv1[0]);

            thread_local_buff.sein[sein_index + 0 * _kl] = cpe_tmp_fv[0];
            thread_local_buff.sein[sein_index + 1 * _kl] = cpe_tmp_fv[1];
            thread_local_buff.sein[sein_index + 2 * _kl] = cpe_tmp_fv[2];
            thread_local_buff.sein[sein_index + 3 * _kl] = cpe_tmp_fv[3];
            thread_local_buff.dse[se_index + 0] += cpe_tmp_fv[0];
            thread_local_buff.dse[se_index + 1] += cpe_tmp_fv[1];
            thread_local_buff.dse[se_index + 2] += cpe_tmp_fv[2];
            thread_local_buff.dse[se_index + 3] += cpe_tmp_fv[3];
            thread_local_buff.se[se_index + 0] += cpe_tmp_fv1[0];
            thread_local_buff.se[se_index + 1] += cpe_tmp_fv1[1];
            thread_local_buff.se[se_index + 2] += cpe_tmp_fv1[2];
            thread_local_buff.se[se_index + 3] += cpe_tmp_fv1[3];
        }
    }
/*
    for (j = 1; j <= _jl; j++)
    {
        sinth  = cpe_ips_sin_theta[j - 1];
        costh  = cpe_ips_cos_theta[j - 1];
        wl     = vx * costh + vy * sinth;
        wlstar = wl * cd;

        sein_index = _sexx_index(0, j);
        se_index = _se_index(0, j);
        e_index  = (j - 1) * _kl - 1;

        for (k = 1; k <= ks; k++)
        {
            sein_index++;
            se_index++;
            e_index++;
            wk0 = cpe_ips_wk[k - 1];
            wf0 = cpe_ips_wf[k - 1];
            ws0 = _zpi * wf0;
            bett = _beta10;
            beta = bett * (wk0 * 28.0 * wlstar - ws0);
            beta = max(0.0, beta);
            thread_local_buff.sein[sein_index] = beta * cpe_ips_e[e_index];
            thread_local_buff.se[se_index] += thread_local_buff.sein[sein_index];
            thread_local_buff.dse[se_index] += beta;
        }
    }*/

}

void cpe_c_sdissip()
{
    int j, k, ks, se_index, seds_index, e_index;
    float arkss, ekspm, sds, ssds, tmpf;

    ks = thread_local_buff.ks0;
    arkss = thread_local_buff.ark;
    ekspm = thread_local_buff.ae * arkss * arkss / 0.0030162;
    sds = _d1 * thread_local_buff.asi / arkss * sqrt(ekspm) * exp(-_d2 * 0.64 / ekspm);

    for (j = 1; j <= _jl; j++)
    {
        se_index = j - 1;
        //seds_index = _sexx_index(0, j);
        e_index = j - 1;
        seds_index = _sexx_index(0, j);
        for (k = 1; k <= ks; k++)
        {
            seds_index++;
            ssds = -_ads * sds * cpe_ips_wk[k - 1];
            tmpf = ssds * cpe_ips_e[e_index];
            thread_local_buff.seds[seds_index] = tmpf;
            thread_local_buff.se[se_index]  += tmpf;
            thread_local_buff.dse[se_index] += ssds;
            se_index += _jl;
            //seds_index++;
            e_index += _jl;
        }
    }
}

void cpe_c_sbottom()
{
    int j, k, ks, se_index, sebo_index;
    float sbo, d0, wk0, dk, ssbo, tmpf;

    ks = thread_local_buff.ks0;
    d0 = cpe_ips_packed_data[D_OFFSET];
    sbo = 0.038 * 2.0 / _g;

    for (k = 1; k <= ks; k++)
    {
//        wk0 = cpe_ips_wk[k - 1];
//        dk  = d0 * wk0;
//        if (dk >= 30.0)
//        {
//            ssbo = 0.0;
//        } else {
//            ssbo = -_abo * sbo * wk0 / sinh(2.0 * dk);//pre tag
//        }

        ssbo = cpe_ips_ssbos[k - 1];

        for (j = 1; j <= _jl; j++)
        {
            se_index = (k - 1) * _jl + j - 1;
            sebo_index = _sexx_index(k, j);
            tmpf = ssbo * cpe_ips_e[j - 1 + (k - 1) * _jl];
            thread_local_buff.sebo[sebo_index] = tmpf;
            thread_local_buff.se[se_index]  += tmpf;
            thread_local_buff.dse[se_index] += ssbo;
        }
    }
}

void cpe_c_scurrent()
{
    int j, k, ks, se_index, e_index;
    floatv4 duxdx0, duxdy0, duydx0, duydy0, costh, sinth;
    floatv4 cost2, sint2, cost2_p1, sint2_p1, ws0, cgdc, cu1, cu2, cu3, sscu;

    floatv4 p5_v4, p5xx0, p5yy0, cgdc_acu, acu_v4 = _acu, onev4, macu_v4 = -_acu;
    floatv4 cgdc_dxx, cgdc_dyy, cgdc_dxyyx, esv4, zpi_v4, wf_v4;

    ks     = thread_local_buff.ks0;
    duxdx0 = cpe_ips_packed_data[UXX_OFFSET];
    duxdy0 = cpe_ips_packed_data[UXY_OFFSET];
    duydx0 = cpe_ips_packed_data[UYX_OFFSET];
    duydy0 = cpe_ips_packed_data[UYY_OFFSET];

    p5_v4 = 0.5;
    onev4 = 1.0;
    macu_v4 = -_acu;
    acu_v4 = _acu;
    zpi_v4 = _zpi;

    p5xx0 = p5_v4 * duxdx0;
    p5yy0 = p5_v4 * duydy0;

    for (k = 1; k <= ks; k++)
    {
        wf_v4 = cpe_ips_wf[k - 1];
        ws0 = zpi_v4 * wf_v4;
        cgdc = simd_set_floatv4(cpe_ips_wk[k - 1] * cpe_ips_ccg[k - 1],
                                cpe_ips_wk[k - 1] * cpe_ips_ccg[k - 1],
                                cpe_ips_wk[k - 1] * cpe_ips_ccg[k - 1],
                                cpe_ips_wk[k - 1] * cpe_ips_ccg[k - 1]);
        cgdc = cgdc / ws0;

        cgdc_acu = cgdc * acu_v4;
        cgdc_dxx = cgdc * duxdx0;
        cgdc_dyy = cgdc * duydy0;
        cgdc_dxyyx = duxdy0 + duydx0;
        cgdc_dxyyx = cgdc_dxyyx * cgdc;
        for (j = 1; j <= _jl; j += 4)
        {
            se_index = e_index  = j - 1 + (k - 1) * _jl;
            sinth  = simd_set_floatv4(cpe_ips_sin_theta[j - 1],cpe_ips_sin_theta[j],
                                        cpe_ips_sin_theta[j + 1],cpe_ips_sin_theta[j + 2]);

            costh  = simd_set_floatv4(cpe_ips_cos_theta[j - 1],cpe_ips_cos_theta[j],
                                        cpe_ips_cos_theta[j + 1],cpe_ips_cos_theta[j + 2]);

            esv4 = simd_set_floatv4(cpe_ips_e[e_index]
                               ,cpe_ips_e[e_index + 1]
                               ,cpe_ips_e[e_index + 2]
                               ,cpe_ips_e[e_index + 3]);

            cost2 = costh * costh;
            sint2 = sinth * sinth;
            cost2_p1 = cost2 + onev4;
            sint2_p1 = sint2 + onev4;

            cu1 = cgdc_dxx * cost2_p1 - p5xx0;
            cu2 = sinth * costh * cgdc_dxyyx;
            cu3 = cgdc_dyy * sint2_p1 - p5yy0;
            sscu = macu_v4 * (cu1 + cu2 + cu3);
            esv4 = esv4 * sscu;
            simd_store(sscu, &cpe_tmp_fv[0]);
            simd_store(esv4, &cpe_tmp_fv1[0]);

            thread_local_buff.dse[se_index] += cpe_tmp_fv[0];
            thread_local_buff.dse[se_index + 1] += cpe_tmp_fv[1];
            thread_local_buff.dse[se_index + 2] += cpe_tmp_fv[2];
            thread_local_buff.dse[se_index + 3] += cpe_tmp_fv[3];
            thread_local_buff.se[se_index] += cpe_tmp_fv1[0];
            thread_local_buff.se[se_index + 1] += cpe_tmp_fv1[1];
            thread_local_buff.se[se_index + 2] += cpe_tmp_fv1[2];
            thread_local_buff.se[se_index + 3] += cpe_tmp_fv1[3];

        }
    }

/*
    for (j = 1; j <= _jl; j++)
    {
        sinth = cpe_ips_sin_theta[j - 1];
        costh = cpe_ips_cos_theta[j - 1];
        cost2 = costh * costh;
        sint2 = sinth * sinth;
        se_index = _se_index(0, j);
        for (k = 1; k <= ks; k++)
        {
            se_index++;
            ws0 = _zpi * cpe_ips_wf[k - 1];
            cgdc = cpe_ips_wk[k - 1] * cpe_ips_ccg[k - 1] / ws0;
            cu1 = (cgdc * ( 1.0 + cost2) - 0.5) * duxdx0;
            cu2 = cgdc * sinth * costh * (duxdy0 + duydx0);
            cu3 = (cgdc * (1.0 + sint2) - 0.5) * duydy0;
            sscu = -_acu * (cu1 + cu2 + cu3);
            thread_local_buff.se[se_index]  += sscu * cpe_ips_e[(j - 1) * _kl + (k - 1)];
            thread_local_buff.dse[se_index] += sscu;
        }
    }*/
}

void cpe_c_update_ee_pe()
{
    int i, j, k, ks, jk_index, e_index, ee_index, se_index, es_index;
    floatv4 pein_tmp, peds_tmp,  pebo_tmp, wstar, gadiag, eef, eefab, sig, deltee, tmp;
    float cd;

    ks = thread_local_buff.ks0;
    cd = (0.8 + 0.065 * cpe_ips_packed_data[W_OFFSET]) * 0.001;
    cd = cpe_ips_packed_data[W_OFFSET] * sqrt(cd);
    wstar = cd;

    floatv4 eees, wkhi, tmpv4, ems, smv4 = _small, eev;
    floatv4 one_e19, one_v4, dses, ses, tmp_minus, grolim_v4, zero_v4;

    one_e19 = 1e-19;
    one_v4 = 1.0;
    zero_v4 = 0.0;

    for (k = 1; k <= ks; k++)
    {
        grolim_v4 = cpe_ips_grolim[k - 1];
        for (j = 1; j <= _jl; j += 4)
        {
            se_index = e_index = (k - 1) * _jl + (j - 1);
            ee_index = (j - 1) * _kl + (k - 1);

            dses = simd_set_floatv4(thread_local_buff.dse[se_index],
                                    thread_local_buff.dse[se_index + 1],
                                    thread_local_buff.dse[se_index + 2],
                                    thread_local_buff.dse[se_index + 3]);
            ses = simd_set_floatv4(thread_local_buff.se[se_index],
                                    thread_local_buff.se[se_index + 1],
                                    thread_local_buff.se[se_index + 2],
                                    thread_local_buff.se[se_index + 3]);
            gadiag = one_v4 - cpe_deltt5 * dses;
            tmp_minus = one_v4 - gadiag;
            gadiag = __builtin_sw64_sllt(tmp_minus, gadiag, one_v4);
            //gadiag = max(1.0, gadiag);

            eef = cpe_deltt * ses / gadiag;
            eefab = __builtin_sw64_vabss(eef);
            tmp_minus = one_e19 - eefab;
            eefab = __builtin_sw64_sllt(tmp_minus, eefab, one_e19);
            sig = eef / eefab;
            deltee = wstar * grolim_v4;
            tmp_minus = deltee - eefab;
            eefab = __builtin_sw64_sllt(tmp_minus, deltee, eefab);
            eees = simd_set_floatv4(cpe_ips_e[e_index]
                               ,cpe_ips_e[e_index + 1]
                               ,cpe_ips_e[e_index + 2]
                               ,cpe_ips_e[e_index + 3]);
            eef = eees + sig * eefab;
            tmp_minus = smv4 - eef;
            tmp = __builtin_sw64_sllt(tmp_minus, eef, smv4);
            simd_store(tmp, &cpe_tmp_fv[0]);
            cpe_ips_ee[ee_index] = cpe_tmp_fv[0];
            cpe_ips_ee[ee_index + 1 * _kl] = cpe_tmp_fv[1];
            cpe_ips_ee[ee_index + 2 * _kl] = cpe_tmp_fv[2];
            cpe_ips_ee[ee_index + 3 * _kl] = cpe_tmp_fv[3];
        }
    }


    for (k = ks + 1; k <= _kl; k++)
    {
        i = k - ks + 1;
        
        wkhi = cpe_ips_wkh[i - 1];

        for (j = 1; j <= _jl; j += 4)
        {
            e_index = k - 1 + (j - 1) * _kl;
            es_index = ks - 1 + (j - 1) * _kl;
            eees = simd_set_floatv4(cpe_ips_ee[es_index]
                               ,cpe_ips_ee[es_index + _kl]
                               ,cpe_ips_ee[es_index + 2 * _kl]
                               ,cpe_ips_ee[es_index + 3 * _kl]);

            tmpv4 = eees * wkhi;
            ems = smv4 - tmpv4;

            eev = __builtin_sw64_sllt(ems, tmpv4, smv4);
            simd_store(eev, &cpe_tmp_fv[0]);

            cpe_ips_ee[e_index] = cpe_tmp_fv[0];
            cpe_ips_ee[e_index + 1 * _kl] = cpe_tmp_fv[1];
            cpe_ips_ee[e_index + 2 * _kl] = cpe_tmp_fv[2];
            cpe_ips_ee[e_index + 3 * _kl] = cpe_tmp_fv[3];

            //cpe_ips_ee[e_index] = max(_small, tmp);
        }
    }

    pein_tmp = 0.0;
    pebo_tmp = 0.0;
    peds_tmp = 0.0;

    for (j = 1; j <= _jl; j++)
    {
        jk_index = (j - 1) * _kl + (0 - 1);
        for (k = 1; k <= _kl; k++)
        {
            jk_index++;
            pein_tmp += 2.0 * cpe_ips_dwk[k - 1] * thread_local_buff.sein[jk_index];
            pebo_tmp += 2.0 * cpe_ips_dwk[k - 1] * thread_local_buff.sebo[jk_index];
            peds_tmp += 2.0 * cpe_ips_dwk[k - 1] * thread_local_buff.seds[jk_index];
        }
    }

    cpe_pexx_output_buf[PEIN_OFFSET] = pein_tmp;
    cpe_pexx_output_buf[PEBO_OFFSET] = pebo_tmp;
    cpe_pexx_output_buf[PEDS_OFFSET] = peds_tmp;
}

void cpe_implsch_kernel()
{
    cpe_c_mean2();
    RPCC(edcc);
    out_flag[PF_MEAN2] += edcc - stcc;
    stcc = edcc;

    cpe_c_get_thresholds();

    RPCC(edcc);
    out_flag[PF_GETTH] += edcc - stcc;
    stcc = edcc;

    memset(thread_local_buff.se,  0, sizeof(float) * _jl * _klp1);
    memset(thread_local_buff.dse, 0, sizeof(float) * _jl * _klp1);

    cpe_c_snonlin ();

    RPCC(edcc);
    out_flag[PF_NOLIN] += edcc - stcc;
    stcc = edcc;

    cpe_c_sinput  ();

    RPCC(edcc);
    out_flag[PF_SINPUT] += edcc - stcc;
    stcc = edcc;

    cpe_c_sdissip ();

    RPCC(edcc);
    out_flag[PF_DISSIP] += edcc - stcc;
    stcc = edcc;

    cpe_c_sbottom ();

    RPCC(edcc);
    out_flag[PF_BOTTOM] += edcc - stcc;
    stcc = edcc;

    cpe_c_scurrent();

    RPCC(edcc);
    out_flag[PF_SCURRENT] += edcc - stcc;
    stcc = edcc;

    cpe_c_update_ee_pe();

    RPCC(edcc);
    out_flag[PF_UPDATEEE] += edcc - stcc;
    stcc = edcc;
}

void seq_write_cpe_mean1_em_data_to_mpe(int my_core, int compute_unit, int write_flag, float *mean1_out_buf)
{
    int wcid, em_offset;
    em_offset = (cpe_mean1_ic - cpe_iys) * cpe_dim4 + (cpe_mean1_ia - cpe_ixs) * cpe_dim3;
    for (wcid = 63; wcid >= 0; wcid--)
    {
        if (my_core == wcid && write_flag != 0)
        {
            //put_reply = 0;
            //athread_put(
            //    PE_MODE, &cpe_mean1_em[0], param.host_ee + em_offset,
            //    sizeof(float) * _jl * _kl, (void*)(&put_reply), 0, 0
            //);
            //while (put_reply != 1);

            put_reply = 0;
            athread_put(
                PE_MODE, &cpe_mean1_data_out[0], mean1_out_buf + compute_unit * GETEM_OUT_DATA,
                sizeof(float) * GETEM_OUT_DATA, (void*)(&put_reply), 0, 0
            );
            while (put_reply != 1);
        }
        athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
    }

    for (wcid = 2; wcid >= 0; wcid--)
    {
        if (my_core % 2 == wcid && write_flag != 0)
        {
            put_reply = 0;
            athread_put(
                PE_MODE, &cpe_mean1_em[0], param.host_ee + em_offset,
                sizeof(float) * _jl * _kl, (void*)(&put_reply), 0, 0
            );
            while (put_reply != 1);
        }
        athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
    }
}

void get_mean1_iaicnsp(int compute_unit, int* ia_set, int* ic_set, int* nsp_set, int* ia, int* ic, int* nsp, float* m1di, int* m1ii, float* m1pd, int* m1pi)
{
    int mean1_packed_data_offset = (compute_unit) * GETEM_IN_FLOAT;
    int mean1_packed_index_offset = (compute_unit) * GETEM_IN_INT;

    get_reply = 0;
    athread_get(
            PE_MODE, ia_set + compute_unit, ia,
            sizeof(int) * 1, (void*)(&get_reply), 0, 0, 0
        );
    while(get_reply != 1);
    get_reply = 0;
    athread_get(
            PE_MODE, ic_set + compute_unit, ic,
            sizeof(int) * 1, (void*)(&get_reply), 0, 0, 0
        );
    while(get_reply != 1);
    get_reply = 0;
    athread_get(
            PE_MODE, nsp_set + compute_unit, nsp,
            sizeof(int) * 1, (void*)(&get_reply), 0, 0, 0
        );
    while(get_reply != 1);

    get_reply = 0;
    athread_get(
        PE_MODE, m1pd + mean1_packed_data_offset, m1di,
        sizeof(float) * GETEM_IN_FLOAT, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);

    get_reply = 0;
    athread_get(
        PE_MODE, m1pi + mean1_packed_index_offset, m1ii,
        sizeof(int) * GETEM_IN_INT, (void*)(&get_reply), 0, 0, 0
    );
    while (get_reply != 1);
}

void cpe_propagat(int my_core)
{
    int ig, compute_unit, packed_data_offset,
            output_e_offset, loop_length, _4d_icia_index;

    int remain_length = local_flag[REMAIN_POINT];
    float *mpe_halo_ips_packed_data, *mpe_halo_ips_buffer;
    int *ppg_nsp_set_ptr;
    int my_start, my_len, next_ia, next_ic, ia_diff, ic_diff, other_start, cur_start;
    int* ia_set, *ic_set;
    float tmp_nsp[9];
    float* total_nsp_ptr;

    loop_length = local_flag[GROUP_SIZE];
    mpe_halo_ips_packed_data = (float*)(local_flag[PACKED_PTR]);
    mpe_halo_ips_buffer = (float*)(local_flag[OUT_DATA_PTR]);
    remain_length = local_flag[REMAIN_POINT];
    ia_set = (int*)(local_flag[IAS_PTR]);
    ic_set = (int*)(local_flag[ICS_PTR]);
    ppg_nsp_set_ptr = (int*)(local_flag[NSPS_PTR]);
    total_nsp_ptr = (float*)(local_flag[TOTAL_NSP_PTR]);

    //char str[100];

    int ee_offset, target_reply, i, j, k, tmp_ia, tmp_ic, now_offset;

    my_start = BLOCK_LOW(my_core, 64, loop_length * 64 + remain_length);
    my_len = BLOCK_SIZE(my_core, 64, loop_length * 64 + remain_length);
/*
    if(param.host_mpi_rank == 0)
    {
        for(ig = 0; ig < 64; ig++)
        {
            if(my_core == ig)
            {
                sprintf(&strlocal[0], "    ####    core == %d, my_start = \n", ppg_nsp_set_ptr);
                printf("%s", strlocal);
            }
            athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
        }
    }*/

    for(ig = 0; ig < loop_length + 1; ig++)
    {
        if(ig != 0)
        {
            standard_wait_flag(my_core);
        }

        if(my_len == ig)
        {
            two_step_write_data(param.host_e + output_e_offset, &cpe_ips_ee[0], sizeof(float) * _jl * _kl, my_core, 0);
            //seq_write_data(mpe_halo_ips_buffer + compute_unit * NUM_PEXX, &cpe_pexx_output_buf[0], sizeof(float) * NUM_PEXX, my_core, 0);

            standard_write_flag(my_core);
            break;
        }

        RPCC(stcc);

        ppgid = my_start + ig;

        if(ig == 0)
        {
            get_reply = 0;

            athread_get(
                        PE_MODE, ia_set + ppgid, &cpe_ppg_ia,
                        sizeof(int) * 1, (void*)(&get_reply), 0, 0, 0
                        );
            athread_get(
                        PE_MODE, ic_set + ppgid, &cpe_ppg_ic,
                        sizeof(int) * 1, (void*)(&get_reply), 0, 0, 0
                        );

            while(get_reply != 2);
        }
        else
        {
            cpe_ppg_ia = next_ia;
            cpe_ppg_ic = next_ic;
        }

        if(ig < my_len - 1)//eg, ll = 1, then when i = 0, no need.
        {
            get_reply = 0;

            athread_get(
                        PE_MODE, ia_set + (ppgid + 1), &next_ia,
                        sizeof(int) * 1, (void*)(&get_reply), 0, 0, 0
                        );
            athread_get(
                        PE_MODE, ic_set + (ppgid + 1), &next_ic,
                        sizeof(int) * 1, (void*)(&get_reply), 0, 0, 0
                        );

            while(get_reply != 2);


            target_reply = 0;
            for(i = -1; i <= 1; i++)
            {
                for(j = -1; j <=1; j++)
                {
                    tmp_ic = next_ic + i;
                    tmp_ia = next_ia + j;

                    ee_offset = (tmp_ic - cpe_iys) * cpe_ix_size + (tmp_ia - cpe_ixs);

                    get_reply = 0;
                    athread_get(
                        PE_MODE, total_nsp_ptr + ee_offset, &tmp_nsp[target_reply],
                    sizeof(float) * 1, (void*)(&get_reply), 0, 0, 0
                    );
                    while(get_reply != 1);
                    target_reply++;
                }
            }
        }


        get_reply = 0;
        athread_get(
                    PE_MODE, ppg_nsp_set_ptr + ppgid, &cpe_ppg_nsp,
                    sizeof(int) * 1, (void*)(&get_reply), 0, 0, 0
                    );
        while(get_reply != 1);

        _4d_icia_index = (cpe_ppg_ic - cpe_iys) * cpe_dim4 + (cpe_ppg_ia - cpe_ixs) * cpe_dim3;

        //if(my_core == 63 && param.host_mpi_rank == 0)
        //    printf("    ####    %d, %d, %d, %d, %d, %d\n", my_core, ig, ppgid, cpe_ppg_ia, cpe_ppg_ic, _4d_icia_index);

        if(ig == 0)//no prefetch.
        {
            get_reply = 0;
            athread_get(
                        PE_MODE, param.host_ppg_values + 6 * _4d_icia_index, &cpe_ppg_values[0],
                    sizeof(float) * 6 * _kl * _jl, (void*)(&get_reply), 0, 0, 0
                    );
            while(get_reply != 1);
            get_reply = 0;
            athread_get(
                        PE_MODE, param.host_ppg_idxs + 5 * _4d_icia_index, &cpe_ppg_idxs[0],
                    sizeof(int) * 5 * _kl * _jl, (void*)(&get_reply), 0, 0, 0
                    );
            while(get_reply != 1);

            target_reply = 0;
            get_reply = 0;
            now_offset = 0;
            for(i = -1; i <= 1; i++)
            {
                for(j = -1; j <=1; j++)
                {
                    tmp_ic = cpe_ppg_ic + i;
                    tmp_ia = cpe_ppg_ia + j;
                    //if(tmp_ia >= cpe_ixs && tmp_ia <= cpe_ixl && tmp_ic >= cpe_iys && tmp_ic <= cpe_iyl)
                    //{
                    target_reply++;
                    ee_offset = (tmp_ic - cpe_iys) * cpe_dim4 + (tmp_ia - cpe_ixs) * cpe_dim3;
                    athread_get(
                                PE_MODE, param.host_ee + ee_offset, &cpe_ppg_ee_now[now_offset],
                                sizeof(float) * _kl * _jl, (void*)(&get_reply), 0, 0, 0
                                );
                    //}

                    now_offset += _kl * _jl;
                }
            }
            while(get_reply < target_reply);
        }

        RPCC(edcc);
        out_flag[PF_PPG_DATA] += edcc - stcc;
        stcc = edcc;

        //output_e_offset    = (ig * 64 + my_core) * _jl * _kl;
        output_e_offset = _4d_icia_index;
        //packed_data_offset = output_e_offset * PPG_BI_DATA;

        // ???? propagat
        cpe_propagat_kernel(cpe_ppg_ia, cpe_ppg_ic, my_core);

        RPCC(edcc);
        out_flag[PF_PPG_KERNEL] += edcc - stcc;
        stcc = edcc;

        cpe_ips_ia = cpe_ppg_ia;
        cpe_ips_ic = cpe_ppg_ic;
        compute_unit = ppgid;

        // ?????????????????????
        load_mpe_ips_packed_data(compute_unit, mpe_halo_ips_packed_data);
        load_mpe_ips_e_wf_ccg();

        for(ee_offset = 0; ee_offset < _jl * _kl; ee_offset++)
            cpe_ips_e[ee_offset] = cpe_ppg_e[ee_offset];

        RPCC(edcc);
        out_flag[PF_IPS_DATA] += edcc - stcc;
        stcc = edcc;

        target_reply = 0;

        if(ig < my_len - 1)//prefetch
        {
            ppg_copy_count = 0;
            memset(&ppg_input_mask[0], 0, sizeof(int) * 9);

            //check this and the next

            if(abs(cpe_ppg_ia - next_ia) < 3 && abs(cpe_ppg_ic - next_ic))
            {
                //two cases: right and below move
                //ignore left and above
                ia_diff = next_ia - cpe_ppg_ia;
                ic_diff = next_ic - cpe_ppg_ic;
                other_start = 0;
                cur_start = 0;
                if(cpe_ppg_ia < next_ia)//right. move to the left. ia_diff > 0
                {
                    if(ic_diff < 0)
                        other_start = -ic_diff;
                    else
                        cur_start = ic_diff;
                    for(i = 0; i < 3 - ia_diff; i++)
                    {
                        for(j = 0; j < 3 - abs(ic_diff); j++)
                        {
                            ppg_dst_x[ppg_copy_count] = i;
                            ppg_dst_y[ppg_copy_count] = j + other_start;
                            ppg_src_x[ppg_copy_count] = i + ia_diff;
                            ppg_src_y[ppg_copy_count] = j + cur_start;
                            ppg_copy_count++;
                        }
                    }
                }
                else if(cpe_ppg_ic < next_ic)//move to the top. ic_diff > 0
                {
                    if(ia_diff < 0)
                        other_start = -ia_diff;
                    else// only 0 :)
                        cur_start = ia_diff;
                    for(i = 0; i < 3 - ic_diff; i++)
                    {
                        for(j = 0; j < 3 - abs(ia_diff); j++)
                        {
                            ppg_dst_y[ppg_copy_count] = i;
                            ppg_dst_x[ppg_copy_count] = j + other_start;
                            ppg_src_y[ppg_copy_count] = i + ic_diff;
                            ppg_src_x[ppg_copy_count] = j + cur_start;
                            ppg_copy_count++;
                        }
                    }
                }
            }

            for(i = 0; i < ppg_copy_count; i++)
            {
                memcpy(&cpe_ppg_ee_now[(ppg_dst_x[i] + ppg_dst_y[i] * 3) * _klmjl], &cpe_ppg_ee_now[(ppg_src_x[i] + ppg_src_y[i] * 3) * _klmjl], sizeof(float) * _klmjl);
                ppg_input_mask[ppg_dst_x[i] + ppg_dst_y[i] * 3] = 1;
            }

            _4d_icia_index = (next_ic - cpe_iys) * cpe_dim4 + (next_ia - cpe_ixs) * cpe_dim3;
            get_reply = 0;

            athread_get(
                        PE_MODE, param.host_ppg_values + 6 * _4d_icia_index, &cpe_ppg_values[0],
                    sizeof(float) * 6 * _kl * _jl, (void*)(&get_reply), 0, 0, 0
                    );
            //while(get_reply != 1);
            //get_reply = 0;
            athread_get(
                        PE_MODE, param.host_ppg_idxs + 5 * _4d_icia_index, &cpe_ppg_idxs[0],
                    sizeof(int) * 5 * _kl * _jl, (void*)(&get_reply), 0, 0, 0
                    );
            //while(get_reply != 1);

            target_reply = 2;
            //get_reply = 0;
            now_offset = 0;
            for(i = -1; i <= 1; i++)
            {
                for(j = -1; j <=1; j++)
                {
                    tmp_ic = next_ic + i;
                    tmp_ia = next_ia + j;
                    //if(tmp_ia >= cpe_ixs && tmp_ia <= cpe_ixl && tmp_ic >= cpe_iys && tmp_ic <= cpe_iyl)
                    //{
                    if(ppg_input_mask[(j + 1) + (i + 1) * 3] == 1)
                    {
                        now_offset += _kl * _jl;
                        continue;
                    }
                    //{
                    if(tmp_nsp[(j + 1) + (i + 1) * 3] == 0.0)
                    {
                        memcpy(&cpe_ppg_ee_now[now_offset], &cpe_nomarine_es[0], sizeof(float) * _jl * _kl);
                        now_offset += _kl * _jl;
                        continue;
                    }

                    target_reply++;
                    ee_offset = (tmp_ic - cpe_iys) * cpe_dim4 + (tmp_ia - cpe_ixs) * cpe_dim3;
                    athread_get(
                                PE_MODE, param.host_ee + ee_offset, &cpe_ppg_ee_now[now_offset],
                                sizeof(float) * _kl * _jl, (void*)(&get_reply), 0, 0, 0
                                );

                    now_offset += _kl * _jl;
                }
            }
        }

        RPCC(edcc);
        out_flag[PF_PPG_DATA] += edcc - stcc;
        stcc = edcc;

        cpe_implsch_kernel();

        while(get_reply < target_reply);

        RPCC(edcc);
        out_flag[PF_PPG_DATA] += edcc - stcc;
        stcc = edcc;


        // ????????д????????
        two_step_write_data(param.host_e + output_e_offset, &cpe_ips_ee[0], sizeof(float) * _jl * _kl, my_core, 1);
        //seq_write_data(mpe_halo_ips_buffer + compute_unit * NUM_PEXX, &cpe_pexx_output_buf[0], sizeof(float) * NUM_PEXX, my_core, 1);

        // ????д?????
        standard_write_flag(my_core);

        RPCC(edcc);
        out_flag[PF_IPS_WRITE] += edcc - stcc;
        stcc = edcc;
    }
}

void cpe_mean1(int my_core)//??????halo_ips????
{
    int ig, compute_unit, loop_length, remain_length;
    int my_max_ig, next_ia, next_ic;
    volatile int far_get_flag;
    float *mean1_packed_data, *mean1_buffer;
    int *ia_set, *ic_set;
    float *next_ee_in = &cpe_ppg_ee_now[0];
    int *nsp_set_ptr;
    int *mean1_packed_index;
    int mean1_packed_data_offset, mean1_packed_index_offset;
    int e_offset[5];
    int wf_dwf_offset;

    int k, j, kj, k1, i1, i;
    float n, a, d;
    int index_tmp;
    //float factor;

    float aets, aetc, thmax, akmax, eemax, eformax;
    floatv4 dwkk, wfk, wfk1, wsk, wsk1, wkk, wkk1, eef0, eekj, eekj1, eekjth;
    floatv4 sinth, costh; 
    float aett, chbh;
    float aet;

    float asi, awk, ae, ark, awf, hb, hbb, ape, tpf, h1_3;
    floatv4 asiv4, awkv4, aev4, awfv4, apev4, aetsv4, aetcv4;
    float a_temp[4];
    floatv4 deltth={_zpi/_jl, _zpi/_jl, _zpi/_jl, _zpi/_jl};

    floatv4 *fv4p1, *fv4p2, *fv4p3, *fv4p4, *fv4p5;
    const floatv4 afv4 = {24.0, 24.0, 24.0, 24.0};
    floatv4 nfv4; float* nfvp = &nfv4;

    get_reply_god = 0;
    get_reply_target = 0;

    floatv4 tmp_fv4_11, tmp_fv4_12, tmp_fv4_13, tmp_fv4_14, tmp_fv4_15, tmp_fv4_16;
    floatv4 tmp_fv4_21, tmp_fv4_22, tmp_fv4_23, tmp_fv4_24, tmp_fv4_25, tmp_fv4_26;
    floatv4 tmp_fv4_31, tmp_fv4_32, tmp_fv4_33, tmp_fv4_34, tmp_fv4_35, tmp_fv4_36;
    floatv4 tmp_fv4_41, tmp_fv4_42, tmp_fv4_43, tmp_fv4_44, tmp_fv4_45, tmp_fv4_46;
    floatv4 tmp_fv4_51, tmp_fv4_52, tmp_fv4_53, tmp_fv4_54, tmp_fv4_55, tmp_fv4_56;
    floatv4 efv4_11, efv4_12, efv4_13, efv4_14, efv4_15;
    floatv4 efv4_21, efv4_22, efv4_23, efv4_24, efv4_25;
    floatv4 efv4_31, efv4_32, efv4_33, efv4_34, efv4_35;
    floatv4 efv4_41, efv4_42, efv4_43, efv4_44, efv4_45;
    floatv4 efv4_51, efv4_52, efv4_53, efv4_54, efv4_55;

    loop_length = local_flag[GROUP_SIZE];
    mean1_packed_data = (float*)(local_flag[PACKED_PTR]);
    mean1_buffer = (float*)(local_flag[OUT_DATA_PTR]);
    mean1_packed_index = (int*)(local_flag[INDEXS_PTR]);
    remain_length = local_flag[REMAIN_POINT];
    ia_set = (int*)(local_flag[IAS_PTR]);
    ic_set = (int*)(local_flag[ICS_PTR]);
    nsp_set_ptr = (int*)(local_flag[NSPS_PTR]);

    my_max_ig = loop_length;
    if(remain_length <= my_core)
        my_max_ig = my_max_ig - 1;

    int my_start = BLOCK_LOW(my_core, 64, loop_length * 64 + remain_length);
    int my_len = BLOCK_SIZE(my_core, 64, loop_length * 64 + remain_length);

    cpe_mean1_ees = &cpe_ppg_ee_now[0];
    cpe_mean1_ees_next = &cpe_ppg_ee_now[_jl * _kl * 5];
    float* fp_tmp;

    double tmpd1, tmpd2, tmpd3, tmpd4;

    for(ig = 0; ig < loop_length + 1; ig++)
    {
        fp_tmp = cpe_mean1_ees;
        cpe_mean1_ees = cpe_mean1_ees_next;
        cpe_mean1_ees_next = fp_tmp;
        if(ig != 0)
        {
            standard_wait_flag(my_core);
        }

        if(my_len == ig)
        {
            seq_write_cpe_mean1_em_data_to_mpe(my_core, compute_unit, 0, mean1_buffer);
            standard_write_flag(my_core);
            break;
        }

        RPCC(stcc);

        compute_unit = my_start + ig;

        if(ig == 0)
        {
            get_mean1_iaicnsp(compute_unit, ia_set, ic_set, nsp_set_ptr, &cpe_mean1_ia, &cpe_mean1_ic, &cpe_mean1_nsp, &cpe_mean1_data_in[0], &cpe_mean1_index_in[0], &mean1_packed_data[0], &mean1_packed_index[0]);

            wf_dwf_offset = (cpe_mean1_ic - cpe_iys) * (_kldp1 * cpe_ix_size) + (cpe_mean1_ia - cpe_ixs) * _kldp1;
            get_reply = 0;
            athread_get(
                PE_MODE, param.host_wf + wf_dwf_offset, &cpe_mean1_wf[0],
                sizeof(float) * _kldp1, (void*)(&get_reply), 0, 0, 0
            );
            while (get_reply != 1);

            get_reply = 0;
            athread_get(
                PE_MODE, param.host_dwf + wf_dwf_offset, &cpe_mean1_dwf[0],
                sizeof(float) * _kldp1, (void*)(&get_reply), 0, 0, 0
            );
            while (get_reply != 1);

            RPCC(edcc);
            out_flag[PF_MEAN1_DATA] += edcc - stcc;
            stcc = edcc;

            e_offset[0] = (cpe_mean1_ic - cpe_iys) * cpe_dim4 + (cpe_mean1_ia - cpe_ixs) * cpe_dim3;
            e_offset[1] = e_offset[0] + cpe_mean1_index_in[LEFT_OFFSET];
            e_offset[2] = e_offset[0] + cpe_mean1_index_in[RIGHT_OFFSET];
            e_offset[3] = e_offset[0] + cpe_mean1_index_in[ABOVE_OFFSET];
            e_offset[4] = e_offset[0] + cpe_mean1_index_in[BELOW_OFFSET];

            for(i = 0; i < 5; i++)
            {
                get_reply = 0;
                athread_get(
                    PE_MODE, param.host_e + e_offset[i], &cpe_mean1_ees[i * _jl * _kl],
                    sizeof(float) * _kl * _jl, (void*)(&get_reply), 0, 0, 0
                );
                while (get_reply != 1);
            }

            RPCC(edcc);
            out_flag[PF_MEAN1_E] += edcc - stcc;
            stcc = edcc;
        }
        else
        {
            cpe_mean1_ia = cpe_mean1_next_ia;
            cpe_mean1_ic = cpe_mean1_next_ic;
            cpe_mean1_nsp = cpe_mean1_next_nsp;

            memcpy(&cpe_mean1_data_in[0], &cpe_mean1_data_next[0], GETEM_IN_FLOAT * sizeof(float));
            memcpy(&cpe_mean1_index_in[0], &cpe_mean1_index_next[0], GETEM_IN_INT * sizeof(int));
            memcpy(&cpe_mean1_wf[0], &cpe_mean1_next_wf[0], sizeof(float) * _kldp1);
            memcpy(&cpe_mean1_dwf[0], &cpe_mean1_next_dwf[0], sizeof(float) * _kldp1);
            //memcpy(&cpe_mean1_ees[0], &cpe_mean1_ees_next[0], sizeof(float) * _kl * _jl * 5);

            RPCC(edcc);
            out_flag[PF_MEAN1_DATA] += edcc - stcc;
            stcc = edcc;
        }

        if(ig < my_len - 1)
        {
            get_mean1_iaicnsp(compute_unit + 1, ia_set, ic_set, nsp_set_ptr, &cpe_mean1_next_ia, &cpe_mean1_next_ic, &cpe_mean1_next_nsp, &cpe_mean1_data_next[0], &cpe_mean1_index_next[0], &mean1_packed_data[0], &mean1_packed_index[0]);

            wf_dwf_offset = (cpe_mean1_next_ic - cpe_iys) * (_kldp1 * cpe_ix_size) + (cpe_mean1_next_ia - cpe_ixs) * _kldp1;
            
            athread_get(
                PE_MODE, param.host_wf + wf_dwf_offset, &cpe_mean1_next_wf[0],
                sizeof(float) * _kldp1, (void*)(&get_reply_god), 0, 0, 0
            );
            get_reply_target++;
            //while (get_reply_god != 1);

            //get_reply = 0;
            athread_get(
                PE_MODE, param.host_dwf + wf_dwf_offset, &cpe_mean1_next_dwf[0],
                sizeof(float) * _kldp1, (void*)(&get_reply_god), 0, 0, 0
            );
            get_reply_target++;
            //while (get_reply_god != 1);


            e_offset[0] = (cpe_mean1_next_ic - cpe_iys) * cpe_dim4 + (cpe_mean1_next_ia - cpe_ixs) * cpe_dim3;
            e_offset[1] = e_offset[0] + cpe_mean1_index_next[LEFT_OFFSET];
            e_offset[2] = e_offset[0] + cpe_mean1_index_next[RIGHT_OFFSET];
            e_offset[3] = e_offset[0] + cpe_mean1_index_next[ABOVE_OFFSET];
            e_offset[4] = e_offset[0] + cpe_mean1_index_next[BELOW_OFFSET];

            for(i = 0; i < 5; i++)
            {
                //get_reply = 0;
                athread_get(
                    PE_MODE, param.host_e + e_offset[i], &cpe_mean1_ees_next[i * _jl * _kl],
                    sizeof(float) * _kl * _jl, (void*)(&get_reply_god), 0, 0, 0
                );
                get_reply_target++;
                //while (get_reply != 1);
            }

            RPCC(edcc);
            out_flag[PF_MEAN1_DATA2] += edcc - stcc;
            stcc = edcc;
        }

        //if(param.host_mpi_rank == 0 && my_core == 0)
        //    printf("    ####    load mean1 data finished\n");

        //cpe_mean1_kernel(my_core);

        ae = 0.0;     aev4 = 0.0;
        asi = 0.0;    asiv4 = 0.0;
        awf = 0.0;    apev4 = 0.0;
        awk = 0.0;    awfv4 = 0.0;
        ark = 0.0;    awkv4 = 0.0;
        hb = 0.0;
        hbb = 0.0;

        ape = 0.0;
        n = cpe_mean1_data_in[EM_N_OFFSET];
        d = cpe_mean1_data_in[EM_D_OFFSET];

        nfvp[0] = n; nfvp[1] = n; nfvp[2] = n; nfvp[3] = n;
        asm volatile ("#nop":::"memory");

        //if(param.host_mpi_rank == 0 && my_core == 0 && cpe_mean1_ia == 2 && cpe_mean1_ic == 2)
        //    printf("    ####    offsets: %d, %d, %d, %d, %d,   %d, %d, %d, %d  n = %.8f, d = %.8f, index offset = %d\n",
        //           e_offset[0], e_offset[1], e_offset[2], e_offset[3], e_offset[4],
        //            cpe_mean1_index_in[LEFT_OFFSET], cpe_mean1_index_in[RIGHT_OFFSET],
        //            cpe_mean1_index_in[ABOVE_OFFSET], cpe_mean1_index_in[BELOW_OFFSET],
        //            n, d, mean1_packed_index_offset);

        ae = 0.0;     aev4 = 0.0;
        asi = 0.0;    asiv4 = 0.0;
        awf = 0.0;    apev4 = 0.0;
        awk = 0.0;    awfv4 = 0.0;
        ark = 0.0;    awkv4 = 0.0;

        //factor = 1.0 / n;
        a = 24;
/*
        for(kj = 0; kj < _jl * _kl; kj++)
        {
            cpe_mean1_em[kj] = (a * cpe_mean1_ees[kj] +
                                cpe_mean1_ees[_jl * _kl * 1 + kj] +
                                cpe_mean1_ees[_jl * _kl * 2 + kj] +
                                cpe_mean1_ees[_jl * _kl * 3 + kj] +
                                cpe_mean1_ees[_jl * _kl * 4 + kj]) / n;
        }*/

        //****
        memcpy(&cpe_mean1_em[0], cpe_mean1_ees + 300, sizeof(float) * 300);

        for(kj = 0; kj < 300; kj += 20)//dirty
        {
            efv4_11 = __builtin_sw64_loadu (efv4_11, cpe_mean1_ees + kj);
            efv4_12 = __builtin_sw64_loadu (efv4_12, cpe_mean1_ees + 300 + kj);
            efv4_13 = __builtin_sw64_loadu (efv4_13, cpe_mean1_ees + 600 + kj);
            efv4_14 = __builtin_sw64_loadu (efv4_14, cpe_mean1_ees + 900 + kj);
            efv4_15 = __builtin_sw64_loadu (efv4_15, cpe_mean1_ees + 1200 + kj);
            efv4_21 = __builtin_sw64_loadu (efv4_21, cpe_mean1_ees + 4    + kj);
            efv4_22 = __builtin_sw64_loadu (efv4_22, cpe_mean1_ees + 304  + kj);
            efv4_23 = __builtin_sw64_loadu (efv4_23, cpe_mean1_ees + 604  + kj);
            efv4_24 = __builtin_sw64_loadu (efv4_24, cpe_mean1_ees + 904  + kj);
            efv4_25 = __builtin_sw64_loadu (efv4_25, cpe_mean1_ees + 1204 + kj);
            efv4_31 = __builtin_sw64_loadu (efv4_31, cpe_mean1_ees + 8    + kj);
            efv4_32 = __builtin_sw64_loadu (efv4_32, cpe_mean1_ees + 308  + kj);
            efv4_33 = __builtin_sw64_loadu (efv4_33, cpe_mean1_ees + 608  + kj);
            efv4_34 = __builtin_sw64_loadu (efv4_34, cpe_mean1_ees + 908  + kj);
            efv4_35 = __builtin_sw64_loadu (efv4_35, cpe_mean1_ees + 1208 + kj);
            efv4_41 = __builtin_sw64_loadu (efv4_41, cpe_mean1_ees + 12   + kj);
            efv4_42 = __builtin_sw64_loadu (efv4_42, cpe_mean1_ees + 312  + kj);
            efv4_43 = __builtin_sw64_loadu (efv4_43, cpe_mean1_ees + 612  + kj);
            efv4_44 = __builtin_sw64_loadu (efv4_44, cpe_mean1_ees + 912  + kj);
            efv4_45 = __builtin_sw64_loadu (efv4_45, cpe_mean1_ees + 1212 + kj);
            efv4_51 = __builtin_sw64_loadu (efv4_51, cpe_mean1_ees + 16   + kj);
            efv4_52 = __builtin_sw64_loadu (efv4_52, cpe_mean1_ees + 316  + kj);
            efv4_53 = __builtin_sw64_loadu (efv4_53, cpe_mean1_ees + 616  + kj);
            efv4_54 = __builtin_sw64_loadu (efv4_54, cpe_mean1_ees + 916  + kj);
            efv4_55 = __builtin_sw64_loadu (efv4_55, cpe_mean1_ees + 1216 + kj);

            tmp_fv4_11 = efv4_11 * afv4 + efv4_12;
            tmp_fv4_21 = efv4_21 * afv4 + efv4_22;
            tmp_fv4_31 = efv4_31 * afv4 + efv4_32;
            tmp_fv4_41 = efv4_41 * afv4 + efv4_42;
            tmp_fv4_51 = efv4_51 * afv4 + efv4_52;
            tmp_fv4_12 = efv4_13 + efv4_14;
            tmp_fv4_22 = efv4_23 + efv4_24;
            tmp_fv4_32 = efv4_33 + efv4_34;
            tmp_fv4_42 = efv4_43 + efv4_44;
            tmp_fv4_52 = efv4_53 + efv4_54;
            tmp_fv4_11 = tmp_fv4_11 + efv4_15;
            tmp_fv4_21 = tmp_fv4_21 + efv4_25;
            tmp_fv4_31 = tmp_fv4_31 + efv4_35;
            tmp_fv4_41 = tmp_fv4_41 + efv4_45;
            tmp_fv4_51 = tmp_fv4_51 + efv4_55;
            tmp_fv4_12 = (tmp_fv4_12 + tmp_fv4_11) / nfv4;
            tmp_fv4_22 = (tmp_fv4_22 + tmp_fv4_21) / nfv4;
            tmp_fv4_32 = (tmp_fv4_32 + tmp_fv4_31) / nfv4;
            tmp_fv4_42 = (tmp_fv4_42 + tmp_fv4_41) / nfv4;
            tmp_fv4_52 = (tmp_fv4_52 + tmp_fv4_51) / nfv4;

            simd_store(tmp_fv4_12, &cpe_mean1_em[kj + 0 ]);
            simd_store(tmp_fv4_22, &cpe_mean1_em[kj + 4 ]);
            simd_store(tmp_fv4_32, &cpe_mean1_em[kj + 8 ]);
            simd_store(tmp_fv4_42, &cpe_mean1_em[kj + 12]);
            simd_store(tmp_fv4_52, &cpe_mean1_em[kj + 16]);
        }
        //****

        RPCC(edcc);
        out_flag[PF_EM] += edcc - stcc;
        stcc = edcc;

        #pragma frequency_hint frequent
        if(cpe_mean1_nsp == 1)
        {

            //if(param.host_mpi_rank == 0 && my_core == 0)
            //    printf("    ####    em calculation finished\n");

            thmax=0.0;
            akmax=0.0;
            eemax=-999;
            eformax=-999;

            float zpiv4=_zpi;
            for(k = 1; k <= _kl; k++)
            {
                k1=k+1;
                i=k-_kl+1;
                i1=i+1;

                dwkk = cpe_ips_dwk[k - 1];
                //index_tmp = (ic - _iys) * (_kldp1 * ix_size) + (ia - _ixs) * _kldp1;
                //dwfk  = cpe_mean1_dwf[k - 1];
                wfk  = cpe_mean1_wf[k - 1];
                wfk1 = cpe_mean1_wf[k1 - 1];

                wsk  = zpiv4 * wfk;
                wsk1 = zpiv4 * wfk1;
                wkk  = cpe_ips_wk[k - 1];
                wkk1 = cpe_ips_wk[k1 - 1];
                eef0 = 0.0;
                //if(param.host_mpi_rank == 0 && my_core == 0 && cpe_mean1_ia == 2 && cpe_mean1_ic == 2)
                //    printf("\n    ####    k = %d,  dwkk = %.8f,  dwfk = %.8f,  wfk = %.8f,  wfk1 = %.8f,  wkk = %.8f,  wkk1 = %.8f\n", k, dwkk, dwfk, wfk, wfk1, wkk, wkk1);


                for (j = 1; j <= _jl; j+=4)
                {
                    int _index_tmp[4] = {(j - 1) * _kl, j * _kl, (j + 1) * _kl, (j + 2) * _kl};

                    if (k < _kl)
                    {
                        eekj  = simd_set_floatv4(cpe_mean1_em[_index_tmp[0] + k - 1]* 1.0
                                                 ,cpe_mean1_em[_index_tmp[1] + k - 1]* 1.0
                                                 ,cpe_mean1_em[_index_tmp[2] + k - 1]* 1.0
                                                 ,cpe_mean1_em[_index_tmp[3] + k - 1] * 1.0);
                        eekj1 = simd_set_floatv4(cpe_mean1_em[_index_tmp[0] + k1 - 1]* 1.0
                                                 ,cpe_mean1_em[_index_tmp[1] + k1 - 1]* 1.0
                                                 ,cpe_mean1_em[_index_tmp[2] + k1 - 1]* 1.0
                                                 ,cpe_mean1_em[_index_tmp[3] + k1 - 1]* 1.0);
                    } else {
                        eekj  = simd_set_floatv4(cpe_mean1_em[_index_tmp[0] + _kl - 1] * cpe_ips_wkh[i - 1]
                                                ,cpe_mean1_em[_index_tmp[1] + _kl - 1] * cpe_ips_wkh[i - 1]
                                                ,cpe_mean1_em[_index_tmp[2] + _kl - 1] * cpe_ips_wkh[i - 1]
                                                ,cpe_mean1_em[_index_tmp[3] + _kl - 1] * cpe_ips_wkh[i - 1]);
                        eekj1 = simd_set_floatv4(cpe_mean1_em[_index_tmp[0] + _kl - 1] * cpe_ips_wkh[i1 - 1]
                                                ,cpe_mean1_em[_index_tmp[1] + _kl - 1] * cpe_ips_wkh[i1 - 1]
                                                ,cpe_mean1_em[_index_tmp[2] + _kl - 1] * cpe_ips_wkh[i1 - 1]
                                                ,cpe_mean1_em[_index_tmp[3] + _kl - 1] * cpe_ips_wkh[i1 - 1]);//implement from mean2

                    }

                    eekjth=eekj*deltth;
                    eef0=eef0+eekjth;
                    sinth = simd_set_floatv4(cpe_ips_sin_theta[j - 1]
                                            ,cpe_ips_sin_theta[ j ]
                                            ,cpe_ips_sin_theta[j + 1]
                                            ,cpe_ips_sin_theta[j + 2]);
                    costh = simd_set_floatv4(cpe_ips_cos_theta[j - 1]
                                            ,cpe_ips_cos_theta[ j ]
                                            ,cpe_ips_cos_theta[j + 1]
                                            ,cpe_ips_cos_theta[j + 2]);

                    //if(my_rank == 0)
                    //    printf("    @@@@    AAAAAA\n");
                    
                    
                    aev4 = aev4 + (eekj+eekj1)*dwkk;
                    asiv4 = asiv4 + (eekj/(wfk*wfk)+eekj1/(wfk1*wfk1))*dwkk;
                    apev4 = apev4 + (eekj*(wsk*wsk)+eekj1*(wsk1*wsk1))*dwkk;
                    awfv4 = awfv4 + (eekj*wfk+eekj1*wfk1)*dwkk;
                    awkv4 = awkv4 + (eekj*wkk+eekj1*wkk1)*dwkk;
                    floatv4 cpe_ips_wkv4=cpe_ips_wk[k - 1];
                    aetsv4 = aetsv4 + (eekj+eekj1)*cpe_ips_wkv4*sinth*dwkk;
                    aetcv4 = aetcv4 + (eekj+eekj1)*cpe_ips_wkv4*costh*dwkk;
            

                    //if(my_rank == 0)
                    //    printf("    @@@@    AAAAAB\n");
                    //if(param.host_mpi_rank == 0 && my_core == 0 && cpe_mean1_ia == 2 && cpe_mean1_ic == 2)
                    //    printf("    ####    j = %d,  %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f\n",
                    //           j, eekj, eekj1, eekjth, eef0, sinth, costh, ae, asi, cpe_mean1_data_out[APE_OFFSET], awf, awk, aets, aetc);
                }
            }
            simd_store(aev4,&a_temp[0]); 
            tmpd1 = ae;
            tmpd1 += a_temp[0];
            tmpd2 = a_temp[1];
            tmpd3 = a_temp[2];
            tmpd4 = a_temp[3];
            tmpd1 = tmpd1 + tmpd2 + tmpd3 + tmpd4;
            //ae=ae+a_temp[0]+a_temp[1]+a_temp[2]+a_temp[3];
            ae = tmpd1;
                    
            simd_store(asiv4,&a_temp[0]);
            asi=asi+a_temp[0]+a_temp[1]+a_temp[2]+a_temp[3];
                    
            simd_store(apev4,&a_temp[0]);
            ape=ape+a_temp[0]+a_temp[1]+a_temp[2]+a_temp[3];
                    
            simd_store(awfv4,&a_temp[0]);
            awf=awf+a_temp[0]+a_temp[1]+a_temp[2]+a_temp[3];
                    
            simd_store(awkv4,&a_temp[0]);
            awk=awk+a_temp[0]+a_temp[1]+a_temp[2]+a_temp[3];
                    
            simd_store(aetsv4,&a_temp[0]);
            aets=aets+a_temp[0]+a_temp[1]+a_temp[2]+a_temp[3];
                    
            simd_store(aetcv4,&a_temp[0]);
            aetc=aetc+a_temp[0]+a_temp[1]+a_temp[2]+a_temp[3];


            RPCC(edcc);
            out_flag[PF_MEAN1_1] += edcc - stcc;
            stcc = edcc;

            ape=_tztz*_zpi/sqrt(ape/ae);
            tpf=(asi/ae)*(awf/ae);


            if(abs(aetc) <= 0.000001)
                aetc=0.00001;

            aett=atan2(aets,aetc);

            if (aett<=0.) aett=360.+aett;
            aet=aett;
            awk=awk/ae;
            tmpd1 = ae;
            tmpd2 = sqrt(tmpd1);
            tmpd3 = 4.0 * tmpd2;
            h1_3 = tmpd3;
            //h1_3=4.*sqrt(ae);

            //hb=_zpi/awk*0.142*tanh(d*awk);
            tmpd1 = _zpi;
            tmpd2 = awk;
            tmpd3 = d;
            tmpd4 = tmpd1 / tmpd2 * 0.142 * tanh(tmpd2 * tmpd3);
            hb = tmpd4;
            hbb=0.78125*d/1.6726;
            if(hb>=hbb) hb=hbb;

            if(h1_3>=hb)
            {
                chbh=(hb/h1_3)*(hb/h1_3);
                //if(my_rank == 0)
                //    printf("    @@@@    ia = %d, ic = %d, chbh = %.8f\n", ia, ic, chbh);
                index_tmp = 0;
                for(j = 1; j <= _jl; j++)
                {
                    for(k = 1; k <= _kl; k++)
                    {
                        //???
                        cpe_mean1_em[index_tmp]=max(chbh*cpe_mean1_em[index_tmp], _small);
                        index_tmp++;
                        //if(my_rank == 0)
                        //    printf("    @@@@  %d  %d  %.8f  \n", j, k, _em[index_tmp]);
                    }
                }
                h1_3=hb;
            }

            tpf = min(1.41*ape, tpf);
            /*
                if(my_rank == 0 && ia == 2 && ic == 2)
                    printf("    @@@@    %.8f, %.8f, %.8f, %.8f, %.8f, %.8f\n",
                           _ape[iaic_index], _tpf[iaic_index],
                           awk, _h1_3[iaic_index], hb, hbb);*/

            if(h1_3 <= 0.05)
            {
                h1_3=0.05;
            }

            cpe_mean1_data_out[APE_OFFSET] = ape;
            cpe_mean1_data_out[H1_3_OFFSET] = h1_3;
            cpe_mean1_data_out[TPF_OFFSET] = tpf;
            cpe_mean1_data_out[AET_OFFSET] = aet;

            RPCC(edcc);
            out_flag[PF_MEAN1_2] += edcc - stcc;
            stcc = edcc;

        }

        //if(param.host_mpi_rank == 0 && my_core == 0)
        //        printf("    ####    %d, %d, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, offset = %d\n",
        //               cpe_mean1_ia, cpe_mean1_ic, asi, awf, aetc, aets, awk, ae, cpe_mean1_data_out[APE_OFFSET], compute_unit * GETEM_OUT_DATA);

        // ????????д????????
        while (get_reply_god != get_reply_target);

        seq_write_cpe_mean1_em_data_to_mpe(my_core, compute_unit, 1, mean1_buffer);

        standard_write_flag(my_core);

        RPCC(edcc);
        out_flag[PF_MEAN1_WRITE] += edcc - stcc;
        stcc = edcc;
    }
}

void cpe_athread_daemon(void *_param)
{
    cpe_mean1_next_dwf = &cpe_ips_wf[0];
    cpe_mean1_next_wf = &cpe_ips_ccg[0];

    int my_core, i;
    double tmpd;
    param = *((struct cpe_init_param*) _param);

    my_core = athread_get_id(-1);

    if(param.host_mpi_rank == 0 && my_core == 0)
        printf(" ", param.host_packed_ranges, param.host_ikp, param.host_ikp1, param.host_ikm, param.host_ikm1, param.host_jp1,
               param.host_jp2, param.host_jm1, param.host_jm2, param.host_packed_consts, param.host_wp, param.host_wm, param.host_wks17, param.host_wk,
               param.host_dwk, param.host_grolim, param.host_wkh, param.host_cos_theta, param.host_sin_theta, param.host_wf, param.host_ccg, param.host_e,
               param.host_ee, param.host_pexx_output_buf, param.host_flag, param.slave_flag, param.host_mpi_rank, param.host_ppg_ia_set, param.host_ppg_ic_set, param.host_ppg_inner_groups,
               param.host_ppg_inner_left, param.host_ppg_idxs, param.host_ppg_values, param.host_halo_ppg_ia_set, param.host_halo_ppg_ic_set, param.host_mean1_ia_set, param.host_mean1_ic_set, param.host_mean1_in_data,
               param.host_mean1_out_buf, param.host_dwf, param.host_glbflag, param.host_ix2, param.host_ssbos);

    wild_my_core = my_core;
    for (i = 0; i < FLAG_SIZE; i++){ local_flag[i] = 0; out_flag[i] = 0;}
    slave_flag_to_wait = 1;

    copy_ips_permanent_arrays(my_core);
    for (i = 0;i < _kldp1; i++)cpe_ips_sqrt_wk [i] = sqrt(cpe_ips_wk [i]); 
    for (i = 0;i < _kldp1; i++)
    {
        tmpd = cpe_ips_wk[i];
        tmpd = 1.0 / sqrt(tmpd);
        cpe_ips_sqrt_div_wk[i] = tmpd;
    }
    for(i = 0; i < _jl * _kl; i++)
        cpe_nomarine_es[i] = _small;

    int jk;
    for (jk = 0; jk < _jl * _kl; jk+=4)
    {
        for(i = 0; i < 4; i++)
        {
            cpe_ppg_jv[jk + i] = (jk + i) / _kl;
            cpe_ppg_kv[jk + i] = ((jk + i) % _kl)  * _jl;
        }
    }

    int index = 0;
    int j11, j12, j21, j22, jpjm_index, mr, k, j;
    
    for (mr = 1; mr <= 2; mr++)
    {
        for (j = 1; j <= _jl; j++)
        {
            jpjm_index = _jpjm_index(mr, j);
            j11 = (cpe_ips_jp1[jpjm_index] - 1);
            j12 = (cpe_ips_jp2[jpjm_index] - 1);
            j21 = (cpe_ips_jm1[jpjm_index] - 1);
            j22 = (cpe_ips_jm2[jpjm_index] - 1);
            cpe_ips_jpms[index] = j11;
            cpe_ips_jpms[index + 1] = j12;
            cpe_ips_jpms[index + 2] = j21;
            cpe_ips_jpms[index + 3] = j22;
            index+=4;      
        }
    }

/*
    int jse11, jse12;
    int jse21, jse22;
    int kp2, kp3, im, im1, kp, kp1;
    int jpms_index2 = 0, se_index;

    int touch_counter[_klp1 * _jl], counter_max;
    memset(&touch_counter[0], 0, sizeof(int) * _klp1 * _jl);

    for (k = 1; k <= _kl; k++)
    {
        im     = cpe_ips_ikm [k - 1];
        im1    = cpe_ips_ikm1[k - 1];
        kp     = cpe_ips_ikp [k - 1];
        kp1    = cpe_ips_ikp1[k - 1];
        kp2    = cpe_ips_ikp [k - 1];
        kp3    = cpe_ips_ikp1[k - 1];

        if (kp >= _kl)
        {
            kp2 = _kl + 1;
            if (kp == _kl) kp2 = _kl;
            kp  = _kl;
            kp1 = _kl;
            kp3 = _kl + 1;
        }

        kp = _jl * (kp - 1);
        kp1 = _jl * (kp1 - 1);
        kp2 = _jl * (kp2 - 1);
        kp3 = _jl * (kp3 - 1);
        im = _jl * (im - 1);
        im1 = _jl * (im1 - 1);

        jpms_index2 = 0;
        for (mr = 1; mr <= 2; mr++)
        {
            for (j = 1; j <= _jl; j++)
            {
                jse11 = cpe_ips_jpms[jpms_index2];
                jse12 = cpe_ips_jpms[jpms_index2 + 1];
                jse21 = cpe_ips_jpms[jpms_index2 + 2];
                jse22 = cpe_ips_jpms[jpms_index2 + 3];

                se_index = kp2 + jse11;
                touch_counter[se_index]++;
                se_index = kp2 + jse12;
                touch_counter[se_index]++;
                se_index = kp3 + jse11;
                touch_counter[se_index]++;
                se_index = kp3 + jse12;
                touch_counter[se_index]++;
                se_index = im + jse21;
                touch_counter[se_index]++;
                se_index = im + jse22;
                touch_counter[se_index]++;
                se_index = im1 + jse21;
                touch_counter[se_index]++;
                se_index = im1 + jse22;
                touch_counter[se_index]++;

                jpms_index2 += 4;
            }
        }
    }

    if(param.host_mpi_rank == 0 && my_core == 0)
    {
        counter_max = 0;
        for(i = 0; i < _klp1 * _jl; i++)
        {
            printf(" **%d** ", touch_counter[i]);
            counter_max = max(counter_max, touch_counter[i]);
        }
        printf("\n counter_max = %d\n", counter_max);
    }*/

    //if(my_core == 0)
    //    printf("    ####    cpe 0 of host mpi %d\n", param.host_mpi_rank);
    if(param.host_mpi_rank == 0 && my_core == 0)
        printf("    ####    ppg left = %d\n", param.host_ppg_inner_left);

    volatile long ppg_cc, ips_cc, halo_ppg_cc, halo_ips_cc, mean1_cc;
    volatile long ppg_time, ips_time, halo_ppg_time, halo_ips_time, mean1_time;
    double ppg_t, ips_t;

    ppg_cc = ips_cc = halo_ppg_cc = halo_ips_cc = mean1_cc = 0;
    ppg_time = ips_time = halo_ppg_time = halo_ips_time = mean1_time = 0;

    while (1)
    {
        standard_wait_flag(my_core);

        if (local_flag[KERNEL_ACTION] == PROPAGAT_FLAG)
        {
            RPCC(stcc);
            cpe_propagat(my_core);
            RPCC(edcc);
            ppg_cc += edcc - stcc;
            ppg_time++;
        }

        if (local_flag[KERNEL_ACTION] == MEAN1_FLAG)
        {
            RPCC(stcc);
            cpe_mean1(my_core);
            RPCC(edcc);
            mean1_cc += edcc - stcc;
            mean1_time++;
        }

        if (local_flag[KERNEL_ACTION] == EXIT_FLAG)
        {
            standard_write_flag(my_core);
            standard_write_flag(my_core);
            break;
        }
    }
}
