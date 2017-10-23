#ifndef C_PUBLIC_VAR_H
#define C_PUBLIC_VAR_H

#include "slave_kernel.h"
#include <stdio.h>

/* !!!!! 变量名以下划线开始的对应原 Fortran 代码里的同名常量 !!!!! */

extern float _deltth;
extern int _glbflag;

/* ========== 以下两个区的由于需要共用，所以要用 extern ========== */
/* ---------- 范围变量区 ---------- */
extern int   _ixs, _ixl, _iys, _iyl;
//extern int   _ix2;
extern int   dim4, dim3, dim2, ix_size, iy_size, ixiy_size;
extern float _deltt5, _deltt, _cong, _al31, _al21, _al13, _al23, _al11;

/* ---------- 数组指针区 ---------- */
extern int   *_ikp, *_ikp1, *_kpmt0, *_kakt0, *_ks0;
extern int   *_ikm, *_ikm1, *_jp1, *_jp2, *_jm1, *_jm2;
extern float *_wx, *_wy, *_wk, *_nsp, *_fconst0, *_awk, *_e, *_w, *_ee;
extern float *_ae, *_asi, *_awf, *_awk, *_ark, *_wf, *_dwk, *_wkh, *_dse;
extern float *_ssbos, *_grolim, *_pein, *_peds, *_pebo, *_sein, *_sebo, *_seds;
extern float *_uxx, *_uxy, *_uyx, *_uyy, *_thet, *_ccg, *_d, *_enh;
extern float *_wp, *_wm, *_wks17, *_x, *_y;
extern float *_ape, *_aet, *_hb, *_hbb, *_h1_3, *_dwf, *_tpf;

// propagat 中使用的缓存起来的坐标和中间变量
extern int   *_idxs;
extern float *_tmp_values;

extern int   *ppg_ia_set, *ppg_ic_set, *ppg_nsp_set;  // propagat inner 实际计算格点
//extern int   *ips_ia_set, *ips_ic_set;  // implsch inner  实际计算格点
extern int   *halo_ppg_ia_set, *halo_ppg_ic_set, *halo_ppg_nsp_set, halo_ppg_marine_count;
//extern int   *halo_ips_ia_set, *halo_ips_ic_set, halo_ips_marine_count;
extern int   halo_ppg_group_count, halo_ppg_group_remain;
//extern int   halo_ips_group_count, halo_ips_group_remain;
extern volatile int   halo_ppg_collect_flag, halo_ips_collect_flag;
//extern float* halo_ips_packed_data, *halo_ips_out_buffer;
extern float* halo_ppg_packed_data, *halo_ppg_out_buffer;
extern int   *mean1_ia_set, *mean1_ic_set, *mean1_nsp_set, mean1_marine_count;
extern int  mean1_group_count, mean1_group_remain;
extern float* mean1_packed_float, *mean1_out_buffer;
extern int*   mean1_packed_int;
extern float* mean1_n;
extern int mean1_halo_marine_count;
extern int ppg_inner2_marine_count;

extern float *cos_theta, *sin_theta;

//extern float *_em;
extern volatile long host_flag[FLAG_SIZE];
extern volatile long slave_flag[FLAG_SIZE];
extern volatile long flag_to_wait;

extern long ips_pack_total, ips_write_total, ips_remain_total, ips_call_total;
extern long mean1_wait_total, mean1_remain_total, mean1_copy_total, mean1_call_total;

extern int _ix2;

extern int ppg_inner_total_cu;

// 将原 fortran 程序里的值和数组指针拷贝到本地
void c_implement_init_once_(int   *ixs_, int   *ixl_,   int   *iys_, int   *iyl_,
    float *wx_,  float *wy_,    float *wk_,  float *nsp_,
    int *kpmt0_, int   *kakt0_, int   *ks0_, float *fconst0_,
    float *awk_, float *ae_,    float *asi_, float *awf_,
    float *ark_, float *wf_,    float *dwk_, float *wkh_,
    float *e_,   float *w_,     float *ee_,  float *dse_,
    float *ssbos_,  float *deltt_, float *deltt5_, float *grolim_,
    float *pein_, float *peds_, float *pebo_, float *sein_,
    float *sebo_, float *seds_, float *uxx_,  float *uxy_,
    float *uyx_,  float *uyy_,  float *thet_, float *ccg_,
    float *d_,    float *enh_,  float *wp_,   float *wm_,
    int *ikp_, int *ikp1_, int *ikm_, int *ikm1_, float *wks17_,
    int *jp1_, int *jp2_, int *jm1_, int *jm2_, float *cong_,
    float *al31_, float *al21_, float *al13_, float *al23_,
    float *al11_, float *x_, float *y_,
    float *ape_, float *aet_, float *hb_, float *hbb_, float *h1_3_, float *dwf_, float *tpf_,
    float *deltth_, int *glbflag_, int *ix2_);

void wait_slave_flag();

static inline unsigned long rpcc()
{
    unsigned long time;
    asm("rtc %0": "=r" (time) : );
    return time;
}


#include "c_public_const.h"

//Huang Chenghuan, 20170318
extern volatile int my_rank;

//void cpu_mean1_emulator();

#endif
