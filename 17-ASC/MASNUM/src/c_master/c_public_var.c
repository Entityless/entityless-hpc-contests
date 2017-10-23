#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "../c_header/c_public_var.h"

/* ========== 公用变量区 ========== */
float _deltth;
int _glbflag;

/* ---------- 范围变量区 ---------- */
int   _ixs, _ixl, _iys, _iyl;
float _deltt5, _deltt, _cong, _al31, _al21, _al13, _al23, _al11;
int dim4, dim3, dim2, ix_size, iy_size, ixiy_size;

/* ---------- 数组指针区 ---------- */
int   *_ikp, *_ikp1, *_kpmt0, *_kakt0, *_ks0;
int   *_ikm, *_ikm1, *_jp1, *_jp2, *_jm1, *_jm2;
float *_wx, *_wy, *_wk, *_nsp, *_fconst0, *_awk, *_e, *_w, *_ee;
float *_ae, *_asi, *_awf, *_awk, *_ark, *_wf, *_dwk, *_wkh, *_dse;
float *_ssbos, *_grolim, *_pein, *_peds, *_pebo, *_sein, *_sebo, *_seds;
float *_uxx, *_uxy, *_uyx, *_uyy, *_thet, *_ccg, *_d, *_enh;
float *_wp, *_wm, *_wks17, *_x, *_y;
float *_ape, *_aet, *_hb, *_hbb, *_h1_3, *_dwf, *_tpf;

// 缓存起来的坐标和中间变量
int   *_idxs;
float *_tmp_values;

int   *ppg_ia_set, *ppg_ic_set, *ppg_nsp_set;  // propagat inner 实际计算格点
//int   *ips_ia_set, *ips_ic_set;  // implsch inner  实际计算格点
int   *halo_ppg_ia_set, *halo_ppg_ic_set, *halo_ppg_nsp_set, halo_ppg_marine_count;
//int   *halo_ips_ia_set, *halo_ips_ic_set, halo_ips_marine_count;
int   halo_ppg_group_count, halo_ppg_group_remain;
//int   halo_ips_group_count, halo_ips_group_remain;
volatile int   halo_ppg_collect_flag, halo_ips_collect_flag;
//float* halo_ips_packed_data, *halo_ips_out_buffer;
float* halo_ppg_packed_data, *halo_ppg_out_buffer;
int   *mean1_ia_set, *mean1_ic_set, *mean1_nsp_set, mean1_marine_count;
int  mean1_group_count, mean1_group_remain;
float* mean1_packed_float, *mean1_out_buffer;
int* mean1_packed_int;
float* mean1_n;
int mean1_halo_marine_count;
int ppg_inner2_marine_count;

float *cos_theta, *sin_theta;

float* mean1_out_buffer;//no input buffer :)

/* ================================ */

int has_inited = 0;

//float *_em;
volatile long host_flag[FLAG_SIZE];
volatile long slave_flag[FLAG_SIZE];
volatile long flag_to_wait;
long ips_pack_total, ips_write_total, ips_remain_total, ips_call_total;
long mean1_wait_total, mean1_remain_total, mean1_copy_total, mean1_call_total;

int _ix2;

int ppg_inner_total_cu;

// 将原 fortran 程序里的值和数组指针拷贝到本地
void c_implement_init_once_(
	int   *ixs_, int   *ixl_,   int   *iys_, int   *iyl_,
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
    float *deltth_, int *glbflag_, int *ix2_
)
{
	int j;
	
	if (has_inited > 0) return;
	has_inited++;
	
	// 循环范围
	_ixs = (*ixs_); _ixl = (*ixl_);
	_iys = (*iys_); _iyl = (*iyl_);

    if(my_rank == -1)
        printf("", _ixs, _ixl, _iys, _iyl);
	
	ix_size   = _ixl - _ixs + 1;
	iy_size   = _iyl - _iys + 1;
	ixiy_size = ix_size * iy_size;
	dim4      = _kl * _jl * ix_size;
	dim3      = _kl * _jl;
	dim2      = _kl;
	
	// 原程序的全局数组
	_wx = wx_; _wy = wy_; _wk = wk_; _nsp = nsp_; _kpmt0 = kpmt0_; 
	_kakt0 = kakt0_; _ks0 = ks0_; _fconst0 = fconst0_; _awk = awk_;
	_ae = ae_; _asi = asi_; _awf = awf_; _awk = awk_; 
	_ark = ark_; _wf = wf_; _dwk = dwk_; _wkh = wkh_; 
	_e = e_; _w = w_; _ee = ee_; _dse = dse_;// _se = se_; 
	_grolim = grolim_; _pein = pein_; _peds = peds_; _pebo = pebo_; 
	_sein = sein_; _sebo = sebo_; _seds = seds_; _ccg = ccg_;
	_uxx = uxx_; _uxy = uxy_; _uyx = uyx_; _uyy = uyy_; _thet = thet_;
	_d = d_; _enh = enh_; _wp = wp_; _wm = wm_; _wks17 = wks17_;
	_ikp = ikp_; _ikp1 = ikp1_; _ikm = ikm_; _ikm1 = ikm1_;
	_jp1 = jp1_; _jp2 = jp2_; _jm1 = jm1_; _jm2 = jm2_; _x = x_; _y = y_;

    _ape = ape_, _aet = aet_, _hb = hb_, _hbb = hbb_, _h1_3 = h1_3_, _dwf = dwf_, _tpf = tpf_;
	
	_deltt  = (*deltt_);
	_deltt5 = (*deltt5_);
	_cong   = (*cong_);
	_al31 = (*al31_);
	_al21 = (*al21_);
	_al13 = (*al13_);
	_al23 = (*al23_);
	_al11 = (*al11_);
	
	ppg_ia_set = (int*) malloc(sizeof(int) * ixiy_size);
    ppg_ic_set = (int*) malloc(sizeof(int) * ixiy_size);
    ppg_nsp_set = (int*) malloc(sizeof(int) * ixiy_size);
    //ips_ia_set = (int*) malloc(sizeof(int) * ixiy_size);
    //ips_ic_set = (int*) malloc(sizeof(int) * ixiy_size);
    assert(ppg_ia_set != NULL && ppg_ic_set != NULL && ppg_nsp_set != NULL);
    //assert(ips_ia_set != NULL && ips_ic_set != NULL);

    halo_ppg_collect_flag = 0;
    halo_ppg_marine_count = 0;
    //halo_ips_collect_flag = 0;
    //halo_ips_marine_count = 0;
    halo_ppg_ia_set = (int*) malloc(sizeof(int) * ixiy_size);
    halo_ppg_ic_set = (int*) malloc(sizeof(int) * ixiy_size);
    halo_ppg_nsp_set = (int*) malloc(sizeof(int) * ixiy_size);
    //halo_ips_ia_set = (int*) malloc(sizeof(int) * ixiy_size);
    //halo_ips_ic_set = (int*) malloc(sizeof(int) * ixiy_size);
    mean1_ia_set = (int*) malloc(sizeof(int) * ixiy_size);
    mean1_ic_set = (int*) malloc(sizeof(int) * ixiy_size);
    mean1_nsp_set = (int*) malloc(sizeof(int) * ixiy_size);
    assert(halo_ppg_ia_set != NULL && halo_ppg_ic_set != NULL && halo_ppg_nsp_set != NULL);
    //assert(halo_ips_ia_set != NULL && halo_ips_ic_set != NULL);
    assert(mean1_ia_set != NULL && mean1_ic_set != NULL && mean1_nsp_set != NULL);
	
	cos_theta = (float*) malloc(sizeof(float) * _jlp1);
	sin_theta = (float*) malloc(sizeof(float) * _jlp1);
	assert(cos_theta != NULL && sin_theta != NULL);


	for (j = 0; j < _jlp1; j++)
	{
		cos_theta[j] = cos(_thet[j]);
		sin_theta[j] = sin(_thet[j]);
	}

    //Huang Chenghuan 20170319

    _deltth = (*deltth_);
    _glbflag = (*glbflag_);
    _ix2 = (*ix2_);

    _ssbos = ssbos_;

    if(my_rank == -1)
        printf(" ", _wx, _wy, _wk, _nsp, _kpmt0, _kakt0,
               _ks0, _fconst0, _awk, _ae, _asi, _awf, _awk,
               _ark, _wf, _dwk, _wkh, _e, _w, _ee, _dse, _grolim,
               _pein, _peds, _pebo, _sein, _sebo, _seds, _ccg,
               _uxx, _uxy, _uyx, _uyy, _thet, _d, _enh, _wp, _wm, _wks17,
               _ikp, _ikp1, _ikm, _ikm1, _jp1, _jp2, _jm1, _jm2, _x, _y,
               _ape, _aet, _hb, _hbb, _h1_3, _dwf, _tpf,
               _deltt, _deltt5, _cong, _al31, _al21, _al13, _al23, _al11,
               _deltth, _glbflag, _ix2, _ssbos
               );
}

//Huang Chenghuan, 20170318
volatile int my_rank;

void wait_slave_flag()
{
    while(slave_flag[0] < flag_to_wait);
    flag_to_wait++;
}
