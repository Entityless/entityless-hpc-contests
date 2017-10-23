#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../c_header/c_public_var.h"
#include "../c_header/slave_kernel.h"
#include "../c_header/athread_emulator.h"


int   has_inited_implsch_param = 0;
int   *emu_cpe_ips_ikp, *emu_cpe_ips_ikp1, *emu_cpe_ips_ikm, *emu_cpe_ips_ikm1;
int   *emu_cpe_ips_jp1, *emu_cpe_ips_jp2, *emu_cpe_ips_jm1, *emu_cpe_ips_jm2;
int   *emu_cpe_ips_ia_set, *emu_cpe_ips_ic_set;
float *emu_cpe_ips_wp, *emu_cpe_ips_wm, *emu_cpe_ips_wks17, *emu_cpe_ips_wk;
float *emu_cpe_ips_dwk, *emu_cpe_ips_grolim, *emu_cpe_ips_wkh, *emu_cpe_ips_cos_theta, *emu_cpe_ips_sin_theta;
float *emu_cpe_ips_wf, *emu_cpe_ips_ccg;
float *emu_cpe_ips_packed_data, *emu_cpe_pexx_output_buf;
float *emu_cpe_ips_e, *emu_cpe_ips_ee;

struct cpe_init_param _cpe_param;

struct emu_c_implsch_local_tmp  // 每一个线程使用的先写后用的中间变量
{
	float asi, awk, ae, ark, awf, ks0, enh;
	float fconst0[_kl];
	float   se[_klp1 * _jl];
	float  dse[_klp1 * _jl];
	float sein[_kl   * _jl];
	float sebo[_kl   * _jl];
	float seds[_kl   * _jl];
}emu_thread_local_buff;

// for debug, MPE 模拟 CPE 进行双线性插值计算
float emu_bilinear_interpolation_qr(float aa, float bb, float cc, float dd, float q, float r)
{
    float phia = (1 - q) * (1 - r); 
    float phib = q * (1 - r); 
    float phic = (1 - q) * r;
    float phid = q * r;
    return (aa * phia + bb * phib + cc * phic + dd * phid);
}

// for debug, MPE 模拟 CPE 计算一组 propagat
void emu_cpe_propagat_kernel(float *cpe_ppg_packed_data, float *cpe_e)
{
	int j, k, cnt = 0;
	float *packed_data_ptr;
	float es[4], exxyy;
	
	packed_data_ptr = cpe_ppg_packed_data;
	for (j = 1; j <= _jl; j++)
	{
		for (k = 1; k <= _kl; k++)
		{
			packed_data_ptr[E10_OFFSET] *= packed_data_ptr[FIEN_OFFSET];
			packed_data_ptr[E20_OFFSET] *= packed_data_ptr[FIEN1_OFFSET];
			packed_data_ptr[E30_OFFSET] *= packed_data_ptr[FIEN_OFFSET];
			packed_data_ptr[E40_OFFSET] *= packed_data_ptr[FIEN1_OFFSET];
			
			packed_data_ptr[E11_OFFSET] *= packed_data_ptr[FIEN_OFFSET];
			packed_data_ptr[E21_OFFSET] *= packed_data_ptr[FIEN1_OFFSET];
			packed_data_ptr[E31_OFFSET] *= packed_data_ptr[FIEN_OFFSET];
			packed_data_ptr[E41_OFFSET] *= packed_data_ptr[FIEN1_OFFSET];
			
			packed_data_ptr[E12_OFFSET] *= packed_data_ptr[FIEN_OFFSET];
			packed_data_ptr[E22_OFFSET] *= packed_data_ptr[FIEN1_OFFSET];
			packed_data_ptr[E32_OFFSET] *= packed_data_ptr[FIEN_OFFSET];
			packed_data_ptr[E42_OFFSET] *= packed_data_ptr[FIEN1_OFFSET];
			
			packed_data_ptr[E13_OFFSET] *= packed_data_ptr[FIEN_OFFSET];
			packed_data_ptr[E23_OFFSET] *= packed_data_ptr[FIEN1_OFFSET];
			packed_data_ptr[E33_OFFSET] *= packed_data_ptr[FIEN_OFFSET];
			packed_data_ptr[E43_OFFSET] *= packed_data_ptr[FIEN1_OFFSET];
			
			es[0] = emu_bilinear_interpolation_qr(
				packed_data_ptr[E10_OFFSET], packed_data_ptr[E20_OFFSET], packed_data_ptr[E30_OFFSET],
				packed_data_ptr[E40_OFFSET], packed_data_ptr[Q_4_OFFSET], packed_data_ptr[R_4_OFFSET]
			);
			
			es[1] = emu_bilinear_interpolation_qr(
				packed_data_ptr[E11_OFFSET], packed_data_ptr[E21_OFFSET], packed_data_ptr[E31_OFFSET],
				packed_data_ptr[E41_OFFSET], packed_data_ptr[Q_4_OFFSET], packed_data_ptr[R_4_OFFSET]
			);
			
			es[2] = emu_bilinear_interpolation_qr(
				packed_data_ptr[E12_OFFSET], packed_data_ptr[E22_OFFSET], packed_data_ptr[E32_OFFSET],
				packed_data_ptr[E42_OFFSET], packed_data_ptr[Q_4_OFFSET], packed_data_ptr[R_4_OFFSET]
			);
			
			es[3] = emu_bilinear_interpolation_qr(
				packed_data_ptr[E13_OFFSET], packed_data_ptr[E23_OFFSET], packed_data_ptr[E33_OFFSET],
				packed_data_ptr[E43_OFFSET], packed_data_ptr[Q_4_OFFSET], packed_data_ptr[R_4_OFFSET]
			);
			
			exxyy = emu_bilinear_interpolation_qr(
				es[0], es[1], es[2], es[3],
				packed_data_ptr[Q_OFFSET], packed_data_ptr[R_OFFSET]
			);
			
			cpe_e[cnt++] = max(exxyy, _small);
			
			packed_data_ptr += PPG_BI_DATA;
		}
	}
}

// 很奇怪，sunway 上这里传 cpe_param 的话 host_ips_{ia,ic}_set 这两个指针的值会变
// 传 *cpe_param 就没有问题
void emp_cpe_implsch_init_param(struct cpe_init_param *cpe_param) 
{
	if (has_inited_implsch_param > 0) return;
	has_inited_implsch_param++;
	
	// 在从核上应该用拷贝，这里用指针复制来代替
	emu_cpe_ips_ikp  = cpe_param->host_ikp;
	emu_cpe_ips_ikp1 = cpe_param->host_ikp1;
	emu_cpe_ips_ikm  = cpe_param->host_ikm;
	emu_cpe_ips_ikm1 = cpe_param->host_ikm1;
	emu_cpe_ips_jp1  = cpe_param->host_jp1;
	emu_cpe_ips_jp2  = cpe_param->host_jp2;
	emu_cpe_ips_jm1  = cpe_param->host_jm1;
	emu_cpe_ips_jm2  = cpe_param->host_jm2;
	emu_cpe_ips_wp   = cpe_param->host_wp;
	emu_cpe_ips_wm   = cpe_param->host_wm;
	emu_cpe_ips_wk   = cpe_param->host_wk;
	emu_cpe_ips_dwk  = cpe_param->host_dwk;
	emu_cpe_ips_wkh  = cpe_param->host_wkh;
	emu_cpe_ips_wf   = cpe_param->host_wf;
	emu_cpe_ips_ccg  = cpe_param->host_ccg;
	emu_cpe_ips_e    = cpe_param->host_e; 
	emu_cpe_ips_ee   = cpe_param->host_ee;
	emu_cpe_ips_wks17     = cpe_param->host_wks17;
	emu_cpe_ips_grolim    = cpe_param->host_grolim;
	emu_cpe_ips_cos_theta = cpe_param->host_cos_theta;
	emu_cpe_ips_sin_theta = cpe_param->host_sin_theta;
	emu_cpe_ips_ia_set    = cpe_param->host_ips_ia_set;
	emu_cpe_ips_ic_set    = cpe_param->host_ips_ic_set;
	emu_cpe_ips_packed_data = cpe_param->host_ips_packed_data;
	emu_cpe_pexx_output_buf = cpe_param->host_pexx_output_buf;
	
	_cpe_param = *cpe_param;
}

void emu_cpe_c_mean2(int core_id) // passed
{
	int j, k, k1, index_tmp;
	float dwkk, wfk, wfk1, wsk, wsk1, wkk, wkk1, ekj, ekj1;

	emu_thread_local_buff.ae  = 0.0;
	emu_thread_local_buff.asi = 0.0;
	emu_thread_local_buff.awf = 0.0;
	emu_thread_local_buff.awk = 0.0;
	emu_thread_local_buff.ark = 0.0;
	
	for (k = 1; k <= _kld; k++)
	{
		k1   = k + 1;
		dwkk = emu_cpe_ips_dwk[k - 1];
		wfk  = emu_cpe_ips_wf[k - 1];
		wfk1 = emu_cpe_ips_wf[k1 - 1];
		wsk  = _zpi * wfk;
		wsk1 = _zpi * wfk1;
		wkk  = emu_cpe_ips_wk[k - 1];
		wkk1 = emu_cpe_ips_wk[k1 - 1];
		for (j = 1; j <= _jl; j++)
		{
			index_tmp = (j - 1) * _kl;
			if (k < _kl)
			{
				ekj  = emu_cpe_ips_e[index_tmp + k - 1];
				ekj1 = emu_cpe_ips_e[index_tmp + k1 - 1];
			} else {
				ekj  = emu_cpe_ips_e[index_tmp + _kl - 1] * emu_cpe_ips_wkh[k - _kl + 1 - 1];
				ekj1 = emu_cpe_ips_e[index_tmp + _kl - 1] * emu_cpe_ips_wkh[k - _kl + 2 - 1];
			}
			emu_thread_local_buff.ae  =  emu_thread_local_buff.ae + (ekj + ekj1) * dwkk;
			emu_thread_local_buff.awf = emu_thread_local_buff.awf + (ekj * wfk + ekj1 * wfk1) * dwkk;
			emu_thread_local_buff.asi = emu_thread_local_buff.asi + (ekj / wsk + ekj1 / wsk1) * dwkk;
			emu_thread_local_buff.awk = emu_thread_local_buff.awk + (ekj * wkk + ekj1 * wkk1) * dwkk;
			emu_thread_local_buff.ark = emu_thread_local_buff.ark + (ekj / sqrt(wkk) + ekj1 / sqrt(wkk1)) * dwkk;
		}
	}
	emu_thread_local_buff.asi =  emu_thread_local_buff.ae / emu_thread_local_buff.asi;
	emu_thread_local_buff.awf = emu_thread_local_buff.awf /  emu_thread_local_buff.ae;
	emu_thread_local_buff.awk = emu_thread_local_buff.awk /  emu_thread_local_buff.ae;
	emu_thread_local_buff.ark =  emu_thread_local_buff.ae / emu_thread_local_buff.ark;
	emu_thread_local_buff.ark = emu_thread_local_buff.ark * emu_thread_local_buff.ark;
}

void emu_cpe_c_get_thresholds() //passed
{
	int i, kpmt, kakt, ks1, ks;
	float vx, vy, ww, wkpmt, wakt;

	vx = emu_cpe_ips_packed_data[WX_OFFSET];
	vy = emu_cpe_ips_packed_data[WY_OFFSET];
	ww = vx * vx + vy * vy;
	if (ww == 0.0) ww = 0.5;
	wkpmt = _cksp * _gc2 / ww;
	kpmt = (int) (log10(wkpmt / emu_cpe_ips_wk[0]) / _alog10pwk);
	kpmt += 1;
	
	wakt = _cksa * emu_thread_local_buff.awk;
	kakt = (int) (log10(wakt / emu_cpe_ips_wk[0]) / _alog10pwk);
	kakt += 1;
	
	ks1 = max(kpmt, kakt);
	ks  = min(ks1, _kl);

	emu_thread_local_buff.ks0   = ks;
	
	for (i = 1; i <= ks; i++) emu_thread_local_buff.fconst0[i - 1] = 1.0;
	for (i = ks + 1; i <= _kl; i++) emu_thread_local_buff.fconst0[i - 1] = 0.0;
}

void emu_cpe_c_snonlin() // passed
{
	float xx, ffacp, ffacp1, cwks17, up, up1, um, um1, sap, sam;
	float wp11, wp12, wp21, wp22, wm11, wm12, wm21, wm22;
	float wp112, wp122, wp212, wp222, wm112, wm122, wm212, wm222;
	float ea1, ea2, ea3, ea4, ea5, ea6, ea7, ea8;
	float eij, eij2, zua, ead1, ead2, fcen, ad, adp, adm, delad, deladp, deladm;
	int kh, k, mr, j, ip, ip1, im, im1, kp, kp1, kp2, kp3;
	int j11, j12, j21, j22, jpjm_index, se_index;

	xx = 0.75 * emu_cpe_ips_packed_data[D_OFFSET] * emu_thread_local_buff.awk;
	if (xx < 0.5) xx = 0.5;
	emu_thread_local_buff.enh = 1.0 + (5.5 / xx) * (1.0 - 0.833 * xx) * exp(-1.25 * xx);
	kh = 0;
	
	for (k = 1; k <= _kl; k++)
	{
		wp11 = emu_cpe_ips_wp[_wpwm_index(k, 1, 1)];
		wp12 = emu_cpe_ips_wp[_wpwm_index(k, 1, 2)];
		wp21 = emu_cpe_ips_wp[_wpwm_index(k, 2, 1)];
		wp22 = emu_cpe_ips_wp[_wpwm_index(k, 2, 2)];
		wm11 = emu_cpe_ips_wm[_wpwm_index(k, 1, 1)];
		wm12 = emu_cpe_ips_wm[_wpwm_index(k, 1, 2)];
		wm21 = emu_cpe_ips_wm[_wpwm_index(k, 2, 1)];
		wm22 = emu_cpe_ips_wm[_wpwm_index(k, 2, 2)];
		
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
		ip     = emu_cpe_ips_ikp [k - 1];
		ip1    = emu_cpe_ips_ikp1[k - 1];
		im     = emu_cpe_ips_ikm [k - 1];
		im1    = emu_cpe_ips_ikm1[k - 1];
		kp     = ip;
		kp1    = ip1;
		kp2    = ip;
		kp3    = ip1;
		cwks17 = _cong * emu_cpe_ips_wks17[k - 1];
		
		if (kp >= _kl)
		{
			kh  = kh + 1;
			kp2 = _kl + 1;
			if (kp == _kl) kp2 = _kl;
			kp  = _kl;
			kp1 = _kl;
			kp3 = _kl + 1;
			ffacp  = emu_cpe_ips_wkh[kh - 1];
			ffacp1 = emu_cpe_ips_wkh[kh + 1 - 1];
		}
		
		fcen = emu_thread_local_buff.fconst0[k - 1] * emu_thread_local_buff.enh;
		
		for (mr = 1; mr <= 2; mr++)
		{
			for (j = 1; j <= _jl; j++)
			{
				eij = emu_cpe_ips_e[(j - 1) * _kl + (k - 1)];
				if (eij < 1e-20) continue;
				
				jpjm_index = _jpjm_index(mr, j);
				j11 = emu_cpe_ips_jp1[jpjm_index];
				j12 = emu_cpe_ips_jp2[jpjm_index];
				j21 = emu_cpe_ips_jm1[jpjm_index];
				j22 = emu_cpe_ips_jm2[jpjm_index];
				
				ea1 = emu_cpe_ips_e[(j11 - 1) * _kl + (kp - 1)];
				ea2 = emu_cpe_ips_e[(j12 - 1) * _kl + (kp - 1)];
				ea3 = emu_cpe_ips_e[(j11 - 1) * _kl + (kp1 - 1)];
				ea4 = emu_cpe_ips_e[(j12 - 1) * _kl + (kp1 - 1)];
				ea5 = emu_cpe_ips_e[(j21 - 1) * _kl + (im - 1)];
				ea6 = emu_cpe_ips_e[(j22 - 1) * _kl + (im - 1)];
				ea7 = emu_cpe_ips_e[(j21 - 1) * _kl + (im1 - 1)];
				ea8 = emu_cpe_ips_e[(j22 - 1) * _kl + (im1 - 1)];
				
				up  = (wp11 * ea1 + wp12 * ea2) * ffacp;
				up1 = (wp21 * ea3 + wp22 * ea4) * ffacp1;
				um  =  wm11 * ea5 + wm12 * ea6;
				um1 =  wm21 * ea7 + wm22 * ea8;
				sap = up + up1;
				sam = um + um1;
				
				eij2 = eij * eij;
				zua  = 2.0 * eij / _al31;
				ead1 = sap / _al11 + sam / _al21;
				ead2 = -2.0 * sap * sam / _al31;
				
				ad   = cwks17 * (eij2 * ead1 + ead2 * eij) * fcen;
				adp  = ad / _al13;
				adm  = ad / _al23;
				
				delad  = cwks17 * (eij * 2.0 * ead1 + ead2) * fcen;
				deladp = cwks17 * (eij2 / _al11 - zua * sam) * fcen / _al13;
				deladm = cwks17 * (eij2 / _al21 - zua * sap) * fcen / _al23;
				
				se_index = _se_index(k, j);
				emu_thread_local_buff.se[se_index]  -= 2.0 * ad;
				emu_thread_local_buff.dse[se_index] -= 2.0 * delad;
				
				se_index = _se_index(kp2, j11);
				emu_thread_local_buff.se[se_index]  += adp * wp11;
				emu_thread_local_buff.dse[se_index] += deladp * wp112;
				
				se_index = _se_index(kp2, j12);
				emu_thread_local_buff.se[se_index]  += adp * wp12;
				emu_thread_local_buff.dse[se_index] += deladp * wp122;
				
				se_index = _se_index(kp3, j11);
				emu_thread_local_buff.se[se_index]  += adp * wp21;
				emu_thread_local_buff.dse[se_index] += deladp * wp212;
				
				se_index = _se_index(kp3, j12);
				emu_thread_local_buff.se[se_index]  += adp * wp22;
				emu_thread_local_buff.dse[se_index] += deladp * wp222;
				
				se_index = _se_index(im, j21);
				emu_thread_local_buff.se[se_index]  += adm * wm11;
				emu_thread_local_buff.dse[se_index] += deladm * wm112;
				
				se_index = _se_index(im, j22);
				emu_thread_local_buff.se[se_index]  += adm * wm12;
				emu_thread_local_buff.dse[se_index] += deladm * wm122;
				
				se_index = _se_index(im1, j21);
				emu_thread_local_buff.se[se_index]  += adm * wm21;
				emu_thread_local_buff.dse[se_index] += deladm * wm212;
				
				se_index = _se_index(im1, j22);
				emu_thread_local_buff.se[se_index]  += adm * wm22;
				emu_thread_local_buff.dse[se_index] += deladm * wm222;
			}
		}
	}
}

void emu_cpe_c_sinput() // passed
{
	int j, k, ks, sein_index, se_index, e_index;
	float vx, vy, cd, costh, sinth, wl, wlstar, wk0, wf0, ws0, bett, beta;
	
	ks = emu_thread_local_buff.ks0;
	vx = emu_cpe_ips_packed_data[WX_OFFSET];
	vy = emu_cpe_ips_packed_data[WY_OFFSET];
	cd = (0.8 + 0.065 * emu_cpe_ips_packed_data[W_OFFSET]) * 0.001;
	
	for (j = 1; j <= _jl; j++)
	{
		sinth  = emu_cpe_ips_sin_theta[j - 1]; 
		costh  = emu_cpe_ips_cos_theta[j - 1]; 
		wl     = vx * costh + vy * sinth;
		wlstar = wl * sqrt(cd);
		
		sein_index = _sexx_index(0, j);
		se_index = _se_index(0, j);
		e_index  = (j - 1) * _kl - 1;
		
		for (k = 1; k <= ks; k++)
		{
			sein_index++;
			se_index++;
			e_index++;
			wk0 = emu_cpe_ips_wk[k - 1];
			wf0 = emu_cpe_ips_wf[k - 1];
			ws0 = _zpi * wf0;
			bett = _beta10;
			beta = bett * (wk0 * 28.0 * wlstar - ws0);
			beta = max(0.0, beta);
			emu_thread_local_buff.sein[sein_index] = beta * emu_cpe_ips_e[e_index];
			emu_thread_local_buff.se[se_index] += emu_thread_local_buff.sein[sein_index];
			emu_thread_local_buff.dse[se_index] += beta;
		}
	}
}

void emu_cpe_c_sdissip() // passed
{
	int j, k, ks, se_index, seds_index, e_index;
	float arkss, ekspm, sds, ssds;
	
	ks = emu_thread_local_buff.ks0;
	arkss = emu_thread_local_buff.ark;
	ekspm = emu_thread_local_buff.ae * arkss * arkss / 0.0030162;
	sds = _d1 * emu_thread_local_buff.asi / arkss * sqrt(ekspm) * exp(-_d2 * 0.64 / ekspm);
	
	for (j = 1; j <= _jl; j++)
	{
		se_index = _se_index(0, j);
		seds_index = _sexx_index(0, j);
		e_index = (j - 1) * _kl - 1;
		for (k = 1; k <= ks; k++)
		{
			se_index++;
			seds_index++;
			e_index++;
			ssds = -_ads * sds * emu_cpe_ips_wk[k - 1];
			emu_thread_local_buff.seds[seds_index] = ssds * emu_cpe_ips_e[e_index];
			emu_thread_local_buff.se[se_index]  += emu_thread_local_buff.seds[seds_index];
			emu_thread_local_buff.dse[se_index] += ssds;
		}
	}
}

void emu_cpe_c_sbottom() // passed
{
	int j, k, ks, se_index, sebo_index;
	float sbo, d0, wk0, dk, ssbo;

	ks = emu_thread_local_buff.ks0;
	d0 = emu_cpe_ips_packed_data[D_OFFSET];
	sbo = 0.038 * 2.0 / _g;
	
	for (k = 1; k <= ks; k++)
	{
		wk0 = emu_cpe_ips_wk[k - 1];
		dk  = d0 * wk0;
		if (dk >= 30.0)
		{
			ssbo = 0.0;
		} else {
			ssbo = -_abo * sbo * wk0 / sinh(2.0 * dk);
		}
		for (j = 1; j <= _jl; j++)
		{
			se_index = _se_index(k, j);
			sebo_index = _sexx_index(k, j);
			emu_thread_local_buff.sebo[sebo_index] = ssbo * emu_cpe_ips_e[(j - 1) * _kl + (k - 1)];
			emu_thread_local_buff.se[se_index]  += emu_thread_local_buff.sebo[sebo_index];
			emu_thread_local_buff.dse[se_index] += ssbo;
		}
	}
}

void emu_cpe_c_scurrent() // passed
{
	int j, k, ks, se_index;
	float duxdx0, duxdy0, duydx0, duydy0, costh, sinth;
	float cost2, sint2, ws0, cgdc, cu1, cu2, cu3, sscu;

	ks     = emu_thread_local_buff.ks0;
	duxdx0 = emu_cpe_ips_packed_data[UXX_OFFSET];
	duxdy0 = emu_cpe_ips_packed_data[UXY_OFFSET];
	duydx0 = emu_cpe_ips_packed_data[UYX_OFFSET];
	duydy0 = emu_cpe_ips_packed_data[UYY_OFFSET];
	for (j = 1; j <= _jl; j++)
	{
		sinth = emu_cpe_ips_sin_theta[j - 1]; 
		costh = emu_cpe_ips_cos_theta[j - 1]; 
		cost2 = costh * costh;
		sint2 = sinth * sinth;
		se_index = _se_index(0, j);
		for (k = 1; k <= ks; k++)
		{
			se_index++;
			ws0 = _zpi * emu_cpe_ips_wf[k - 1];
			cgdc = _wk[k - 1] * emu_cpe_ips_ccg[k - 1] / ws0;
			cu1 = (cgdc * ( 1.0 + cost2) - 0.5) * duxdx0;
			cu2 = cgdc * sinth * costh * (duxdy0 + duydx0);
			cu3 = (cgdc * (1.0 + sint2) - 0.5) * duydy0;
			sscu = -_acu * (cu1 + cu2 + cu3);
			emu_thread_local_buff.se[se_index]  += sscu * emu_cpe_ips_e[(j - 1) * _kl + (k - 1)];
			emu_thread_local_buff.dse[se_index] += sscu;
		}
	}
}

void emu_cpe_c_update_ee_pe() // passed
{
	int i, j, k, ks, jk_index, e_index, se_index;
	float pein_tmp, peds_tmp,  pebo_tmp, cd, wstar, gadiag, eef, eefab, sig, deltee, tmp;

	ks = emu_thread_local_buff.ks0;
	cd = (0.8 + 0.065 * emu_cpe_ips_packed_data[W_OFFSET]) * 0.001;
	wstar = emu_cpe_ips_packed_data[W_OFFSET] * sqrt(cd);
	
	for (j = 1; j <= _jl; j++)
	{
		e_index  = (j - 1) * _kl - 1;
		se_index = _se_index(0, j);
		for (k = 1; k <= ks; k++)
		{
			e_index++;
			se_index++;
			gadiag = 1.0 - _deltt5 * emu_thread_local_buff.dse[se_index];
			gadiag = max(1.0, gadiag);
			eef = _deltt * emu_thread_local_buff.se[se_index] / gadiag;
			eefab = fabs(eef);
			eefab = max(eefab, 1e-19);
			sig = eef / eefab;
			deltee = wstar * emu_cpe_ips_grolim[k - 1];
			eefab = min(eefab, deltee);
			eef = emu_cpe_ips_e[e_index] + sig * eefab;
			tmp = max(eef, 0.0);
			emu_cpe_ips_ee[e_index] = max(_small, tmp);
		}
		
		e_index = (j - 1) * _kl;
		for (k = ks + 1; k <= _kl; k++)
		{
			i = k - ks + 1;
			tmp = emu_cpe_ips_ee[e_index + ks - 1] * emu_cpe_ips_wkh[i - 1];
			emu_cpe_ips_ee[e_index + k - 1] = max(_small, tmp);
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
			pein_tmp += 2.0 * _dwk[k - 1] * emu_thread_local_buff.sein[jk_index];
			pebo_tmp += 2.0 * _dwk[k - 1] * emu_thread_local_buff.sebo[jk_index];
			peds_tmp += 2.0 * _dwk[k - 1] * emu_thread_local_buff.seds[jk_index];
		}
	}
	
	emu_cpe_pexx_output_buf[PEIN_OFFSET] = pein_tmp;
	emu_cpe_pexx_output_buf[PEBO_OFFSET] = pebo_tmp;
	emu_cpe_pexx_output_buf[PEDS_OFFSET] = peds_tmp;
}

void emu_cpe_implsch_kernel(int ig, int core_id)
{
	int ia, ic, compute_unit, e_offset, wf_ccg_offset;
	
	compute_unit = ig * 64 + core_id;
	ia = emu_cpe_ips_ia_set[compute_unit];
	ic = emu_cpe_ips_ic_set[compute_unit];
	
	e_offset = _e_4d_index(1, 1, ia, ic);
	wf_ccg_offset = (ic - _iys) * (_kldp1 * ix_size) + (ia - _ixs) * _kldp1;
	
	// 每一轮使用 _e(:, :, ia, ic), _ee(:, :, ia, ic), _wf(:, ia, ic), _ccg(:, ia, ic)
	emu_cpe_ips_e    = _cpe_param.host_e   + e_offset;
	emu_cpe_ips_ee   = _cpe_param.host_ee  + e_offset;
	emu_cpe_ips_wf   = _cpe_param.host_wf  + wf_ccg_offset;
	emu_cpe_ips_ccg  = _cpe_param.host_ccg + wf_ccg_offset;
	emu_cpe_ips_packed_data = _cpe_param.host_ips_packed_data + compute_unit * IPS_SV_DATA;
	
	emu_cpe_c_mean2(core_id);
	emu_cpe_c_get_thresholds();
	
	memset(emu_thread_local_buff.se,  0, sizeof(float) * _jl * _klp1);
	memset(emu_thread_local_buff.dse, 0, sizeof(float) * _jl * _klp1);
	
	emu_cpe_c_snonlin ();
	emu_cpe_c_sinput  ();
	emu_cpe_c_sdissip ();
	emu_cpe_c_sbottom ();
	emu_cpe_c_scurrent();
	
	emu_cpe_c_update_ee_pe();
}